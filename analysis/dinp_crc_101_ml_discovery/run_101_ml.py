from __future__ import annotations

from pathlib import Path
import json
import platform
import warnings

import numpy as np
import pandas as pd
import sklearn
from sklearn.base import clone
from sklearn.feature_selection import SelectKBest, f_classif, mutual_info_classif
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, confusion_matrix, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT
SEED = 20260909
PRIMARY_EXTERNAL_DATASETS = ["GSE10950", "GSE74602"]


def make_selectors(n_features: int) -> list[tuple[str, object, int | str]]:
    selectors: list[tuple[str, object, int | str]] = [
        ("all_features", SelectKBest(score_func=f_classif, k="all"), "all")
    ]
    for k in [10, 20, 40, 60, 80, 95]:
        if k <= n_features:
            selectors.append((f"f_classif_top{k}", SelectKBest(score_func=f_classif, k=k), k))
    for k in [10, 20, 40, 60]:
        if k <= n_features:
            selectors.append(
                (
                    f"mutual_info_top{k}",
                    SelectKBest(score_func=lambda X, y: mutual_info_classif(X, y, random_state=SEED), k=k),
                    k,
                )
            )
    return selectors


def make_classifiers() -> list[tuple[str, object]]:
    return [
        (
            "logistic_l2",
            LogisticRegression(C=1.0, penalty="l2", solver="lbfgs", class_weight="balanced", max_iter=5000, random_state=SEED),
        ),
        (
            "logistic_l1",
            LogisticRegression(C=0.5, penalty="l1", solver="liblinear", class_weight="balanced", max_iter=5000, random_state=SEED),
        ),
        (
            "linear_svm",
            SVC(C=1.0, kernel="linear", class_weight="balanced", probability=True, random_state=SEED),
        ),
        (
            "rbf_svm",
            SVC(C=1.0, kernel="rbf", gamma="scale", class_weight="balanced", probability=True, random_state=SEED),
        ),
        (
            "random_forest",
            RandomForestClassifier(
                n_estimators=300,
                max_features="sqrt",
                class_weight="balanced",
                random_state=SEED,
                n_jobs=1,
            ),
        ),
        (
            "extra_trees",
            ExtraTreesClassifier(
                n_estimators=300,
                max_features="sqrt",
                class_weight="balanced",
                random_state=SEED,
                n_jobs=1,
            ),
        ),
        (
            "hist_gradient_boosting",
            HistGradientBoostingClassifier(
                max_iter=200,
                learning_rate=0.05,
                max_leaf_nodes=15,
                l2_regularization=1.0,
                random_state=SEED,
            ),
        ),
        (
            "gradient_boosting",
            GradientBoostingClassifier(
                n_estimators=150,
                learning_rate=0.05,
                max_depth=2,
                random_state=SEED,
            ),
        ),
        (
            "knn",
            KNeighborsClassifier(n_neighbors=7, weights="distance"),
        ),
    ]


def make_model_catalog(n_features: int) -> list[dict]:
    selectors = make_selectors(n_features)
    classifiers = make_classifiers()
    catalog: list[dict] = []
    model_number = 1
    for selector_name, selector, requested_k in selectors:
        for classifier_name, classifier in classifiers:
            catalog.append(
                {
                    "model_id": f"M{model_number:03d}_{selector_name}__{classifier_name}",
                    "selector_name": selector_name,
                    "classifier_name": classifier_name,
                    "requested_features": requested_k,
                    "selective_selector": selector_name != "all_features",
                    "selector": clone(selector),
                    "classifier": clone(classifier),
                }
            )
            model_number += 1
    # Two independent full-feature baselines make the catalog exactly 101 models.
    for classifier_name, classifier in [
        ("lda_shrinkage", LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")),
        ("gaussian_nb", GaussianNB()),
    ]:
        catalog.append(
            {
                "model_id": f"M{model_number:03d}_all_features__{classifier_name}",
                "selector_name": "all_features",
                "classifier_name": classifier_name,
                "requested_features": "all",
                "selective_selector": False,
                "selector": SelectKBest(score_func=f_classif, k="all"),
                "classifier": clone(classifier),
            }
        )
        model_number += 1
    if len(catalog) != 101:
        raise AssertionError(f"Expected 101 models, got {len(catalog)}")
    return catalog


def load_cohorts() -> tuple[pd.DataFrame, dict[str, pd.DataFrame], list[str], list[str]]:
    overlap = pd.read_csv(OUTPUTS / "DINP_CRC_overlap.csv", dtype=str).fillna("")
    overlap["gene_symbol"] = overlap["gene_symbol"].str.strip().str.upper()
    overlap = overlap.drop_duplicates("gene_symbol")
    genes = overlap["gene_symbol"].tolist()

    manifest = pd.read_csv(OUTPUTS / "DINP_CRC_TCGA_GEO_sample_manifest.csv", dtype=str).fillna("")
    expression = pd.read_csv(OUTPUTS / "DINP_CRC_TCGA_GEO_target_expression_long.csv")
    expression["gene_symbol"] = expression["gene_symbol"].str.strip().str.upper()
    expression["sample_id"] = expression["sample_id"].astype(str)
    sample_manifest = manifest[manifest["dataset_id"].isin(["TCGA-COAD", "GSE10950", "GSE74602"])].copy()
    sample_ids = set(sample_manifest["sample_id"])
    expression = expression[expression["sample_id"].isin(sample_ids)].copy()
    duplicate_keys = expression.duplicated(["sample_id", "gene_symbol"]).sum()
    if duplicate_keys:
        raise ValueError(f"Expression matrix has {duplicate_keys} duplicate sample/gene keys")
    wide = expression.pivot(index="sample_id", columns="gene_symbol", values="expression_value")
    wide = wide.apply(pd.to_numeric, errors="coerce")
    measured_genes = [gene for gene in genes if gene in wide.columns]
    missing_genes = [gene for gene in genes if gene not in wide.columns]
    if len(measured_genes) < 2:
        raise ValueError("Too few measured DINP-CRC genes for ML")

    cohorts: dict[str, pd.DataFrame] = {}
    for dataset_id in ["TCGA-COAD", "GSE10950", "GSE74602"]:
        current = sample_manifest[sample_manifest["dataset_id"].eq(dataset_id)].copy()
        current = current[current["sample_id"].isin(wide.index)]
        current = current[current["group"].isin(["normal", "tumor"])]
        current = current.drop_duplicates("sample_id")
        current = current.set_index("sample_id")
        current["label"] = (current["group"] == "tumor").astype(int)
        data = wide.loc[current.index, measured_genes].copy()
        data.insert(0, "label", current["label"].astype(int))
        cohorts[dataset_id] = data

    return overlap, cohorts, measured_genes, missing_genes


def make_pipeline(selector: object, classifier: object) -> Pipeline:
    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("select", clone(selector)),
            ("model", clone(classifier)),
        ]
    )


def score_values(model: Pipeline, X: pd.DataFrame) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        return np.asarray(model.predict_proba(X)[:, 1], dtype=float)
    if hasattr(model, "decision_function"):
        return np.asarray(model.decision_function(X), dtype=float)
    return np.asarray(model.predict(X), dtype=float)


def classification_metrics(y_true: np.ndarray, scores: np.ndarray, predictions: np.ndarray) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=int)
    predictions = np.asarray(predictions, dtype=int)
    output: dict[str, float] = {}
    output["roc_auc"] = float(roc_auc_score(y_true, scores)) if len(np.unique(y_true)) == 2 else np.nan
    output["pr_auc"] = float(average_precision_score(y_true, scores)) if len(np.unique(y_true)) == 2 else np.nan
    output["balanced_accuracy"] = float(balanced_accuracy_score(y_true, predictions))
    output["f1"] = float(f1_score(y_true, predictions, zero_division=0))
    tn, fp, fn, tp = confusion_matrix(y_true, predictions, labels=[0, 1]).ravel()
    output["sensitivity"] = float(tp / (tp + fn)) if (tp + fn) else np.nan
    output["specificity"] = float(tn / (tn + fp)) if (tn + fp) else np.nan
    return output


def selected_features(model: Pipeline, feature_names: list[str]) -> list[str]:
    selector = model.named_steps["select"]
    support = selector.get_support()
    return [feature for feature, keep in zip(feature_names, support) if keep]


def safe_mean(values: list[float]) -> float:
    values = [float(value) for value in values if pd.notna(value)]
    return float(np.mean(values)) if values else np.nan


def main() -> None:
    overlap, cohorts, measured_genes, missing_genes = load_cohorts()
    train = cohorts["TCGA-COAD"]
    X_train = train[measured_genes]
    y_train = train["label"].to_numpy(dtype=int)
    if set(np.unique(y_train)) != {0, 1}:
        raise ValueError("TCGA-COAD training cohort must contain both normal and tumor labels")

    catalog = make_model_catalog(len(measured_genes))
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    model_rows: list[dict] = []
    fold_selection_counts = {gene: 0 for gene in measured_genes}
    selective_fold_selection_counts = {gene: 0 for gene in measured_genes}
    full_selection_counts = {gene: 0 for gene in measured_genes}
    selective_full_selection_counts = {gene: 0 for gene in measured_genes}
    gene_selected_models: dict[str, list[str]] = {gene: [] for gene in measured_genes}
    gene_external_metrics: dict[str, list[dict[str, float]]] = {gene: [] for gene in measured_genes}

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=UserWarning)
        for index, spec in enumerate(catalog, start=1):
            pipeline = make_pipeline(spec["selector"], spec["classifier"])
            fold_metrics: list[dict[str, float]] = []
            fold_selected: list[list[str]] = []
            for fold_index, (train_idx, valid_idx) in enumerate(cv.split(X_train, y_train), start=1):
                fold_model = clone(pipeline)
                fold_model.fit(X_train.iloc[train_idx], y_train[train_idx])
                fold_X = X_train.iloc[valid_idx]
                fold_y = y_train[valid_idx]
                fold_scores = score_values(fold_model, fold_X)
                fold_predictions = fold_model.predict(fold_X).astype(int)
                fold_metrics.append(classification_metrics(fold_y, fold_scores, fold_predictions))
                selected = selected_features(fold_model, measured_genes)
                fold_selected.append(selected)
                for gene in selected:
                    fold_selection_counts[gene] += 1
                    if spec["selective_selector"]:
                        selective_fold_selection_counts[gene] += 1

            full_model = clone(pipeline)
            full_model.fit(X_train, y_train)
            full_selected = selected_features(full_model, measured_genes)
            for gene in full_selected:
                full_selection_counts[gene] += 1
                gene_selected_models[gene].append(spec["model_id"])
                if spec["selective_selector"]:
                    selective_full_selection_counts[gene] += 1

            external_metrics: dict[str, dict[str, float]] = {}
            for dataset_id in PRIMARY_EXTERNAL_DATASETS:
                external = cohorts[dataset_id]
                X_external = external[measured_genes]
                y_external = external["label"].to_numpy(dtype=int)
                external_scores = score_values(full_model, X_external)
                external_predictions = full_model.predict(X_external).astype(int)
                metrics = classification_metrics(y_external, external_scores, external_predictions)
                external_metrics[dataset_id] = metrics
                for gene in full_selected:
                    gene_external_metrics[gene].append(metrics)

            row = {
                "model_id": spec["model_id"],
                "model_number": index,
                "selector_name": spec["selector_name"],
                "classifier_name": spec["classifier_name"],
                "requested_features": spec["requested_features"],
                "selective_selector": spec["selective_selector"],
                "full_train_selected_feature_count": len(full_selected),
                "full_train_selected_features": ";".join(full_selected),
                "internal_cv_folds": len(fold_metrics),
            }
            for metric in ["roc_auc", "pr_auc", "balanced_accuracy", "f1", "sensitivity", "specificity"]:
                values = [metrics[metric] for metrics in fold_metrics]
                row[f"internal_{metric}_mean"] = safe_mean(values)
                row[f"internal_{metric}_sd"] = float(np.nanstd(values, ddof=1)) if len(values) > 1 else np.nan
            for dataset_id, metrics in external_metrics.items():
                for metric, value in metrics.items():
                    row[f"{dataset_id}_{metric}"] = value
            row["external_mean_roc_auc"] = safe_mean(
                [external_metrics[dataset_id]["roc_auc"] for dataset_id in PRIMARY_EXTERNAL_DATASETS]
            )
            row["external_mean_pr_auc"] = safe_mean(
                [external_metrics[dataset_id]["pr_auc"] for dataset_id in PRIMARY_EXTERNAL_DATASETS]
            )
            row["external_mean_balanced_accuracy"] = safe_mean(
                [external_metrics[dataset_id]["balanced_accuracy"] for dataset_id in PRIMARY_EXTERNAL_DATASETS]
            )
            row["external_auc_ge_0_75_both_datasets"] = all(
                external_metrics[dataset_id]["roc_auc"] >= 0.75 for dataset_id in PRIMARY_EXTERNAL_DATASETS
            )
            model_rows.append(row)
            if index % 10 == 0 or index == len(catalog):
                print(f"completed {index}/101 models", flush=True)

    model_results = pd.DataFrame(model_rows)
    selective_count = int(model_results["selective_selector"].sum())
    all_count = int(len(model_results))
    if selective_count != 90:
        raise AssertionError(f"Expected 90 selective models, got {selective_count}")

    # External-validation-qualified models are used only for the second stability
    # component. The qualification threshold is fixed before comparing against
    # the existing PPI/pathway/Tier 1 ranking.
    model_results["external_validation_qualified"] = (
        model_results["selective_selector"]
        & model_results["GSE10950_roc_auc"].ge(0.75)
        & model_results["GSE74602_roc_auc"].ge(0.75)
    )
    qualified_models = model_results[model_results["external_validation_qualified"]].copy()
    qualified_model_count = int(len(qualified_models))
    if qualified_model_count == 0:
        raise AssertionError("No external-validation-qualified models")
    qualified_selection_counts = {gene: 0 for gene in measured_genes}
    qualified_external_metrics: dict[str, list[dict[str, float]]] = {gene: [] for gene in measured_genes}
    for _, model_row in qualified_models.iterrows():
        selected = [gene for gene in str(model_row["full_train_selected_features"]).split(";") if gene]
        for gene in selected:
            qualified_selection_counts[gene] += 1
            qualified_external_metrics[gene].append(
                {
                    "roc_auc": float(model_row["external_mean_roc_auc"]),
                    "pr_auc": float(model_row["external_mean_pr_auc"]),
                    "balanced_accuracy": float(model_row["external_mean_balanced_accuracy"]),
                }
            )

    gene_rows: list[dict] = []
    cross_rank_path = OUTPUTS / "DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv"
    cross_rank = pd.read_csv(cross_rank_path, dtype=str).fillna("") if cross_rank_path.exists() else pd.DataFrame()
    if not cross_rank.empty:
        cross_rank["gene_symbol"] = cross_rank["gene_symbol"].str.upper()
        cross_rank = cross_rank.drop_duplicates("gene_symbol")

    for gene in overlap["gene_symbol"]:
        measured = gene in measured_genes
        full_count = full_selection_counts.get(gene, 0)
        selective_full_count = selective_full_selection_counts.get(gene, 0)
        fold_count = fold_selection_counts.get(gene, 0)
        selective_fold_count = selective_fold_selection_counts.get(gene, 0)
        selected_model_metrics = gene_external_metrics.get(gene, [])
        selected_qualified_metrics = qualified_external_metrics.get(gene, [])
        mean_external_auc = safe_mean([metrics["roc_auc"] for metrics in selected_model_metrics])
        mean_external_pr_auc = safe_mean([metrics["pr_auc"] for metrics in selected_model_metrics])
        mean_external_bal_acc = safe_mean([metrics["balanced_accuracy"] for metrics in selected_model_metrics])
        mean_qualified_external_auc = safe_mean([metrics["roc_auc"] for metrics in selected_qualified_metrics])
        stability = selective_full_count / selective_count if measured else 0.0
        fold_stability = selective_fold_count / (selective_count * 5) if measured else 0.0
        qualified_stability = qualified_selection_counts.get(gene, 0) / qualified_model_count if measured else 0.0
        external_auc_score = float(np.clip(2.0 * (mean_external_auc - 0.5), 0.0, 1.0)) if pd.notna(mean_external_auc) else 0.0
        qualified_auc_score = float(np.clip(2.0 * (mean_qualified_external_auc - 0.5), 0.0, 1.0)) if pd.notna(mean_qualified_external_auc) else 0.0
        ml_priority = float(np.cbrt(stability * qualified_stability * qualified_auc_score)) if measured else 0.0
        row = {
            "gene_symbol": gene,
            "expression_measured": measured,
            "full_train_selection_count_all_101": full_count,
            "full_train_selection_frequency_all_101": full_count / all_count if measured else 0.0,
            "full_train_selection_count_selective_90": selective_full_count,
            "model_selection_frequency_selective_90": stability,
            "cv_fold_selection_count_all_101x5": fold_count,
            "cv_fold_selection_frequency_all_101x5": fold_count / (all_count * 5) if measured else 0.0,
            "cv_fold_selection_count_selective_90x5": selective_fold_count,
            "cv_fold_selection_frequency_selective_90x5": fold_stability,
            "selected_model_count_with_external_metrics": len(selected_model_metrics) // len(PRIMARY_EXTERNAL_DATASETS),
            "mean_external_roc_auc_when_selected": mean_external_auc,
            "mean_external_pr_auc_when_selected": mean_external_pr_auc,
            "mean_external_balanced_accuracy_when_selected": mean_external_bal_acc,
            "external_auc_score_0_to_1": external_auc_score,
            "external_qualified_model_count": qualified_selection_counts.get(gene, 0),
            "external_qualified_model_selection_frequency": qualified_stability,
            "mean_external_roc_auc_in_qualified_models": mean_qualified_external_auc,
            "external_qualified_auc_score_0_to_1": qualified_auc_score,
            "ml_priority_score": ml_priority,
            "selected_model_ids": ";".join(gene_selected_models.get(gene, [])),
            "stable_ml_flag": bool(measured and stability >= 0.80 and qualified_stability >= 0.80),
        }
        if not cross_rank.empty and gene in set(cross_rank["gene_symbol"]):
            cross_row = cross_rank[cross_rank["gene_symbol"].eq(gene)].iloc[0]
            for column in ["cross_rank", "candidate_tier", "cross_support_flag", "cross_rank_score", "hub_rank", "representative_term_count", "independent_datasets_significant_fdr_lt_0_05", "consensus_direction"]:
                row[f"{column}_cross_rank"] = cross_row.get(column, "")
        else:
            for column in ["cross_rank", "candidate_tier", "cross_support_flag", "cross_rank_score", "hub_rank", "representative_term_count", "independent_datasets_significant_fdr_lt_0_05", "consensus_direction"]:
                row[f"{column}_cross_rank"] = ""
        gene_rows.append(row)

    gene_stability = pd.DataFrame(gene_rows)
    gene_stability = gene_stability.sort_values(
        ["ml_priority_score", "model_selection_frequency_selective_90", "cv_fold_selection_frequency_selective_90x5", "gene_symbol"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    gene_stability.insert(0, "ml_rank", np.arange(1, len(gene_stability) + 1))

    model_results = model_results.sort_values(
        ["external_mean_roc_auc", "internal_roc_auc_mean", "external_mean_balanced_accuracy", "model_id"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    model_results.insert(0, "external_priority_rank", np.arange(1, len(model_results) + 1))

    model_results.to_csv(OUTPUTS / "DINP_CRC_101ML_model_results.csv", index=False)
    gene_stability.to_csv(OUTPUTS / "DINP_CRC_101ML_gene_stability.csv", index=False)
    gene_stability[gene_stability["stable_ml_flag"]].head(20).to_csv(
        OUTPUTS / "DINP_CRC_101ML_stable_genes_top20.csv", index=False
    )

    catalog_output = pd.DataFrame(
        [
            {
                "model_id": spec["model_id"],
                "model_number": index,
                "selector_name": spec["selector_name"],
                "classifier_name": spec["classifier_name"],
                "requested_features": spec["requested_features"],
                "selective_selector": spec["selective_selector"],
                "classifier_parameters": json.dumps(spec["classifier"].get_params(), sort_keys=True, default=str),
            }
            for index, spec in enumerate(catalog, start=1)
        ]
    )
    catalog_output.to_csv(OUTPUTS / "DINP_CRC_101ML_model_catalog.csv", index=False)

    stable = gene_stability[gene_stability["stable_ml_flag"]].copy()
    tier1 = gene_stability[gene_stability["cross_support_flag_cross_rank"].astype(str).str.lower().eq("true")]
    overlap_stable_tier1 = stable[stable["cross_support_flag_cross_rank"].astype(str).str.lower().eq("true")]
    top_for_report = gene_stability.head(20)
    summary_lines = [
        "# DINP–CRC 101-model machine-learning discovery",
        "",
        "## Frozen design",
        "",
        "- Input universe: all 97 DINP–CRC overlap genes; Tier 1/cross-ranking labels were not used to select features or fit models.",
        "- Training/discovery cohort: TCGA-COAD (tumor vs normal), stratified 5-fold cross-validation.",
        "- Primary external validation: GSE10950 and GSE74602. GSE156355 was not used in this discovery/validation run.",
        "- Models: 101 total = 11 feature-selection configurations × 9 classifiers + 2 full-feature baselines.",
        "- Feature selection and scaling were fitted inside each training fold; the final model was refit on all TCGA-COAD samples before external validation.",
        "",
        "## 101-model catalog",
        "",
        "Feature selectors: all features; ANOVA F-test top 10/20/40/60/80/95; mutual-information top 10/20/40/60.",
        "Classifiers: logistic L2, logistic L1, linear SVM, RBF SVM, random forest, extra-trees, histogram gradient boosting, gradient boosting, and distance-weighted kNN. Additional baselines: shrinkage LDA and Gaussian naive Bayes.",
        "",
        "## Stability and external-validation ranking",
        "",
        "`model_selection_frequency_selective_90` is the fraction of the 90 non-all-feature models whose full-TCGA fitted selector retained a gene. The CV-fold frequency is reported separately. External-validation-qualified models are selective models with ROC-AUC ≥0.75 in both GSE10950 and GSE74602. The ML priority score is the geometric mean of all-selective-model frequency, qualified-model frequency, and qualified-model external ROC-AUC score; it does not use Tier 1 labels.",
        "",
        f"- Frozen input genes: {len(overlap)}",
        f"- Measured expression features used for fitting: {len(measured_genes)}",
        f"- Input genes without expression values: {len(missing_genes)} ({', '.join(missing_genes)})",
        f"- Models completed: {len(model_results)}",
        f"- External-validation-qualified selective models (ROC-AUC ≥0.75 in both GEO cohorts): {qualified_model_count}",
        f"- Default stable-ML flag (selection frequency ≥80% in all selective models and ≥80% in qualified models): {len(stable)}",
        f"- Stable ML genes overlapping pre-existing Tier 1: {len(overlap_stable_tier1)}",
        "",
        "## Top 20 ML-priority genes",
        "",
        "| ML rank | Gene | Stable flag | Selective frequency | Qualified frequency | Qualified mean external AUC | Existing cross-rank | Tier 1 overlap | Direction |",
        "|---:|---|---|---:|---:|---:|---:|---|---|",
    ]
    for _, row in top_for_report.iterrows():
        external_auc = "" if pd.isna(row["mean_external_roc_auc_in_qualified_models"]) else f"{row['mean_external_roc_auc_in_qualified_models']:.3f}"
        cross_rank_value = row.get("cross_rank_cross_rank", "")
        tier1_flag = str(row.get("cross_support_flag_cross_rank", "")).lower() == "true"
        summary_lines.append(
            f"| {int(row['ml_rank'])} | {row['gene_symbol']} | {'yes' if row['stable_ml_flag'] else 'no'} | {row['model_selection_frequency_selective_90']:.3f} | {row['external_qualified_model_selection_frequency']:.3f} | {external_auc} | {cross_rank_value} | {'yes' if tier1_flag else 'no'} | {row.get('consensus_direction_cross_rank', '')} |"
        )
    summary_lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This stage is discovery plus external classification validation for the tumor/normal expression phenotype. It does not establish DINP causality, and the high separability of tumor versus normal should not be confused with exposure-response evidence.",
            "",
            "Tier 1/cross-rank overlap is an independent post hoc concordance check. It is not part of the ML feature-selection objective and should be reviewed before any MR instrument decision.",
        ]
    )
    (OUTPUTS / "DINP_CRC_101ML_summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    log = {
        "seed": SEED,
        "python": platform.python_version(),
        "scikit_learn": sklearn.__version__,
        "input_gene_count": int(len(overlap)),
        "measured_feature_count": int(len(measured_genes)),
        "unmeasured_input_genes": missing_genes,
        "training_dataset": "TCGA-COAD",
        "external_validation_datasets": PRIMARY_EXTERNAL_DATASETS,
        "excluded_from_this_run": ["GSE156355", "TCGA paired sensitivity", "TCGA-vs-GTEx sensitivity"],
        "training_sample_counts": train["label"].value_counts().sort_index().rename({0: "normal", 1: "tumor"}).to_dict(),
        "external_sample_counts": {
            dataset_id: cohorts[dataset_id]["label"].value_counts().sort_index().rename({0: "normal", 1: "tumor"}).to_dict()
            for dataset_id in PRIMARY_EXTERNAL_DATASETS
        },
        "model_count": int(len(model_results)),
        "selective_model_count": selective_count,
        "external_validation_qualified_model_count": qualified_model_count,
        "cv": {"method": "StratifiedKFold", "n_splits": 5, "shuffle": True, "random_state": SEED},
        "feature_selection_leakage_control": "selectors and StandardScaler fitted within each training fold",
        "external_validation_qualified_definition": "selective model with GSE10950 ROC-AUC >= 0.75 and GSE74602 ROC-AUC >= 0.75",
        "stable_ml_definition": "full-TCGA selection frequency >= 0.80 in all 90 selective models and >= 0.80 in external-validation-qualified selective models",
        "tier1_used_in_model_fit": False,
        "qc": {
            "model_ids_unique": int(model_results["model_id"].nunique()) == 101,
            "gene_symbols_unique": int(gene_stability["gene_symbol"].nunique()) == len(overlap),
            "output_gene_rows": int(len(gene_stability)),
            "external_auc_range": [float(model_results["external_mean_roc_auc"].min()), float(model_results["external_mean_roc_auc"].max())],
            "stable_ml_count": int(len(stable)),
            "stable_ml_tier1_overlap_count": int(len(overlap_stable_tier1)),
        },
    }
    (OUTPUTS / "DINP_CRC_101ML_log.md").write_text(
        "# DINP–CRC 101-model ML audit log\n\n```json\n"
        + json.dumps(log, indent=2, ensure_ascii=False)
        + "\n```\n",
        encoding="utf-8",
    )

    print(f"Completed 101 models; measured features={len(measured_genes)}; stable ML genes={len(stable)}")
    print("Top 10 ML genes:")
    print(gene_stability[["ml_rank", "gene_symbol", "ml_priority_score", "model_selection_frequency_selective_90", "mean_external_roc_auc_when_selected"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
