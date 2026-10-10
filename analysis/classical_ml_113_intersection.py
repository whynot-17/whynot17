from __future__ import annotations

import importlib.util
import json
import sys
from datetime import date, datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for dependency_path in reversed([ROOT / "work" / "ml-deps", ROOT / "work" / "enrichment-deps"]):
    sys.path.insert(0, str(dependency_path))

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC


ANALYSIS_DIR = ROOT / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR))
import validate_113_across_geo_replacement as geo  # noqa: E402


OUT = ROOT / "results" / "2026-10-11" / "classical_ml_113_4cohorts"
INPUT_113 = ROOT / "results" / "2026-10-10" / "DINP_286_x_TCGA_GDC_FDR_7748.tsv"
COHORTS = ("GSE32323", "GSE44076", "GSE41258", "GSE106582")
ALGORITHMS = ("LASSO_logistic", "Random_Forest", "SVM_RFE")
RANDOM_SEED = 20261011
N_INNER_FOLDS = 3
N_BOOTSTRAP = 1000
RF_TREES = 400
FEATURE_COUNTS = (5, 10, 20, 40, 80)
LASSO_C_GRID = (0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0)


def load_cohort(accession: str, targets: list[str]) -> dict:
    metadata, table, matrix_sha = geo.load_series_matrix(accession)
    records = geo.assign_groups(accession, geo.sample_records(metadata))
    sample_ids = list(table.columns[1:])
    by_gsm = {record["gsm"]: record for record in records}
    if set(sample_ids) - set(by_gsm):
        raise ValueError(f"{accession}: sample metadata do not match expression columns")
    records = [by_gsm[gsm] for gsm in sample_ids]
    platform = geo.PLATFORM_BY_GSE[accession]
    probe_map, annotation_sha = geo.platform_mapping(platform)
    genes, matrix, probe_counts, transformed = geo.map_probe_matrix(table, probe_map, accession)
    selected = [i for i, rec in enumerate(records) if rec["group"] in {"Tumor", "Normal"}]
    pair_ids = geo.complete_pair_ids(records, selected)
    selected = [i for i in selected if records[i]["pair_id"] in set(pair_ids)]
    if not pair_ids:
        raise ValueError(f"{accession}: no complete patient pairs")
    gene_index = {symbol: index for index, symbol in enumerate(genes)}
    represented = {
        gene for gene in targets
        if gene in gene_index and np.isfinite(matrix[gene_index[gene], selected]).any()
    }
    return {
        "accession": accession,
        "platform": platform,
        "records": records,
        "matrix": matrix,
        "genes": genes,
        "gene_index": gene_index,
        "probe_counts": probe_counts,
        "selected_indices": selected,
        "pair_ids": pair_ids,
        "represented": represented,
        "matrix_sha256": matrix_sha,
        "annotation_sha256": annotation_sha,
        "transformed": transformed,
    }


def assemble_samples(cohort: dict, features: list[str]) -> tuple[np.ndarray, pd.DataFrame]:
    accession = cohort["accession"]
    records = cohort["records"]
    matrix = cohort["matrix"]
    feature_indices = [cohort["gene_index"][gene] for gene in features]
    rows, vectors = [], []
    selected = cohort["selected_indices"]
    for patient in cohort["pair_ids"]:
        for group in ("Tumor", "Normal"):
            sample_indices = [
                index for index in selected
                if records[index]["pair_id"] == patient and records[index]["group"] == group
            ]
            if not sample_indices:
                continue
            values = np.nanmean(matrix[feature_indices, :][:, sample_indices], axis=1)
            vectors.append(values)
            rows.append({
                "accession": accession,
                "platform": cohort["platform"],
                "pair_id": patient,
                "group": group,
                "label": int(group == "Tumor"),
                "group_id": f"{accession}::{patient}",
                "source_sample_count_aggregated": len(sample_indices),
            })
    x = np.asarray(vectors, dtype=float)
    if x.ndim != 2 or x.shape[1] != len(features):
        raise ValueError(f"{accession}: invalid assembled expression shape {x.shape}")
    # Impute isolated missing values within cohort, then scale without using labels.
    medians = np.nanmedian(x, axis=0)
    bad = ~np.isfinite(medians)
    if bad.any():
        raise ValueError(f"{accession}: all-NA common genes: {[features[i] for i in np.where(bad)[0]]}")
    missing_row, missing_col = np.where(~np.isfinite(x))
    if len(missing_row):
        x[missing_row, missing_col] = medians[missing_col]
    means = x.mean(axis=0)
    std = x.std(axis=0, ddof=1)
    flat = std <= 1e-12
    std[flat] = 1.0
    x = (x - means) / std
    x[:, flat] = 0.0
    return x, pd.DataFrame(rows)


def inner_splits(x: np.ndarray, y: np.ndarray, groups: np.ndarray, seed: int):
    splitter = StratifiedGroupKFold(
        n_splits=N_INNER_FOLDS, shuffle=True, random_state=seed
    )
    return list(splitter.split(x, y, groups))


def make_logistic(c_value: float, seed: int) -> LogisticRegression:
    return LogisticRegression(
        C=c_value,
        l1_ratio=1.0,
        solver="liblinear",
        class_weight="balanced",
        max_iter=10000,
        random_state=seed,
    )


def fit_lasso(x_train, y_train, groups_train, x_test, seed):
    splits = inner_splits(x_train, y_train, groups_train, seed)
    cv_scores = []
    for c_value in LASSO_C_GRID:
        fold_scores = []
        for train_ix, valid_ix in splits:
            scaler = StandardScaler().fit(x_train[train_ix])
            x_inner_train = scaler.transform(x_train[train_ix])
            x_inner_valid = scaler.transform(x_train[valid_ix])
            model = make_logistic(c_value, seed)
            model.fit(x_inner_train, y_train[train_ix])
            fold_scores.append(roc_auc_score(y_train[valid_ix], model.predict_proba(x_inner_valid)[:, 1]))
        cv_scores.append(float(np.mean(fold_scores)))
    best_score = max(cv_scores)
    best_index = min(
        (i for i, value in enumerate(cv_scores) if np.isclose(value, best_score)),
        key=lambda i: LASSO_C_GRID[i],
    )
    best_c = LASSO_C_GRID[best_index]
    scaler = StandardScaler().fit(x_train)
    model = make_logistic(best_c, seed)
    model.fit(scaler.transform(x_train), y_train)
    scores = model.predict_proba(scaler.transform(x_test))[:, 1]
    coefficients = model.coef_.ravel()
    selected = np.flatnonzero(np.abs(coefficients) > 1e-8)
    order = selected[np.argsort(-np.abs(coefficients[selected]), kind="stable")]
    return scores, order, np.abs(coefficients), {"best_C": best_c, "inner_mean_AUC": best_score}


def make_rf(seed: int) -> RandomForestClassifier:
    return RandomForestClassifier(
        n_estimators=RF_TREES,
        max_features="sqrt",
        min_samples_leaf=2,
        class_weight="balanced_subsample",
        n_jobs=-1,
        random_state=seed,
    )


def fit_random_forest(x_train, y_train, groups_train, x_test, seed):
    splits = inner_splits(x_train, y_train, groups_train, seed)
    candidate_ks = [k for k in FEATURE_COUNTS if k <= x_train.shape[1]]
    score_by_k = {k: [] for k in candidate_ks}
    for fold, (train_ix, valid_ix) in enumerate(splits):
        ranker = make_rf(seed + fold + 1)
        ranker.fit(x_train[train_ix], y_train[train_ix])
        order = np.argsort(-ranker.feature_importances_, kind="stable")
        for k in candidate_ks:
            selected = order[:k]
            model = make_rf(seed + 100 + fold)
            model.fit(x_train[train_ix][:, selected], y_train[train_ix])
            fold_score = roc_auc_score(
                y_train[valid_ix], model.predict_proba(x_train[valid_ix][:, selected])[:, 1]
            )
            score_by_k[k].append(float(fold_score))
    means = {k: float(np.mean(values)) for k, values in score_by_k.items()}
    best_k = max(candidate_ks, key=lambda k: (means[k], -k))
    ranker = make_rf(seed + 200)
    ranker.fit(x_train, y_train)
    order = np.argsort(-ranker.feature_importances_, kind="stable")
    selected = order[:best_k]
    model = make_rf(seed + 201)
    model.fit(x_train[:, selected], y_train)
    scores = model.predict_proba(x_test[:, selected])[:, 1]
    return scores, selected, ranker.feature_importances_, {"best_k": best_k, "inner_mean_AUC": means[best_k]}


def make_linear_svc(seed: int) -> LinearSVC:
    return LinearSVC(
        C=1.0,
        class_weight="balanced",
        dual="auto",
        max_iter=10000,
        random_state=seed,
    )


def rfe_order(x_train_scaled: np.ndarray, y_train: np.ndarray, seed: int) -> np.ndarray:
    min_k = min(5, x_train_scaled.shape[1])
    selector = RFE(
        estimator=make_linear_svc(seed),
        n_features_to_select=min_k,
        step=0.1,
    )
    selector.fit(x_train_scaled, y_train)
    return np.argsort(selector.ranking_, kind="stable")


def fit_svm_rfe(x_train, y_train, groups_train, x_test, seed):
    splits = inner_splits(x_train, y_train, groups_train, seed)
    candidate_ks = [k for k in FEATURE_COUNTS if k <= x_train.shape[1]]
    score_by_k = {k: [] for k in candidate_ks}
    for fold, (train_ix, valid_ix) in enumerate(splits):
        scaler = StandardScaler().fit(x_train[train_ix])
        x_inner_train = scaler.transform(x_train[train_ix])
        x_inner_valid = scaler.transform(x_train[valid_ix])
        order = rfe_order(x_inner_train, y_train[train_ix], seed + fold)
        for k in candidate_ks:
            selected = order[:k]
            model = make_linear_svc(seed + 100 + fold)
            model.fit(x_inner_train[:, selected], y_train[train_ix])
            fold_score = roc_auc_score(
                y_train[valid_ix], model.decision_function(x_inner_valid[:, selected])
            )
            score_by_k[k].append(float(fold_score))
    means = {k: float(np.mean(values)) for k, values in score_by_k.items()}
    best_k = max(candidate_ks, key=lambda k: (means[k], -k))
    scaler = StandardScaler().fit(x_train)
    x_train_scaled = scaler.transform(x_train)
    x_test_scaled = scaler.transform(x_test)
    order = rfe_order(x_train_scaled, y_train, seed + 200)
    selected = order[:best_k]
    model = make_linear_svc(seed + 201)
    model.fit(x_train_scaled[:, selected], y_train)
    scores = model.decision_function(x_test_scaled[:, selected])
    coefficient_scores = np.zeros(x_train.shape[1], dtype=float)
    coefficient_scores[selected] = np.abs(model.coef_.ravel())
    selected = selected[np.argsort(-coefficient_scores[selected], kind="stable")]
    return scores, selected, coefficient_scores, {"best_k": best_k, "inner_mean_AUC": means[best_k]}


def metric_summary(y_true, scores, threshold, groups, seed):
    predictions = (scores >= threshold).astype(int)
    tn = int(np.sum((y_true == 0) & (predictions == 0)))
    fp = int(np.sum((y_true == 0) & (predictions == 1)))
    fn = int(np.sum((y_true == 1) & (predictions == 0)))
    tp = int(np.sum((y_true == 1) & (predictions == 1)))
    observed_auc = float(roc_auc_score(y_true, scores))
    unique_groups = np.unique(groups)
    index_by_group = {group: np.flatnonzero(groups == group) for group in unique_groups}
    rng = np.random.default_rng(seed)
    boot_aucs = []
    for _ in range(N_BOOTSTRAP):
        sampled_groups = rng.choice(unique_groups, size=len(unique_groups), replace=True)
        indices = np.concatenate([index_by_group[group] for group in sampled_groups])
        boot_y = y_true[indices]
        if np.unique(boot_y).size == 2:
            boot_aucs.append(roc_auc_score(boot_y, scores[indices]))
    ci_low, ci_high = np.percentile(boot_aucs, [2.5, 97.5]) if boot_aucs else (np.nan, np.nan)
    return {
        "n_samples": int(len(y_true)),
        "n_patients": int(len(unique_groups)),
        "AUC": observed_auc,
        "AUC_pair_bootstrap_95CI_low": float(ci_low),
        "AUC_pair_bootstrap_95CI_high": float(ci_high),
        "accuracy": float(accuracy_score(y_true, predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, predictions)),
        "sensitivity_tumor": float(recall_score(y_true, predictions, pos_label=1, zero_division=0)),
        "specificity_normal": float(recall_score(y_true, predictions, pos_label=0, zero_division=0)),
        "precision_tumor": float(precision_score(y_true, predictions, pos_label=1, zero_division=0)),
        "F1_tumor": float(f1_score(y_true, predictions, pos_label=1, zero_division=0)),
        "TN": tn,
        "FP": fp,
        "FN": fn,
        "TP": tp,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    query_frame = pd.read_csv(INPUT_113, sep="\t")
    targets = sorted(set(query_frame["gene_symbol"].astype(str).str.strip().str.upper()) - {""})
    if len(targets) != 113:
        raise ValueError(f"Expected 113 unique query genes, got {len(targets)}")

    loaded = {accession: load_cohort(accession, targets) for accession in COHORTS}
    common_features = sorted(set.intersection(*(loaded[a]["represented"] for a in COHORTS)))
    if len(common_features) < 20:
        raise ValueError(f"Too few genes represented in every cohort: {len(common_features)}")
    missing_union = sorted(set(targets) - set(common_features))

    x_parts, metadata_parts, cohort_summary, checksums = [], [], [], {
        "113_intersection_table_sha256": geo.sha256(INPUT_113),
    }
    for accession in COHORTS:
        cohort = loaded[accession]
        x_cohort, samples = assemble_samples(cohort, common_features)
        x_parts.append(x_cohort)
        metadata_parts.append(samples)
        cohort_summary.append({
            "accession": accession,
            "platform": cohort["platform"],
            "complete_pairs_used": len(cohort["pair_ids"]),
            "samples_used": len(samples),
            "tumor_samples": int(samples["label"].sum()),
            "normal_samples": int((samples["label"] == 0).sum()),
            "candidate_genes": len(targets),
            "genes_measured_on_platform": len(cohort["represented"]),
            "common_features_used": len(common_features),
            "expression_transform": "log2(PLIER signal)" if cohort["transformed"] else "GEO series-matrix values retained",
            "normalization": "within-cohort gene-wise z-score over complete-pair samples; labels not used",
        })
        checksums[f"{accession}_series_matrix_sha256"] = cohort["matrix_sha256"]
        checksums[f"{cohort['platform']}_annotation_sha256"] = cohort["annotation_sha256"]

    x = np.vstack(x_parts)
    sample_metadata = pd.concat(metadata_parts, ignore_index=True)
    y = sample_metadata["label"].to_numpy(dtype=int)
    groups = sample_metadata["group_id"].astype(str).to_numpy()
    accessions = sample_metadata["accession"].astype(str).to_numpy()
    feature_ix = {gene: i for i, gene in enumerate(common_features)}

    performance_rows, selection_rows, prediction_rows = [], [], []
    for outer_index, heldout in enumerate(COHORTS):
        test_mask = accessions == heldout
        train_mask = ~test_mask
        x_train, y_train, groups_train = x[train_mask], y[train_mask], groups[train_mask]
        x_test, y_test, groups_test = x[test_mask], y[test_mask], groups[test_mask]
        test_samples = sample_metadata.loc[test_mask].reset_index(drop=True)
        train_pairs = len(np.unique(groups_train))
        for algorithm_index, algorithm in enumerate(ALGORITHMS):
            seed = RANDOM_SEED + outer_index * 1000 + algorithm_index * 100
            if algorithm == "LASSO_logistic":
                scores, selected_ix, importance, tune = fit_lasso(
                    x_train, y_train, groups_train, x_test, seed
                )
                threshold = 0.5
            elif algorithm == "Random_Forest":
                scores, selected_ix, importance, tune = fit_random_forest(
                    x_train, y_train, groups_train, x_test, seed
                )
                threshold = 0.5
            else:
                scores, selected_ix, importance, tune = fit_svm_rfe(
                    x_train, y_train, groups_train, x_test, seed
                )
                threshold = 0.0
            metrics = metric_summary(y_test, np.asarray(scores), threshold, groups_test, seed + 700)
            performance_rows.append({
                "algorithm": algorithm,
                "heldout_cohort": heldout,
                "n_training_patient_pairs": train_pairs,
                "n_selected_features": len(selected_ix),
                "selected_feature_count_tuning": tune,
                **metrics,
            })
            for rank, gene_index in enumerate(selected_ix, start=1):
                selection_rows.append({
                    "algorithm": algorithm,
                    "heldout_cohort": heldout,
                    "gene_symbol": common_features[int(gene_index)],
                    "selection_rank": rank,
                    "feature_score": float(importance[int(gene_index)]),
                    "selected": True,
                })
            for j, row in test_samples.iterrows():
                prediction_rows.append({
                    **row.to_dict(),
                    "algorithm": algorithm,
                    "heldout_cohort": heldout,
                    "prediction_score": float(scores[j]),
                    "predicted_group": "Tumor" if scores[j] >= threshold else "Normal",
                    "threshold": threshold,
                    "n_selected_features": len(selected_ix),
                })

    performance = pd.DataFrame(performance_rows)
    selection = pd.DataFrame(selection_rows)
    predictions = pd.DataFrame(prediction_rows)
    frequency = (
        selection.groupby(["algorithm", "gene_symbol"], as_index=False)
        .agg(selected_in_outer_training_folds=("heldout_cohort", "nunique"),
             mean_selection_rank=("selection_rank", "mean"),
             mean_feature_score=("feature_score", "mean"))
    )
    all_features = pd.MultiIndex.from_product(
        [ALGORITHMS, common_features], names=["algorithm", "gene_symbol"]
    ).to_frame(index=False)
    frequency = all_features.merge(frequency, on=["algorithm", "gene_symbol"], how="left")
    frequency["selected_in_outer_training_folds"] = frequency[
        "selected_in_outer_training_folds"
    ].fillna(0).astype(int)
    frequency["selected_in_at_least_3_of_4_folds"] = frequency[
        "selected_in_outer_training_folds"
    ] >= 3
    consensus = (
        frequency.pivot(index="gene_symbol", columns="algorithm", values="selected_in_outer_training_folds")
        .reset_index()
    )
    consensus["algorithms_selected_in_at_least_3_folds"] = sum(
        (consensus[algorithm] >= 3).astype(int) for algorithm in ALGORITHMS
    )
    consensus["consensus_selected_by_2plus_algorithms"] = (
        consensus["algorithms_selected_in_at_least_3_folds"] >= 2
    )

    pooled_rows = []
    for algorithm in ALGORITHMS:
        part = predictions[predictions["algorithm"] == algorithm].copy()
        part_groups = part["group_id"].astype(str).to_numpy()
        scores = part["prediction_score"].to_numpy(dtype=float)
        labels = part["label"].to_numpy(dtype=int)
        threshold = float(part["threshold"].iloc[0])
        metrics = metric_summary(labels, scores, threshold, part_groups, RANDOM_SEED + 9000 + ALGORITHMS.index(algorithm))
        cohort_auc = performance.loc[performance["algorithm"] == algorithm, "AUC"]
        pooled_rows.append({
            "algorithm": algorithm,
            "heldout_scope": "pooled_leave_one_cohort_out_predictions",
            "mean_cohort_AUC": float(cohort_auc.mean()),
            **metrics,
        })
    pooled = pd.DataFrame(pooled_rows)

    cohort_summary_frame = pd.DataFrame(cohort_summary)
    coverage_rows = []
    for accession in COHORTS:
        for gene in targets:
            coverage_rows.append({
                "accession": accession,
                "gene_symbol": gene,
                "has_usable_probe": gene in loaded[accession]["represented"],
                "included_in_cross_cohort_ML": gene in feature_ix,
                "probe_count": loaded[accession]["probe_counts"].get(gene, 0),
            })
    coverage = pd.DataFrame(coverage_rows)

    OUT.mkdir(parents=True, exist_ok=True)
    cohort_summary_frame.to_csv(OUT / "cohort_sample_summary.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "gene_platform_coverage.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_symbol": common_features}).to_csv(OUT / "features_used_common_coverage.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_symbol": missing_union}).to_csv(OUT / "candidate_genes_not_measured_in_all_cohorts.tsv", sep="\t", index=False)
    sample_metadata.to_csv(OUT / "paired_sample_metadata.tsv", sep="\t", index=False)
    performance.to_csv(OUT / "leave_one_cohort_out_metrics.tsv", sep="\t", index=False)
    pooled.to_csv(OUT / "pooled_out_of_cohort_metrics.tsv", sep="\t", index=False)
    selection.to_csv(OUT / "selected_genes_by_algorithm_and_fold.tsv", sep="\t", index=False)
    frequency.to_csv(OUT / "feature_selection_frequency.tsv", sep="\t", index=False)
    consensus.to_csv(OUT / "algorithm_consensus_genes.tsv", sep="\t", index=False)
    predictions.to_csv(OUT / "leave_one_cohort_out_predictions.tsv", sep="\t", index=False)
    (OUT / "source_checksums.json").write_text(json.dumps(checksums, indent=2), encoding="utf-8")

    summary = {
        "analysis_date": date.today().isoformat(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_gene_count": len(targets),
        "common_measured_gene_count": len(common_features),
        "candidate_genes_not_measured_in_all_cohorts": missing_union,
        "cohorts": COHORTS,
        "paired_samples_used": int(len(sample_metadata)),
        "complete_patient_pairs_total": int(sample_metadata["group_id"].nunique()),
        "algorithms": {
            "LASSO_logistic": "L1 logistic regression; C selected by 3-fold StratifiedGroupKFold within training cohorts",
            "Random_Forest": f"RandomForestClassifier; top-k feature count selected by grouped inner CV; {RF_TREES} trees",
            "SVM_RFE": "linear SVM recursive feature elimination; feature count selected by grouped inner CV",
        },
        "outer_validation": "leave-one-GEO-cohort-out; all samples from the held-out cohort are excluded from training",
        "inner_validation_group": "patient pair; paired tumor and normal samples are kept in the same inner fold",
        "preprocessing": "probe-to-symbol mapping and per-gene median probe aggregation as prior current-project GEO pipeline; use complete patient pairs; within-cohort gene-wise z-score without labels",
        "metrics": "AUC, pair-cluster bootstrap 95% CI, accuracy, balanced accuracy, sensitivity, specificity, precision, F1",
        "random_seed": RANDOM_SEED,
        "checksums": checksums,
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = [
        "# 113个DINP×CRC交集基因的经典三算法机器学习",
        "",
        f"- 分析日期：{summary['analysis_date']}",
        "- 分类任务：CRC肿瘤 vs 正常组织。",
        f"- 候选集：113个DINP×TCGA/GDC FDR-DEG交集基因；跨队列共同可测的特征用于模型：**{len(common_features)}个**。",
        f"- 平台覆盖差异：有{len(missing_union)}个候选基因至少在一个GEO平台没有可用探针，未纳入跨队列模型；详见 `candidate_genes_not_measured_in_all_cohorts.tsv`。",
        "- 队列：GSE32323、GSE44076、GSE41258、GSE106582。只保留完整患者配对；同一患者的肿瘤/正常样本在内层交叉验证中始终作为一个分组。",
        "- 预处理：沿用当前项目的探针注释和多探针中位数汇总；GSE41258按既有流程做log2变换；每个队列按基因做无标签z-score，再进入模型。",
        "- 外层评估：留一队列外部验证（LOCO），每轮整个队列只用于测试，其余三个队列训练。模型特征选择和调参仅在训练数据中完成。",
        "- 三算法：LASSO-logistic（L1）、Random Forest、linear SVM-RFE；LASSO的C及RF/SVM-RFE的特征数使用按患者分组的内层3折CV选择。",
        "- AUC置信区间使用患者配对为单位的bootstrap；性能仅反映这四个GEO队列的肿瘤/正常区分，不等于临床诊断效能。",
        "",
        "## 队列样本",
        "",
        "| 队列 | 完整患者对 | 纳入样本 | 肿瘤 | 正常 | 共同特征数 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in cohort_summary_frame.iterrows():
        lines.append(
            f"| {row['accession']} | {row['complete_pairs_used']} | {row['samples_used']} | "
            f"{row['tumor_samples']} | {row['normal_samples']} | {row['common_features_used']} |"
        )
    lines.extend([
        "",
        "## 留一队列验证性能",
        "",
        "| 算法 | 测试队列 | AUC (患者配对bootstrap 95% CI) | 准确率 | 平衡准确率 | 灵敏度 | 特异度 | 特征数 |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ])
    for _, row in performance.iterrows():
        lines.append(
            f"| {row['algorithm']} | {row['heldout_cohort']} | "
            f"{row['AUC']:.3f} ({row['AUC_pair_bootstrap_95CI_low']:.3f}–{row['AUC_pair_bootstrap_95CI_high']:.3f}) | "
            f"{row['accuracy']:.3f} | {row['balanced_accuracy']:.3f} | "
            f"{row['sensitivity_tumor']:.3f} | {row['specificity_normal']:.3f} | {row['n_selected_features']} |"
        )
    lines.extend(["", "## 三算法共同/稳定特征", ""])
    stable_consensus = consensus[consensus["consensus_selected_by_2plus_algorithms"]].copy()
    if stable_consensus.empty:
        lines.append("没有基因同时满足：至少两个算法在四轮训练中各自有≥3轮选中。完整入选频率见 `feature_selection_frequency.tsv`。")
    else:
        lines.append("以下基因至少被两个算法在四轮外层训练中各自选中≥3次：")
        lines.append(", ".join(stable_consensus["gene_symbol"].astype(str).tolist()))
    lines.extend([
        "",
        "## 结果文件",
        "",
        "- `leave_one_cohort_out_metrics.tsv`：每个算法在每个独立留出队列的指标。",
        "- `pooled_out_of_cohort_metrics.tsv`：四个外部留出队列预测合并后的汇总指标。",
        "- `leave_one_cohort_out_predictions.tsv`：逐样本外部预测分数和分类结果。",
        "- `selected_genes_by_algorithm_and_fold.tsv`：每轮训练中各算法选出的基因及排序。",
        "- `feature_selection_frequency.tsv`、`algorithm_consensus_genes.tsv`：跨轮特征稳定性和算法共识。",
        "- `gene_platform_coverage.tsv`、`cohort_sample_summary.tsv`：平台覆盖及配对样本量。",
        "- `source_checksums.json`、`analysis_summary.json`：源文件校验和及模型配置。",
        "",
    ])
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps({
        "candidate_genes": len(targets),
        "common_features": len(common_features),
        "not_measured_in_all_cohorts": missing_union,
        "samples": len(sample_metadata),
        "pairs": int(sample_metadata["group_id"].nunique()),
        "performance": performance[["algorithm", "heldout_cohort", "AUC", "accuracy", "balanced_accuracy", "n_selected_features"]].to_dict(orient="records"),
        "consensus": stable_consensus["gene_symbol"].astype(str).tolist(),
        "output": str(OUT),
    }, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
