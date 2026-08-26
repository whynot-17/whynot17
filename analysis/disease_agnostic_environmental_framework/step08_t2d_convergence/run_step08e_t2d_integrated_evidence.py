#!/usr/bin/env python3
"""Step 8E: integrated evidence profiles and flagship classification.

This is a read-only synthesis of frozen Step 5--8D outputs.  It does not
re-run an upstream analysis, add a data source, change a Tier A definition,
or combine domains into an opaque numerical score.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STEP5 = ROOT / "step05_t2d_screen"
STEP6 = ROOT / "step06_t2d_robustness"
STEP7 = ROOT / "step07_genecard_convergence"
STEP8D = ROOT / "step08_transcriptomic_convergence"

AXES = {
    "cluster_5": {"biomarkers": ["URXP02"], "axis_label": "URXP02"},
    "cluster_6": {"biomarkers": ["URXUBA", "URXUSR"], "axis_label": "URXUBA + URXUSR"},
    "cluster_8": {"biomarkers": ["URXUPB"], "axis_label": "URXUPB"},
    "cluster_11": {"biomarkers": ["URXUUR"], "axis_label": "URXUUR"},
}

# These are qualitative classifications required by the locked Step 8E plan.
# They are documented in the output and are not a post-hoc numerical score.
CLASSIFICATION = {
    "cluster_5": {
        "final_classification": "Flagship",
        "pathway_concentration": "High",
        "transcriptomic_support": "Moderate, tissue-specific support",
        "classification_rationale": (
            "Most concentrated cross-layer profile: FDR-supported and robust epidemiology, "
            "strong GeneCards convergence, a focused xenobiotic/CYP pathway theme, strong "
            "network enrichment, and interpretable human tissue support. Cycle direction is "
            "less complete than the other axes and the gene set is comparatively small."
        ),
    },
    "cluster_6": {
        "final_classification": "Supported",
        "pathway_concentration": "Moderate",
        "transcriptomic_support": "Moderate, tissue-specific support",
        "classification_rationale": (
            "Strong epidemiology, robustness, GeneCards enrichment, and network structure; "
            "pathway themes center on miRNA/interleukin/IL-17/MAPK biology but remain broad, "
            "and transcriptomic direction differs by tissue. Retained as a two-biomarker axis."
        ),
    },
    "cluster_8": {
        "final_classification": "Exploratory",
        "pathway_concentration": "Low",
        "transcriptomic_support": "Limited and heterogeneous",
        "classification_rationale": (
            "Epidemiology, robustness, GeneCards enrichment, and network enrichment are "
            "present, but the very large input gene set, broad RNA/protein/cell-cycle themes, "
            "and mixed tissue directionality limit a focused biological interpretation."
        ),
    },
    "cluster_11": {
        "final_classification": "Supported",
        "pathway_concentration": "Moderate",
        "transcriptomic_support": "Moderate, tissue-specific support",
        "classification_rationale": (
            "Strong epidemiology, robustness, GeneCards enrichment, and network structure; "
            "pathway results converge on transport/trafficking and protein-related themes but "
            "are not as focused as cluster_5. Transcriptomic support is clearly tissue-specific."
        ),
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fmt(value: object, digits: int = 3) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not math.isfinite(number):
        return "NA"
    return f"{number:.{digits}g}"


def fmt_p(value: object) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not math.isfinite(number):
        return "NA"
    return f"{number:.3g}"


def str_or_blank(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value)


def bool_text(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def read_inputs() -> dict[str, pd.DataFrame]:
    paths = {
        "step5": STEP5 / "t2d_primary_29_tests.csv",
        "step6": STEP6 / "t2d_robustness_results.csv",
        "step6_cycle": STEP6 / "t2d_cycle_heterogeneity.csv",
        "step7_joint": STEP7 / "t2d_step7_joint_prioritization.csv",
        "step7_enrichment": STEP7 / "t2d_cluster_enrichment_primary.csv",
        "step8a_modules": HERE / "t2d_step8_module_summary.csv",
        "step8b_representatives": HERE / "t2d_step8_module_representatives.csv",
        "step8c_summary": HERE / "t2d_step8c_network_summary.csv",
        "step8c_nodes": HERE / "t2d_step8c_network_nodes.csv",
        "step8c_modules": HERE / "t2d_step8c_network_modules.csv",
        "step8d_axis": STEP8D / "t2d_step8d_axis_summary.csv",
        "step8d_dataset_axis": STEP8D / "t2d_step8d_dataset_axis_summary.csv",
        "step8d_synthesis": STEP8D / "t2d_step8d_module_cross_dataset_synthesis.csv",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing frozen input(s): " + "; ".join(missing))
    return {name: pd.read_csv(path) for name, path in paths.items()}


def build_transcriptomic_summary(axis_id: str, data: dict[str, pd.DataFrame], network_modules: pd.DataFrame) -> dict[str, object]:
    axis = data["step8d_axis"].loc[data["step8d_axis"]["cluster_id"] == axis_id]
    dataset_axis = data["step8d_dataset_axis"].loc[data["step8d_dataset_axis"]["cluster_id"] == axis_id]
    synthesis = data["step8d_synthesis"].merge(network_modules[["module_id", "cluster_id"]], on="module_id", how="left")
    synthesis = synthesis.loc[synthesis["cluster_id"] == axis_id]
    n_modules_total = int(network_modules.loc[network_modules["cluster_id"] == axis_id, "module_id"].nunique())
    n_modules_tested = int((pd.to_numeric(synthesis["n_datasets_tested"], errors="coerce") > 0).sum())
    n_estimable = int(pd.to_numeric(synthesis["n_datasets_tested"], errors="coerce").sum())
    n_q = int(pd.to_numeric(synthesis["n_dataset_q_lt_0_05"], errors="coerce").sum())
    positive_axis = int((dataset_axis["n_modules_positive"] > dataset_axis["n_modules_negative"]).sum())
    negative_axis = int((dataset_axis["n_modules_negative"] > dataset_axis["n_modules_positive"]).sum())
    tissue_bits = []
    for row in dataset_axis.sort_values("accession").itertuples(index=False):
        tissue_bits.append(
            f"{row.accession}/{row.tissue}: {int(row.n_modules_positive)}+/{int(row.n_modules_negative)}-; "
            f"median Δ={fmt(row.median_module_delta)}; q<0.05={int(row.n_modules_q_lt_0_05)}"
        )
    return {
        "transcriptomic_frozen_modules": n_modules_total,
        "transcriptomic_modules_tested": n_modules_tested,
        "transcriptomic_estimable_dataset_module_results": n_estimable,
        "transcriptomic_modules_with_dataset_q_lt_0_05": n_q,
        "transcriptomic_dataset_axis_majority_positive": positive_axis,
        "transcriptomic_dataset_axis_majority_negative": negative_axis,
        "transcriptomic_median_axis_delta": float(axis["median_module_delta"].iloc[0]) if len(axis) else math.nan,
        "transcriptomic_tissue_summary": " | ".join(tissue_bits),
    }


def main() -> None:
    data = read_inputs()
    for axis_id in AXES:
        if axis_id not in set(data["step7_joint"]["cluster_id"].astype(str)):
            raise ValueError(f"{axis_id} absent from frozen Step 7 joint prioritization")

    step5 = data["step5"]
    step6 = data["step6"]
    step6_cycle = data["step6_cycle"]
    step7_joint = data["step7_joint"]
    step7_enrichment = data["step7_enrichment"]
    pathway_modules = data["step8a_modules"]
    pathway_reps = data["step8b_representatives"]
    network_summary = data["step8c_summary"]
    network_nodes = data["step8c_nodes"]
    network_modules = data["step8c_modules"]

    profiles: list[dict[str, object]] = []
    final_rows: list[dict[str, object]] = []
    risk_rows: list[dict[str, object]] = []

    risk_catalog = {
        "cross_sectional_temporality": ("High", "Step 5 is survey-based and cross-sectional; association is not temporally ordered.", "Do not use causal or incidence language."),
        "reverse_causation_potential": ("Moderate", "Outcome and biomarker are observed within the public T2D screen rather than prospectively sequenced.", "Retain as an interpretation boundary; no mediation claim."),
        "broad_pathway_annotation": ("Axis-specific", "Annotation enrichment is direction-agnostic and can be redundant; concentration is reported separately.", "Use compact representatives and thematic language only."),
        "gene_set_size_inflation": ("Axis-specific", "Step 7 cluster gene-set sizes range from 65 to 3,016 genes.", "Report query size and GeneCards universe; avoid causal-gene wording."),
        "sparse_gene_instability": ("Axis-specific", "Small CTD/network inputs can make enrichment or modules sensitive to individual genes.", "Treat cluster_5 as focused but comparatively small; no single-gene causal claim."),
        "tissue_specific_transcriptomic_heterogeneity": ("High", "Step 8D directions vary across liver, muscle, adipose, and islet datasets.", "Keep tissue labels and do not call a universal T2D program."),
        "annotation_density_bias": ("Moderate", "Better-annotated genes and larger modules can accumulate more pathway/network evidence.", "Use common frozen rules and report mapping/module coverage."),
        "cluster_level_biomarker_dependence": ("Axis-specific", "cluster_6 contains two correlated biomarker tests and must remain one exposure axis.", "Do not count URXUBA and URXUSR as independent flagship discoveries."),
    }

    for axis_id, axis_config in AXES.items():
        biomarkers = axis_config["biomarkers"]
        e5 = step5.loc[step5["biomarker"].isin(biomarkers)].copy()
        e6 = step6.loc[step6["biomarker"].isin(biomarkers)].copy()
        e6c = step6_cycle.loc[step6_cycle["biomarker"].isin(biomarkers)].copy()
        e7j = step7_joint.loc[step7_joint["cluster_id"] == axis_id].iloc[0]
        e7 = step7_enrichment.loc[step7_enrichment["cluster_id"] == axis_id].iloc[0]
        p8a = pathway_modules.loc[pathway_modules["cluster_id"] == axis_id].copy()
        p8b = pathway_reps.loc[pathway_reps["cluster_id"] == axis_id].copy()
        p8b = p8b.loc[p8b["selected_for_compact_summary"].map(bool_text)].sort_values("selected_rank")
        n8c = network_summary.loc[network_summary["cluster_id"] == axis_id].iloc[0]
        n8cn = network_nodes.loc[network_nodes["cluster_id"] == axis_id].copy()
        n8cn["_top"] = n8cn["top_network_prioritized"].map(bool_text)
        top_genes = n8cn.loc[n8cn["_top"]].sort_values("network_priority_score", ascending=False)["preferred_name"].astype(str).head(6).tolist()
        n8cm = network_modules.loc[network_modules["cluster_id"] == axis_id].copy()
        n8cm["_annotations"] = pd.to_numeric(n8cm["n_significant_annotations"], errors="coerce").fillna(0)
        top_network_modules = [
            f"{row.module_id}: {str_or_blank(row.top_annotation)}"
            for row in n8cm.sort_values("_annotations", ascending=False).itertuples(index=False)
        ][:4]
        transcript = build_transcriptomic_summary(axis_id, data, network_modules)

        epi_details = []
        for row in e5.sort_values("biomarker").itertuples(index=False):
            epi_details.append(
                f"{row.biomarker}: OR {fmt(row.OR)} (95% CI {fmt(row.CI_low)}–{fmt(row.CI_high)}), "
                f"P={fmt_p(row.P)}, q={fmt_p(row.BH_FDR)}; N={int(row.N)}, T2D cases={int(row.T2D_cases)}"
            )
        robustness_details = []
        for row in e6.sort_values("biomarker").itertuples(index=False):
            cycle_p = e6c.loc[e6c["biomarker"] == row.biomarker, "interaction_P_F"]
            cycle_p_text = fmt_p(cycle_p.iloc[0]) if len(cycle_p) else "NA"
            discordant = str_or_blank(row.cycle_discordant) or "none listed"
            robustness_details.append(
                f"{row.biomarker}: LOCO {int(row.loco_same_n)}/{int(row.loco_n)}, "
                f"cycle direction {int(row.cycle_same_n)}/{int(row.cycle_n)}, "
                f"Pinteraction={cycle_p_text}, discordant={discordant}"
            )
        pathway_representative_text = "; ".join(p8b["representative_term"].astype(str).head(8).tolist())

        classification = CLASSIFICATION[axis_id]
        profile = {
            "cluster_id": axis_id,
            "axis_label": axis_config["axis_label"],
            "biomarkers": ";".join(biomarkers),
            "n_biomarkers": len(biomarkers),
            "epidemiology_strength": "Strong",
            "epidemiology_n_fdr_supported_tests": int((pd.to_numeric(e5["BH_FDR"], errors="coerce") < 0.05).sum()),
            "epidemiology_primary_details": " | ".join(epi_details),
            "robustness_strength": "Strong",
            "robustness_step6_priority_tiers": ";".join(sorted(e6["priority_tier"].astype(str).unique())),
            "robustness_details": " | ".join(robustness_details),
            "genecard_strength": "Strong",
            "genecard_universe_n": int(e7["gene_cards_k"]),
            "genecard_ctd_genes": int(e7["n_cluster_ctd_genes"]),
            "genecard_overlap": int(e7["n_overlap"]),
            "genecard_enrichment_or": float(e7["odds_ratio"]) if math.isfinite(float(e7["odds_ratio"])) else math.nan,
            "genecard_bh_fdr": float(e7["bh_fdr"]),
            "pathway_concentration": classification["pathway_concentration"],
            "pathway_significant_terms": int(pd.to_numeric(p8a["n_terms"], errors="coerce").sum()),
            "pathway_reduced_modules": len(p8a),
            "pathway_compact_representatives": len(p8b),
            "pathway_representatives": pathway_representative_text,
            "network_strength": "Strong",
            "network_input_genes": int(n8c["n_step7_overlap_genes"]),
            "network_nodes": int(n8c["n_network_nodes"]),
            "network_edges_score_ge_0_7": int(n8c["n_network_edges_score_ge_700"]),
            "network_observed_expected_edge_ratio": float(n8c["observed_to_expected_edge_ratio"]),
            "network_empirical_p": float(n8c["empirical_p_ge_observed"]),
            "network_louvain_modules": int(n8c["n_louvain_modules"]),
            "network_dominant_modules": "; ".join(top_network_modules),
            "network_prioritized_genes": ";".join(top_genes),
            "transcriptomic_support": classification["transcriptomic_support"],
            **transcript,
            "final_classification": classification["final_classification"],
            "classification_rationale": classification["classification_rationale"],
        }
        profiles.append(profile)
        final_rows.append({
            "cluster_id": axis_id,
            "axis_label": axis_config["axis_label"],
            "biomarkers": ";".join(biomarkers),
            "final_classification": classification["final_classification"],
            "pathway_concentration": classification["pathway_concentration"],
            "transcriptomic_support": classification["transcriptomic_support"],
            "classification_rationale": classification["classification_rationale"],
        })
        for risk_flag, (severity, evidence, handling) in risk_catalog.items():
            if risk_flag == "broad_pathway_annotation":
                evidence = f"{classification['pathway_concentration']} concentration; representatives: {pathway_representative_text}"
            elif risk_flag == "gene_set_size_inflation":
                evidence = f"Step 7 CTD gene set={int(e7['n_cluster_ctd_genes'])}; network input genes={int(n8c['n_step7_overlap_genes'])}."
            elif risk_flag == "sparse_gene_instability":
                evidence = f"Step 7 CTD gene set={int(e7['n_cluster_ctd_genes'])}; network nodes={int(n8c['n_network_nodes'])}."
            elif risk_flag == "cluster_level_biomarker_dependence":
                evidence = "Two frozen biomarker tests in this axis." if len(biomarkers) > 1 else "Single biomarker test in this axis."
            risk_rows.append({
                "cluster_id": axis_id,
                "axis_label": axis_config["axis_label"],
                "risk_flag": risk_flag,
                "severity": severity,
                "evidence": evidence,
                "handling": handling,
            })

    profiles_df = pd.DataFrame(profiles).sort_values("cluster_id")
    final_df = pd.DataFrame(final_rows).sort_values("cluster_id")
    risks_df = pd.DataFrame(risk_rows).sort_values(["cluster_id", "risk_flag"])
    profiles_df.to_csv(HERE / "t2d_step8e_axis_evidence_profiles.csv", index=False)
    final_df.to_csv(HERE / "t2d_step8e_final_classification.csv", index=False)
    risks_df.to_csv(HERE / "t2d_step8e_interpretation_risks.csv", index=False)

    input_paths = {
        "step5": STEP5 / "t2d_primary_29_tests.csv",
        "step6": STEP6 / "t2d_robustness_results.csv",
        "step6_cycle": STEP6 / "t2d_cycle_heterogeneity.csv",
        "step7_joint": STEP7 / "t2d_step7_joint_prioritization.csv",
        "step7_enrichment": STEP7 / "t2d_cluster_enrichment_primary.csv",
        "step8a_modules": HERE / "t2d_step8_module_summary.csv",
        "step8b_representatives": HERE / "t2d_step8_module_representatives.csv",
        "step8c_summary": HERE / "t2d_step8c_network_summary.csv",
        "step8c_nodes": HERE / "t2d_step8c_network_nodes.csv",
        "step8c_modules": HERE / "t2d_step8c_network_modules.csv",
        "step8d_axis": STEP8D / "t2d_step8d_axis_summary.csv",
        "step8d_dataset_axis": STEP8D / "t2d_step8d_dataset_axis_summary.csv",
        "step8d_synthesis": STEP8D / "t2d_step8d_module_cross_dataset_synthesis.csv",
    }
    manifest = {
        "status": "complete_integrated_evidence_profiles",
        "analysis": "Step 8E T2D integrated evidence synthesis and flagship selection",
        "frozen_axes": list(AXES),
        "upstream_chain": "29 frozen tests -> 14 FDR-positive -> 13 robust -> 11 clusters -> 5 GeneCards-enriched -> 4 Tier A -> pathway -> network -> transcriptomics -> integrated evidence",
        "no_new_analysis": True,
        "no_opaque_total_score": True,
        "classification_rule": "evidence-profile classification using the locked Flagship/Supported/Exploratory criteria; no domain scores summed",
        "classification": final_df.to_dict(orient="records"),
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in input_paths.items()},
        "outputs": [
            "t2d_step8e_axis_evidence_profiles.csv",
            "t2d_step8e_final_classification.csv",
            "t2d_step8e_interpretation_risks.csv",
            "STEP8E_T2D_INTEGRATED_EVIDENCE_REPORT.md",
            "STEP8E_MANIFEST.json",
            "run_step08e_t2d_integrated_evidence.py",
        ],
        "interpretation_boundary": "No causal effect, mediation, pathway activation, exposure-induced expression, or universal T2D mechanism is claimed.",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
    }
    (HERE / "STEP8E_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    report = [
        "# Step 8E — T2D integrated evidence synthesis and flagship selection",
        "",
        "- Status: **complete_integrated_evidence_profiles**",
        "- Scope: four frozen Tier A exposure axes only; no upstream result was re-run or changed.",
        "- Method: evidence-profile classification; no opaque 0–100 total score.",
        "",
        "## Final classification",
        "",
        "| Axis | Biomarker(s) | Classification | Pathway concentration | Transcriptomic support |",
        "|---|---|---|---|---|",
    ]
    for row in final_df.itertuples(index=False):
        report.append(f"| {row.cluster_id} | {row.biomarkers} | **{row.final_classification}** | {row.pathway_concentration} | {row.transcriptomic_support} |")
    report.extend([
        "",
        "## Evidence profiles",
        "",
        "| Axis | Epidemiology | Robustness | GeneCards | Pathway | Network | Transcriptomics |",
        "|---|---|---|---|---|---|---|",
    ])
    for row in profiles_df.itertuples(index=False):
        report.append(f"| {row.cluster_id} | {row.epidemiology_strength} | {row.robustness_strength} | {row.genecard_strength} (q={fmt_p(row.genecard_bh_fdr)}) | {row.pathway_concentration} ({row.pathway_significant_terms} terms; {row.pathway_reduced_modules} modules) | {row.network_strength} (O/E={fmt(row.network_observed_expected_edge_ratio)}; empirical P={fmt_p(row.network_empirical_p)}) | {row.transcriptomic_support} ({row.transcriptomic_modules_tested}/{row.transcriptomic_frozen_modules} modules tested) |")
    report.extend([
        "",
        "### Flagship: cluster_5 / URXP02",
        "",
        "The most concentrated profile is cluster_5. Its compact pathway representatives repeatedly identify xenobiotics, xenobiotic metabolism, cytochrome P450 metabolism, and related chemical-response themes. The axis also has significant network enrichment and human transcriptomic support in a tissue-specific rather than universal pattern. This is a flagship positive demonstration of the framework, not a causal T2D mechanism.",
        "",
        "### Supported axes: cluster_6 and cluster_11",
        "",
        "Both axes retain strong epidemiologic, robustness, GeneCards, and network evidence. Their pathway themes and transcriptomic signals are interpretable but broader or more tissue-dependent than cluster_5, so they are retained as supporting discoveries rather than co-equal flagships.",
        "",
        "### Exploratory axis: cluster_8",
        "",
        "Cluster_8 remains a real FDR-supported and network-supported discovery, but its large gene-set input, broad RNA/protein/cell-cycle annotations, and heterogeneous tissue directionality limit a focused biological claim. Exploratory does not mean negative.",
        "",
        "## Answers to the five prespecified questions",
        "",
        "1. **Most concentrated evidence:** cluster_5 / URXP02, because its pathway signal is thematically focused on xenobiotic/CYP biology while retaining the other evidence layers.",
        "2. **Stable but non-flagship:** cluster_6 and cluster_11; both have strong upstream and network support but broader or more tissue-dependent biology.",
        "3. **Effect of tissue specificity:** it does not invalidate the positive branch, but it prevents a claim of one universal T2D transcriptional direction. T2D tissues remain explicitly separated.",
        "4. **Permitted language:** multi-layer biological convergence, network-supported biological convergence, tissue-specific transcriptomic support, and prioritized environmental axis. Not permitted: causal mechanism, mediation, exposure-induced pathway activation, or universal T2D mechanism.",
        "5. **Framework proof-of-concept:** yes, as an outcome-firewalled environmental screening framework that generated a frozen T2D positive branch and then enabled disease-specific prioritization. The biological branch remains hypothesis-generating rather than causal.",
        "",
        "## Interpretation-risk boundary",
        "",
        "Risk flags are provided in `t2d_step8e_interpretation_risks.csv`. The main recurring boundaries are cross-sectional temporality, reverse-causation potential, annotation-density bias, gene-set size effects, and tissue-specific transcriptomic heterogeneity. cluster_6 remains one biomarker axis despite containing two correlated tests.",
        "",
        "## Reproducibility",
        "",
        "All input paths and SHA-256 checksums are recorded in `STEP8E_MANIFEST.json`. This step reads only frozen Step 5–8D outputs and introduces no new data source.",
    ])
    (HERE / "STEP8E_T2D_INTEGRATED_EVIDENCE_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
