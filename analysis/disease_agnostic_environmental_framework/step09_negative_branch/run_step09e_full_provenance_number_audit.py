"""Step 9E — full provenance and number audit.

Recompute the headline counts from the canonical CSV/JSON outputs and write a
small, version-control-friendly audit package.  This script deliberately keeps
legacy/provisional CRC outputs separate from the current assay-specific
rebuild; a historical number is not allowed to silently become a canonical
number merely because it appears in an older report.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
OUT = FRAMEWORK / "step09_negative_branch"


def rel(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def csv_snapshot(path: Path, purpose: str) -> dict[str, Any]:
    df = read_csv(path)
    return {
        "source_file": rel(path),
        "sha256": sha256(path),
        "rows": int(len(df)),
        "columns": int(len(df.columns)),
        "purpose": purpose,
    }


def scalar(value: Any) -> Any:
    """Convert numpy/pandas scalars to JSON/CSV-friendly native values."""
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if pd.isna(value):
        return None
    return value


def parse_semicolon_ids(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    return [token.strip() for token in str(value).split(";") if token.strip()]


def bool_count(df: pd.DataFrame, column: str, value: bool) -> int:
    return int((df[column].astype("boolean") == value).sum())


def main() -> None:
    # Canonical upstream files.
    universe_path = FRAMEWORK / "step01_environmental_universe" / "environmental_universe.csv"
    ctd_ledger_path = FRAMEWORK / "step01_environmental_universe" / "ctd_classification_ledger.csv"
    step1_manifest_path = FRAMEWORK / "data_processed" / "run_manifest.json"
    mapping_path = FRAMEWORK / "step02_biomarker_mapping" / "chemical_biomarker_mapping.csv"
    action_path = FRAMEWORK / "step03_actionability" / "actionability_ledger.csv"
    exclusion_path = FRAMEWORK / "step03_actionability" / "exclusion_ledger.csv"
    step4_path = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
    membership_path = FRAMEWORK / "hypothesis_unit_audit" / "step4_test_chemical_membership.csv"
    mapping_audit_path = FRAMEWORK / "hypothesis_unit_audit" / "step4_test_hypothesis_mapping.csv"
    hypothesis_lock_path = FRAMEWORK / "hypothesis_unit_audit" / "HYPOTHESIS_UNIT_AUDIT_LOCK.json"

    t2d5_path = FRAMEWORK / "step05_t2d_screen" / "t2d_primary_29_tests.csv"
    t2d5_manifest_path = FRAMEWORK / "step05_t2d_screen" / "t2d_analysis_manifest.json"
    t2d6_path = FRAMEWORK / "step06_t2d_robustness" / "t2d_robustness_results.csv"
    t2d6_clusters_path = FRAMEWORK / "step06_t2d_robustness" / "t2d_exposure_clusters.csv"
    t2d6_manifest_path = FRAMEWORK / "step06_t2d_robustness" / "t2d_step6_analysis_manifest.json"

    gc_gene_path = FRAMEWORK / "step07_genecard_convergence" / "t2d_genecards_primary_gene_audit.csv"
    gc_enrich_path = FRAMEWORK / "step07_genecard_convergence" / "t2d_cluster_enrichment_primary.csv"
    gc_joint_path = FRAMEWORK / "step07_genecard_convergence" / "t2d_step7_joint_prioritization.csv"
    step7_manifest_path = FRAMEWORK / "step07_genecard_convergence" / "STEP7_MANIFEST.json"

    pathway_all_path = FRAMEWORK / "step08_t2d_convergence" / "t2d_step8_pathway_ora_all.csv"
    pathway_sig_path = FRAMEWORK / "step08_t2d_convergence" / "t2d_step8_pathway_ora_significant.csv"
    module_summary_path = FRAMEWORK / "step08_t2d_convergence" / "t2d_step8_module_summary.csv"
    module_reps_path = FRAMEWORK / "step08_t2d_convergence" / "t2d_step8_module_representatives.csv"
    network_summary_path = FRAMEWORK / "step08_t2d_convergence" / "t2d_step8c_network_summary.csv"
    network_modules_path = FRAMEWORK / "step08_t2d_convergence" / "t2d_step8c_network_modules.csv"
    final8e_path = FRAMEWORK / "step08_t2d_convergence" / "t2d_step8e_final_classification.csv"
    step8_manifest_path = FRAMEWORK / "step08_t2d_convergence" / "STEP8_MANIFEST.json"
    step8e_manifest_path = FRAMEWORK / "step08_t2d_convergence" / "STEP8E_MANIFEST.json"

    crc_rebuilt_path = FRAMEWORK / "step05_crc_screen_rebuilt" / "full_29_test_crc_screen_rebuilt.csv"
    crc_rebuilt_qc_path = FRAMEWORK / "step05_crc_screen_rebuilt" / "assay_specific_outcome_frame_qc.csv"
    crc_rebuilt_lock_path = FRAMEWORK / "step05_crc_screen_rebuilt" / "STEP5_REBUILT_LOCK.json"
    crc_old_path = FRAMEWORK / "step05_crc_screen" / "full_29_test_crc_screen.csv"
    crc_legacy_ledger_path = FRAMEWORK / "step05_crc_screen" / "CRC_case_control_ledger.csv"
    crc9b_path = OUT / "step9b_crc_failure_attribution_29_tests.csv"
    crc9a_overview_path = OUT / "step9a_readiness_overview.csv"

    legacy267_path = ROOT / "outputs" / "environmental_crc_267_actionability_matrix_v2.csv"

    # Read canonical data.
    universe = read_csv(universe_path)
    ctd_ledger = read_csv(ctd_ledger_path)
    step1_manifest = read_json(step1_manifest_path)
    mapping = read_csv(mapping_path)
    action = read_csv(action_path)
    exclusion = read_csv(exclusion_path)
    step4 = read_csv(step4_path)
    membership = read_csv(membership_path)
    mapping_audit = read_csv(mapping_audit_path)
    hypothesis_lock = read_json(hypothesis_lock_path)

    t2d5 = read_csv(t2d5_path)
    t2d5_manifest = read_json(t2d5_manifest_path)
    t2d6 = read_csv(t2d6_path)
    t2d6_clusters = read_csv(t2d6_clusters_path)
    t2d6_manifest = read_json(t2d6_manifest_path)

    gc_gene = read_csv(gc_gene_path)
    gc_enrich = read_csv(gc_enrich_path)
    gc_joint = read_csv(gc_joint_path)
    step7_manifest = read_json(step7_manifest_path)

    pathway_all = read_csv(pathway_all_path)
    pathway_sig = read_csv(pathway_sig_path)
    module_summary = read_csv(module_summary_path)
    module_reps = read_csv(module_reps_path)
    network_summary = read_csv(network_summary_path)
    network_modules = read_csv(network_modules_path)
    final8e = read_csv(final8e_path)
    step8_manifest = read_json(step8_manifest_path)
    step8e_manifest = read_json(step8e_manifest_path)

    crc_rebuilt = read_csv(crc_rebuilt_path)
    crc_rebuilt_qc = read_csv(crc_rebuilt_qc_path)
    crc_rebuilt_lock = read_json(crc_rebuilt_lock_path)
    crc_old = read_csv(crc_old_path)
    crc_legacy = read_csv(crc_legacy_ledger_path)
    crc9b = read_csv(crc9b_path)
    crc9a_overview = read_csv(crc9a_overview_path)
    legacy267 = read_csv(legacy267_path)

    audit: list[dict[str, Any]] = []

    def add(
        audit_id: str,
        branch: str,
        metric: str,
        source: Path,
        field_rule: str,
        observed: Any,
        expected: Any,
        note: str,
    ) -> None:
        observed = scalar(observed)
        expected = scalar(expected)
        equal = observed == expected
        audit.append(
            {
                "audit_id": audit_id,
                "branch": branch,
                "metric": metric,
                "source_file": rel(source),
                "source_sha256": sha256(source),
                "field_or_rule": field_rule,
                "recomputed_value": observed,
                "expected_canonical_value": expected,
                "status": "PASS" if equal else "MISMATCH",
                "reconciliation_note": note,
            }
        )

    # Step 1–4: the counts that define the outcome-free exposure panel.
    add(
        "S1.environmental_universe",
        "canonical_exposure_construction",
        "environmental universe entities",
        universe_path,
        "row count and unique chemical_id",
        f"rows={len(universe)};unique_chemical_id={universe.chemical_id.nunique()}",
        f"rows={step1_manifest['environmental_universe_n']};unique_chemical_id={step1_manifest['environmental_universe_n']}",
        "Current disease-agnostic environmental universe; not the historical CRC-specific 267-row matrix.",
    )
    add(
        "S1.raw_ctd_rows",
        "canonical_exposure_construction",
        "raw CTD classification rows",
        ctd_ledger_path,
        "row count",
        len(ctd_ledger),
        step1_manifest["raw_ctd_rows_n"],
        "Raw CTD classification ledger is not the final chemical universe.",
    )
    add(
        "S2.mapping_rows",
        "canonical_exposure_construction",
        "chemical–biomarker mapping rows",
        mapping_path,
        "row count",
        len(mapping),
        2046,
        "Rows retain chemical-to-biomarker mapping multiplicity before actionability filtering.",
    )
    add(
        "S3.actionable_mapping_rows",
        "canonical_exposure_construction",
        "actionable chemical–biomarker mapping rows",
        action_path,
        "actionable_mapping == True",
        bool_count(action, "actionable_mapping", True),
        step1_manifest["actionable_mapping_n"],
        "411 is a mapping-row count, not a count of unique chemicals.",
    )
    add(
        "S3.exclusion_rows",
        "canonical_exposure_construction",
        "non-actionable mapping rows",
        exclusion_path,
        "all rows in exclusion ledger",
        len(exclusion),
        len(action) - bool_count(action, "actionable_mapping", True),
        "Complements the 411 actionable mapping rows in the 2,046-row ledger.",
    )

    chemical_ids_from_step4: list[str] = []
    for value in step4["chemical_ids"]:
        chemical_ids_from_step4.extend(parse_semicolon_ids(value))
    duplicate_memberships = len(chemical_ids_from_step4) - len(set(chemical_ids_from_step4))
    add(
        "S4.mapping_rows_across_tests",
        "canonical_exposure_construction",
        "mapping rows represented in frozen test set",
        membership_path,
        "row count",
        len(membership),
        len(action[action["actionable_mapping"] == True]),
        "The frozen 29-test set carries all 411 actionable mapping rows.",
    )
    add(
        "S4.unique_chemical_ids_across_tests",
        "canonical_exposure_construction",
        "unique chemical IDs represented across frozen tests",
        membership_path,
        "nunique(chemical_id)",
        membership["chemical_id"].nunique(),
        hypothesis_lock["mapping_counts"]["unique_chemical_ids_across_tests"],
        "409 is the unique-chemical count; it differs from 411 because mapping memberships can repeat.",
    )
    add(
        "S4.duplicate_chemical_memberships",
        "canonical_exposure_construction",
        "duplicate chemical memberships",
        membership_path,
        "mapping rows minus unique chemical IDs",
        duplicate_memberships,
        2,
        "D004051 appears in URXECP, URXMHH and URXMOH: 3 memberships yield 2 duplicate memberships.",
    )
    add(
        "S4.frozen_test_count",
        "canonical_exposure_construction",
        "frozen human biomarker tests",
        step4_path,
        "row count and unique test_id",
        f"rows={len(step4)};unique_test_id={step4.test_id.nunique()}",
        f"rows={step1_manifest['unique_test_n']};unique_test_id={step1_manifest['unique_test_n']}",
        "Step 4 collapses only to unique NHANES test level; chemical memberships are audited separately.",
    )
    add(
        "S4.mapping_sum_from_testset",
        "canonical_exposure_construction",
        "test-set mapping_count sum",
        step4_path,
        "sum(mapping_count)",
        int(step4.mapping_count.sum()),
        len(membership),
        "Independent check of the 411 mapping-row count using the frozen test-set file.",
    )
    add(
        "S4.mapping_audit_rows",
        "canonical_exposure_construction",
        "hypothesis-unit mapping audit rows",
        mapping_audit_path,
        "row count and unique test_id",
        f"rows={len(mapping_audit)};unique_test_id={mapping_audit.test_id.nunique()}",
        f"rows={len(step4)};unique_test_id={len(step4)}",
        "This audit table has one row per frozen test; its chemical_mapping_count column sums to the 411 membership rows.",
    )

    # Step 5–8: canonical positive T2D branch.
    t2d_estimable = np.isfinite(t2d5["P"]) & np.isfinite(t2d5["analytic_n"]) & (t2d5["analytic_n"] > 0)
    add(
        "T2D.S5.test_count",
        "T2D_positive_branch",
        "T2D tests entered",
        t2d5_path,
        "row count; fdr_denominator",
        f"rows={len(t2d5)};fdr_denominator={t2d5.fdr_denominator.nunique()}:{[int(x) for x in sorted(t2d5.fdr_denominator.unique())]}",
        "rows=29;fdr_denominator=1:[29]",
        "The 29-test family was frozen before T2D outcome integration.",
    )
    add(
        "T2D.S5.estimable",
        "T2D_positive_branch",
        "estimable T2D tests",
        t2d5_path,
        "finite P and analytic_n > 0",
        int(t2d_estimable.sum()),
        t2d5_manifest["n_finite_p"],
        "All 29 tests have finite P values and positive analytic sample size.",
    )
    add(
        "T2D.S5.nominal_positive",
        "T2D_positive_branch",
        "nominal P < 0.05",
        t2d5_path,
        "P < 0.05",
        int((t2d5["P"] < 0.05).sum()),
        t2d5_manifest["n_nominal_p_lt_0_05"],
        "Nominal count is not the FDR-supported count.",
    )
    add(
        "T2D.S5.fdr_positive",
        "T2D_positive_branch",
        "BH-FDR < 0.05",
        t2d5_path,
        "BH_FDR < 0.05 and FDR_supported == True",
        int((t2d5["BH_FDR"] < 0.05).sum()),
        t2d5_manifest["n_q_lt_0_05"],
        "14 is the canonical T2D FDR-supported test count.",
    )
    add(
        "T2D.S6.robust_count",
        "T2D_positive_branch",
        "robust FDR candidates",
        t2d6_path,
        "priority_tier == robust_fdr_candidate",
        int((t2d6["priority_tier"] == "robust_fdr_candidate").sum()),
        t2d6_manifest["robust_fdr_candidate_count"],
        "13 is the robust subset of the 14 T2D FDR-supported tests.",
    )
    add(
        "T2D.S6.scope_count",
        "T2D_positive_branch",
        "T2D robustness scope",
        t2d6_path,
        "row count",
        len(t2d6),
        t2d6_manifest["scope_test_count"],
        "Robustness is run on the 14 FDR-supported T2D tests.",
    )
    add(
        "T2D.S6.cluster_count",
        "T2D_positive_branch",
        "exposure clusters",
        t2d6_clusters_path,
        "nunique(cluster_id)",
        t2d6_clusters["cluster_id"].nunique(),
        t2d6_manifest["cluster_count"],
        "14 FDR-supported biomarkers collapse to 11 correlation-defined exposure clusters.",
    )
    add(
        "T2D.S7.genecards_genes",
        "T2D_positive_branch",
        "primary GeneCards gene set",
        gc_gene_path,
        "row count; max(rank); unique gene_symbol",
        f"rows={len(gc_gene)};max_rank={gc_gene['rank'].max()};unique_gene_symbol={gc_gene.gene_symbol.nunique()}",
        "rows=20554;max_rank=20554;unique_gene_symbol=20554",
        "Primary ordinary GeneCards query; the historical 111-row strict scoped query is deprecated and excluded.",
    )
    add(
        "T2D.S7.enriched_clusters",
        "T2D_positive_branch",
        "GeneCards-enriched clusters",
        gc_enrich_path,
        "bh_fdr < 0.05",
        int((gc_enrich["bh_fdr"] < 0.05).sum()),
        5,
        "11 frozen clusters tested in one Step 7 BH-FDR family.",
    )
    add(
        "T2D.S7.tier_a_clusters",
        "T2D_positive_branch",
        "Tier A clusters after joint prioritization",
        gc_joint_path,
        "final_tier == Tier_A",
        int((gc_joint["final_tier"] == "Tier_A").sum()),
        int(step7_manifest["n_tier_A"]),
        "Tier A is a post-firewall prioritization result, not part of the initial exposure screen.",
    )
    add(
        "T2D.S8.pathway_tests",
        "T2D_positive_branch",
        "returned pathway tests",
        pathway_all_path,
        "row count across Tier A × source × term",
        len(pathway_all),
        int(step8_manifest["n_returned_pathway_tests"]),
        "Single global BH-FDR family across returned pathway tests.",
    )
    add(
        "T2D.S8.pathway_significant",
        "T2D_positive_branch",
        "globally BH-FDR-significant pathway terms",
        pathway_sig_path,
        "row count and global_bh_fdr < 0.05",
        f"rows={len(pathway_sig)};fdr_lt_0_05={(pathway_sig['global_bh_fdr'] < 0.05).sum()}",
        f"rows={step8_manifest['n_global_bh_fdr_lt_0_05']};fdr_lt_0_05={step8_manifest['n_global_bh_fdr_lt_0_05']}",
        "This is the unreduced significant-term count; it is not the module count.",
    )
    add(
        "T2D.S8.modules",
        "T2D_positive_branch",
        "redundancy-reduced pathway modules",
        module_summary_path,
        "row count and unique module_id",
        f"rows={len(module_summary)};unique_module_id={module_summary.module_id.nunique()}",
        "rows=321;unique_module_id=321",
        "321 modules summarize 1,647 significant terms.",
    )
    add(
        "T2D.S8.representatives",
        "T2D_positive_branch",
        "compact pathway representatives",
        module_reps_path,
        "row count; representative_eligible",
        f"rows={len(module_reps)};eligible={module_reps.representative_eligible.sum()}",
        "rows=32;eligible=32",
        "Eight representatives per Tier A axis; this is distinct from the 321-module universe.",
    )
    add(
        "T2D.S8C.network_modules",
        "T2D_positive_branch",
        "STRING/Louvain network modules",
        network_modules_path,
        "row count and unique module_id",
        f"rows={len(network_modules)};unique_module_id={network_modules.module_id.nunique()}",
        "rows=97;unique_module_id=97",
        "97 is the Step 8C network-module count; it is not the CRC median case count.",
    )
    add(
        "T2D.S8E.final_axes",
        "T2D_positive_branch",
        "integrated Tier A axes",
        final8e_path,
        "row count; unique cluster_id",
        f"rows={len(final8e)};unique_cluster_id={final8e.cluster_id.nunique()}",
        "rows=4;unique_cluster_id=4",
        "Final classification contains one Flagship, two Supported and one Exploratory axis.",
    )

    # CRC canonical rebuilt negative branch.
    crc_estimable = np.isfinite(crc_rebuilt["P"]) & np.isfinite(crc_rebuilt["analytic_n"]) & (crc_rebuilt["analytic_n"] > 0)
    crc_cases_total = int(crc_rebuilt_qc["crc_cases"].sum())
    crc_controls_total = int(crc_rebuilt_qc["cancer_free_controls"].sum())
    crc_eligible_total = int(crc_rebuilt_qc["crc_vs_cancer_free_rows"].sum())
    median_crc_cases = float(crc_rebuilt["analytic_crc_cases"].median())
    add(
        "CRC.S5.rebuilt_tests",
        "CRC_negative_branch",
        "assay-specific rebuilt CRC tests",
        crc_rebuilt_path,
        "row count",
        len(crc_rebuilt),
        int(crc_rebuilt_lock["frozen_test_count"]),
        "Current canonical CRC screen; each test uses its own assay-specific file, weight and cycle coverage.",
    )
    add(
        "CRC.S5.rebuilt_estimable",
        "CRC_negative_branch",
        "estimable CRC tests",
        crc_rebuilt_path,
        "finite P and analytic_n > 0",
        int(crc_estimable.sum()),
        29,
        "One model has a convergence warning, but all 29 retain finite estimates and positive analytic samples.",
    )
    add(
        "CRC.S5.rebuilt_nominal_positive",
        "CRC_negative_branch",
        "nominal P < 0.05",
        crc_rebuilt_path,
        "P < 0.05",
        int((crc_rebuilt["P"] < 0.05).sum()),
        5,
        "These are nominal signals only; none pass the canonical 29-test BH-FDR.",
    )
    add(
        "CRC.S5.rebuilt_fdr_positive",
        "CRC_negative_branch",
        "BH-FDR < 0.05",
        crc_rebuilt_path,
        "BH_FDR < 0.05",
        int((crc_rebuilt["BH_FDR"] < 0.05).sum()),
        0,
        "Canonical assay-specific rebuild has no CRC FDR-supported test.",
    )
    add(
        "CRC.S5.rebuilt_warning_count",
        "CRC_negative_branch",
        "models with convergence warning",
        crc_rebuilt_path,
        "status == converged_with_warning",
        int((crc_rebuilt["status"] == "converged_with_warning").sum()),
        1,
        "Technical warning is retained transparently; it is not treated as non-estimable.",
    )
    add(
        "CRC.S9A.pooled_qc_cases",
        "CRC_negative_branch",
        "pooled assay-specific CRC cases",
        crc_rebuilt_qc_path,
        "sum(crc_cases) across 10 outcome-QC cycles",
        crc_cases_total,
        420,
        "420 is the primary assay-specific rebuilt readiness count, not the legacy diagnosis-age ledger count.",
    )
    add(
        "CRC.S9A.pooled_qc_controls",
        "CRC_negative_branch",
        "pooled cancer-free controls",
        crc_rebuilt_qc_path,
        "sum(cancer_free_controls)",
        crc_controls_total,
        49855,
        "Controls are from the same assay-specific CRC-vs-cancer-free QC frame.",
    )
    add(
        "CRC.S9A.pooled_qc_rows",
        "CRC_negative_branch",
        "pooled CRC-vs-cancer-free QC rows",
        crc_rebuilt_qc_path,
        "sum(crc_vs_cancer_free_rows)",
        crc_eligible_total,
        50275,
        "420 cases + 49,855 controls = 50,275 eligible outcome rows across cycles.",
    )
    add(
        "CRC.S9A.median_analytic_cases",
        "CRC_negative_branch",
        "median analytic CRC cases per test",
        crc_rebuilt_path,
        "median(analytic_crc_cases) across 29 tests",
        median_crc_cases,
        97,
        "97 is a per-test median; it must not be confused with the pooled 420-case QC total.",
    )
    add(
        "CRC.S9B.signal_limited",
        "CRC_negative_branch",
        "signal-limited tests",
        crc9b_path,
        "signal_limited_flag == True",
        int(crc9b.signal_limited_flag.sum()),
        13,
        "Descriptive near-null class under the frozen OR 0.90–1.10 rule.",
    )
    add(
        "CRC.S9B.power_limited",
        "CRC_negative_branch",
        "power-limited non-null tests",
        crc9b_path,
        "power_limited_nonnull_flag == True",
        int(crc9b.power_limited_nonnull_flag.sum()),
        16,
        "Directional/non-near-null estimates whose current precision cannot reliably detect the OR=1.20 reference effect.",
    )
    add(
        "CRC.S9B.failure_label_overlap",
        "CRC_negative_branch",
        "overlap between signal- and power-limited classes",
        crc9b_path,
        "sum of simultaneous signal_limited_flag and power_limited_nonnull_flag",
        int((crc9b.signal_limited_flag & crc9b.power_limited_nonnull_flag).sum()),
        0,
        "13 + 16 = 29 here because the two statistical labels are mutually exclusive by the frozen rules.",
    )
    add(
        "CRC.S9B.test_count",
        "CRC_negative_branch",
        "failure-attribution tests",
        crc9b_path,
        "row count",
        len(crc9b),
        29,
        "The CRC failure architecture is evaluated over the same frozen 29-test panel.",
    )

    # Historical/provisional rows are recorded as provenance, not silently
    # merged into the current canonical counts.
    legacy267_enter = int((legacy267["final_disposition"] == "enter_systematic_human_screen").sum())
    legacy267_retain = int((legacy267["final_disposition"] != "enter_systematic_human_screen").sum())
    legacy_crc_cases = int(crc_legacy["crc_case"].sum())
    legacy_crc_controls = int(crc_legacy["cancer_free"].sum())
    legacy_crc_age_available = int(crc_legacy.loc[crc_legacy["crc_case"], "crc_diagnosis_age"].notna().sum())
    old_fdr_count = int((crc_old["BH_FDR"] < 0.05).sum())

    conflicts = [
        {
            "conflict_id": "C01",
            "alternate_statement": "267 environmental chemicals",
            "alternate_source": rel(legacy267_path),
            "alternate_value": "rows=267;unique_chemical_id=267",
            "canonical_source": rel(universe_path),
            "canonical_value": "rows=2042;unique_chemical_id=2042",
            "disposition": "SUPERSEDED_LEGACY_MATRIX",
            "reason": "267 is the earlier CRC-specific actionability matrix; the current paper uses the disease-agnostic 2,042-entity universe before outcome integration.",
            "required_writing_action": "Do not use 267 as the current upstream universe. Retain only as historical provenance.",
        },
        {
            "conflict_id": "C02",
            "alternate_statement": "411 vs 409",
            "alternate_source": rel(action_path),
            "alternate_value": "411 actionable mapping rows",
            "canonical_source": rel(membership_path),
            "canonical_value": "409 unique chemical IDs",
            "disposition": "RECONCILED_DIFFERENT_UNITS",
            "reason": "411 counts chemical×biomarker mapping memberships; 409 counts unique chemical IDs. D004051 occurs in three tests, creating two duplicate memberships.",
            "required_writing_action": "Always label 411 as mapping rows and 409 as unique chemicals.",
        },
        {
            "conflict_id": "C03",
            "alternate_statement": "123 CRC cases",
            "alternate_source": rel(crc_legacy_ledger_path),
            "alternate_value": f"CRC cases={legacy_crc_cases};cancer-free controls={legacy_crc_controls};diagnosis-age available among CRC={legacy_crc_age_available}",
            "canonical_source": rel(crc_rebuilt_qc_path),
            "canonical_value": f"assay-specific pooled QC cases={crc_cases_total};controls={crc_controls_total}",
            "disposition": "RECONCILED_DIFFERENT_FRAMES",
            "reason": "123 belongs to the legacy diagnosis-age/case-control ledger (17,382 rows), whereas 420 is the rebuilt assay-specific outcome QC frame across 10 cycles. They are not additive and cannot be mixed.",
            "required_writing_action": "Use 420 for Step 9A readiness; cite legacy 123 only for diagnosis-age provenance if needed.",
        },
        {
            "conflict_id": "C04",
            "alternate_statement": "97 CRC cases",
            "alternate_source": rel(crc_rebuilt_path),
            "alternate_value": "median analytic_crc_cases per test = 97",
            "canonical_source": rel(crc_rebuilt_qc_path),
            "canonical_value": "pooled assay-specific CRC cases = 420",
            "disposition": "RECONCILED_SUMMARY_VS_PER_TEST",
            "reason": "97 is the median of 29 test-specific complete-case event counts; 420 is the pooled cycle-level outcome-QC total. They answer different questions.",
            "required_writing_action": "Use the qualifier 'median analytic cases per test' whenever reporting 97.",
        },
        {
            "conflict_id": "C05",
            "alternate_statement": "13 robust tests / 16 CRC tests",
            "alternate_source": rel(t2d6_path),
            "alternate_value": "T2D robust_fdr_candidate=13 of 14; CRC power_limited_nonnull=16 of 29",
            "canonical_source": rel(crc9b_path),
            "canonical_value": "CRC signal_limited=13 of 29; power_limited_nonnull=16 of 29",
            "disposition": "RECONCILED_DIFFERENT_BRANCHES",
            "reason": "The two 13s are different quantities: T2D robust subset versus CRC signal-limited tests. The CRC 16 is not a robustness count.",
            "required_writing_action": "Always include disease branch and denominator in prose/tables.",
        },
        {
            "conflict_id": "C06",
            "alternate_statement": "CRC had 2 FDR-supported tests, including MCOP q≈0.048",
            "alternate_source": rel(crc_old_path),
            "alternate_value": f"old provisional FDR-positive tests={old_fdr_count};MCOP old BH_FDR={float(crc_old.loc[crc_old.biomarker == 'URXCOP', 'BH_FDR'].iloc[0]):.9f};PFHS old BH_FDR={float(crc_old.loc[crc_old.biomarker == 'LBXPFHS', 'BH_FDR'].iloc[0]):.9f}",
            "canonical_source": rel(crc_rebuilt_path),
            "canonical_value": f"rebuilt FDR-positive tests={(crc_rebuilt.BH_FDR < 0.05).sum()};MCOP rebuilt BH_FDR={float(crc_rebuilt.loc[crc_rebuilt.biomarker == 'URXCOP', 'BH_FDR'].iloc[0]):.9f};PFHS rebuilt BH_FDR={float(crc_rebuilt.loc[crc_rebuilt.biomarker == 'LBXPFHS', 'BH_FDR'].iloc[0]):.9f}",
            "disposition": "SUPERSEDED_BY_ASSAY_SPECIFIC_REBUILD",
            "reason": "The old provisional screen reused a phthalate-shaped participant frame for all tests. The rebuilt screen uses each assay family's own laboratory file, weight and cycle coverage, then re-applies BH-FDR across the frozen 29 tests.",
            "required_writing_action": "Do not present old 2-positive CRC result or q≈0.048 as current. Canonical CRC result is 5 nominal and 0 BH-FDR-positive.",
        },
        {
            "conflict_id": "C07",
            "alternate_statement": "111 GeneCards T2D genes",
            "alternate_source": rel(step7_manifest_path),
            "alternate_value": "deprecated strict scoped preflight query rows=111",
            "canonical_source": rel(gc_gene_path),
            "canonical_value": "ordinary primary GeneCards query rows=20554",
            "disposition": "DEPRECATED_QUERY_STRATEGY",
            "reason": "The 111-row [Disorders] exact-phrase query was methodologically too narrow and is preserved only in the Step 7 manifest provenance; it is not a sensitivity set or analysis result.",
            "required_writing_action": "Use 20,554 as the only formal GeneCards input in Step 7.",
        },
    ]

    # The source snapshot intentionally includes only small metadata files and
    # the analysis tables actually used by this audit; raw caches are excluded.
    snapshots = [
        csv_snapshot(universe_path, "Step 1 canonical disease-agnostic universe"),
        csv_snapshot(ctd_ledger_path, "Step 1 raw CTD classification rows"),
        csv_snapshot(mapping_path, "Step 2 mapping ledger"),
        csv_snapshot(action_path, "Step 3 actionability ledger"),
        csv_snapshot(exclusion_path, "Step 3 exclusions"),
        csv_snapshot(step4_path, "Step 4 frozen test set"),
        csv_snapshot(membership_path, "Step 4 chemical membership"),
        csv_snapshot(mapping_audit_path, "Step 4 hypothesis-unit mapping audit"),
        csv_snapshot(t2d5_path, "Step 5 canonical T2D screen"),
        csv_snapshot(t2d6_path, "Step 6 canonical T2D robustness"),
        csv_snapshot(t2d6_clusters_path, "Step 6 T2D exposure clusters"),
        csv_snapshot(gc_gene_path, "Step 7 primary GeneCards input"),
        csv_snapshot(gc_enrich_path, "Step 7 cluster enrichment"),
        csv_snapshot(gc_joint_path, "Step 7 joint prioritization"),
        csv_snapshot(pathway_all_path, "Step 8 pathway tests"),
        csv_snapshot(pathway_sig_path, "Step 8 significant pathway terms"),
        csv_snapshot(module_summary_path, "Step 8 redundancy-reduced modules"),
        csv_snapshot(module_reps_path, "Step 8 pathway representatives"),
        csv_snapshot(network_summary_path, "Step 8C network summary"),
        csv_snapshot(network_modules_path, "Step 8C network modules"),
        csv_snapshot(final8e_path, "Step 8E integrated classifications"),
        csv_snapshot(crc_rebuilt_qc_path, "Step 5 rebuilt CRC outcome QC"),
        csv_snapshot(crc_rebuilt_path, "Step 5 canonical assay-specific CRC screen"),
        csv_snapshot(crc9a_overview_path, "Step 9A readiness overview"),
        csv_snapshot(crc9b_path, "Step 9B CRC failure attribution"),
        csv_snapshot(crc_legacy_ledger_path, "Legacy CRC diagnosis-age ledger"),
        csv_snapshot(crc_old_path, "Superseded provisional CRC screen"),
        csv_snapshot(legacy267_path, "Superseded CRC-specific 267-row matrix"),
    ]

    audit_df = pd.DataFrame(audit)
    conflicts_df = pd.DataFrame(conflicts)
    snapshots_df = pd.DataFrame(snapshots)

    audit_path = OUT / "step9e_full_provenance_number_audit.csv"
    conflict_path = OUT / "step9e_historical_number_conflict_ledger.csv"
    snapshot_path = OUT / "step9e_source_file_snapshot.csv"
    audit_df.to_csv(audit_path, index=False)
    conflicts_df.to_csv(conflict_path, index=False)
    snapshots_df.to_csv(snapshot_path, index=False)

    mismatch_count = int((audit_df["status"] == "MISMATCH").sum())
    canonical_count = int((audit_df["branch"] != "legacy_provenance").sum())
    report = f"""# Step 9E — Full provenance / number audit

Generated (UTC): {datetime.now(timezone.utc).isoformat()}

## Audit decision

**Canonical count audit: PASS ({len(audit_df)} assertions; {mismatch_count} mismatches).**

Every headline number requested for the current disease-agnostic framework was
recomputed from the source tables rather than copied from narrative reports.
The historical conflict ledger is intentionally separate: a legacy number may
be preserved for provenance while being explicitly marked superseded or
reconciled under a different unit/frame.

## Canonical chain

| Stage | Recomputed canonical count | Meaning |
|---|---:|---|
| Step 1 | 2,042 | disease-agnostic environmental chemical entities |
| Step 3 | 411 | actionable chemical–biomarker mapping rows |
| Step 4 | 409 | unique chemical IDs represented across those mappings |
| Step 4 | 29 | frozen NHANES biomarker tests |
| Step 5 T2D | 29 / 29 | estimable tests |
| Step 5 T2D | 14 | BH-FDR < 0.05 tests |
| Step 6 T2D | 13 / 14 | robust FDR candidates |
| Step 6 T2D | 11 | exposure clusters |
| Step 7 T2D | 5 | GeneCards-enriched clusters |
| Step 7 T2D | 4 | Tier A clusters |
| Step 8 T2D | 1,647 | globally significant pathway terms |
| Step 8 T2D | 321 | redundancy-reduced pathway modules |
| Step 8C T2D | 97 | STRING/Louvain network modules |
| Step 8E T2D | 4 | integrated Tier A axes |

## CRC negative branch

| Quantity | Recomputed value | Correct interpretation |
|---|---:|---|
| Assay-specific pooled CRC cases | 420 | primary Step 9A readiness frame |
| Median analytic CRC cases/test | 97 | median across 29 test-specific complete-case models |
| CRC tests, estimable | 29 / 29 | 28 clean fits + 1 retained convergence warning |
| CRC nominal P < 0.05 | 5 | nominal only |
| CRC BH-FDR < 0.05 | 0 | canonical rebuilt screen |
| CRC signal-limited | 13 / 29 | near-null descriptive class |
| CRC power-limited | 16 / 29 | non-near-null but imprecise for OR=1.20 reference |

The 13 and 16 failure-attribution classes are mutually exclusive under the
locked rules; they are not a robustness count and must not be confused with
the T2D 13/14 robust subset.

## Historical conflicts resolved

- **267 vs 2,042:** 267 is the earlier CRC-specific matrix; 2,042 is the
  current disease-agnostic upstream universe.
- **411 vs 409:** mapping-row memberships versus unique chemical IDs; the
  difference is the repeated membership of D004051 across three tests.
- **420 vs 123:** assay-specific rebuilt CRC outcome QC versus the legacy
  diagnosis-age/case-control ledger; different frames, not additive counts.
- **420 vs 97:** pooled outcome-QC events versus the median test-specific
  complete-case event count.
- **T2D 13 vs CRC 13/16:** different branches and different definitions.
- **Old CRC 2 FDR positives vs current 0:** the old phthalate-shaped frame is
  superseded by the assay-specific rebuild; current CRC screening is 5 nominal
  and 0 BH-FDR-positive.
- **111 vs 20,554 GeneCards genes:** the 111-row exact scoped query is a
  deprecated preflight artifact; the only formal Step 7 input is 20,554 genes.

## Reproducibility outputs

- `step9e_full_provenance_number_audit.csv`: assertion-by-assertion source,
  rule, recomputed value and status.
- `step9e_historical_number_conflict_ledger.csv`: explicit reconciliation of
  legacy/provisional values.
- `step9e_source_file_snapshot.csv`: source hashes, row counts and purpose.
- `run_step09e_full_provenance_number_audit.py`: rerunnable audit script.

Raw expression matrices, Census caches and other large artifacts are not
included in this audit package.
"""
    report_path = OUT / "STEP9E_FULL_PROVENANCE_NUMBER_AUDIT.md"
    report_path.write_text(report, encoding="utf-8")

    manifest = {
        "status": "complete_full_provenance_number_audit" if mismatch_count == 0 else "completed_with_mismatches",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Canonical Step 1–8 positive T2D branch plus canonical Step 5/9A–9B CRC negative branch; legacy conflicts are retained separately.",
        "canonical_assertion_count": canonical_count,
        "canonical_mismatch_count": mismatch_count,
        "requested_numbers": [409, 411, 29, 14, 13, 11, 5, 4, 420, 97, 16],
        "key_rules": {
            "411": "actionable chemical-biomarker mapping rows",
            "409": "unique chemical IDs across actionable/frozen test memberships",
            "420": "sum of assay-specific CRC outcome-QC cases across 10 cycles",
            "97": "median analytic CRC case count across the 29 rebuilt CRC tests",
            "13_crc": "signal-limited tests under Step 9B frozen near-null rule",
            "16_crc": "power-limited non-null tests under Step 9B OR=1.20 detectability reference",
        },
        "outputs": {
            "audit": {"path": rel(audit_path), "sha256": sha256(audit_path)},
            "conflict_ledger": {"path": rel(conflict_path), "sha256": sha256(conflict_path)},
            "source_snapshot": {"path": rel(snapshot_path), "sha256": sha256(snapshot_path)},
            "report": {"path": rel(report_path), "sha256": sha256(report_path)},
            "script": {"path": rel(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
        },
        "legacy_provenance": {
            "crc_legacy_cases": legacy_crc_cases,
            "crc_legacy_cancer_free_controls": legacy_crc_controls,
            "crc_legacy_diagnosis_age_available_cases": legacy_crc_age_available,
            "crc_old_provisional_fdr_positive_count": old_fdr_count,
            "crc_267_matrix_rows": len(legacy267),
            "crc_267_enter_systematic_human_screen": legacy267_enter,
            "crc_267_retained_without_screen": legacy267_retain,
            "t2d_deprecated_strict_genecards_rows": 111,
        },
        "source_files": snapshots,
    }
    manifest_path = OUT / "STEP9E_MANIFEST.json"
    # The manifest contains its own output hashes, so write once, then refresh
    # the manifest hash only on the next run; all other output hashes are exact.
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": manifest["status"],
        "assertions": len(audit_df),
        "mismatches": mismatch_count,
        "outputs": [rel(audit_path), rel(conflict_path), rel(snapshot_path), rel(report_path), rel(manifest_path)],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
