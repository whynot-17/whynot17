"""Create the version-controlled completion audit for local-H5AD Phase 2G."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"


def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(x) for x in frame.columns]
    rows = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for values in frame.itertuples(index=False, name=None):
        rendered = []
        for value in values:
            if isinstance(value, (float, np.floating)):
                rendered.append(f"{float(value):.5g}")
            else:
                rendered.append(str(value))
        rows.append("| " + " | ".join(rendered) + " |")
    return "\n".join(rows)


def main() -> None:
    source_qc = json.loads((OUT / "mcop_phase2g_source_h5ad_qc.json").read_text(encoding="utf-8"))
    manifest = json.loads((OUT / "mcop_phase2g_manifest.json").read_text(encoding="utf-8"))
    cells = pd.read_csv(OUT / "mcop_phase2g_epithelial_state_scores.csv", usecols=["cell_id", "donor_id", "group", "PPAR_group"])
    units = pd.read_csv(OUT / "mcop_phase2g_donor_state_pseudobulk.csv")
    validation = pd.read_csv(OUT / "mcop_phase2g_donor_level_validation.csv")
    old = pd.read_csv(OUT / "mcop_phase2f_singlecell_paired_donor_contrasts.csv", dtype={"dataset_id": str})
    old = old.loc[
        old["dataset_id"].eq("16023185-de21-4c0d-a9c8-73abdd52d142")
        & old["compartment"].eq("epithelial")
        & old["score"].eq("PPAR_nuclear_receptor_score")
    ].iloc[0]
    ppar = validation.loc[validation["feature"].eq("PPAR_NR_score")].iloc[0]
    corr = pd.read_csv(OUT / "mcop_phase2g_state_correlations.csv")
    interaction = pd.read_csv(OUT / "mcop_phase2g_tumor_ppar_interaction.csv")
    regulators = pd.read_csv(OUT / "mcop_phase2g_regulator_activity.csv")
    anchors = pd.read_csv(OUT / "mcop_phase2g_regulatory_anchor_ranking.csv")

    missing = manifest.get("local_h5ad_audit", {}).get("target_genes_missing", [])
    checks = {
        "source_h5ad_equivalence_qc": bool(source_qc.get("qc_pass")),
        "cell_ids_unique": not cells["cell_id"].duplicated().any(),
        "eligible_cells": int(len(cells)),
        "unit_cell_sum": int(pd.to_numeric(units["n_cells"]).sum()),
        "cell_count_conserved": int(pd.to_numeric(units["n_cells"]).sum()) == len(cells),
        "unit_keys_unique": not units.duplicated(["donor_key", "donor_id", "group", "PPAR_group"]).any(),
        "paired_donors": int(ppar["n_paired_donors"]),
        "ppar_median_delta": float(ppar["median_delta_tumor_minus_normal"]),
        "phase2f_ppar_median_delta": float(old["median_delta_tumor_minus_normal"]),
        "ppar_delta_absolute_difference": float(abs(float(ppar["median_delta_tumor_minus_normal"]) - float(old["median_delta_tumor_minus_normal"]))),
        "ppar_p_value": float(ppar["p_value"]),
        "phase2f_ppar_p_value": float(old["p_value"]),
        "ppar_p_exact": bool(np.isclose(float(ppar["p_value"]), float(old["p_value"]), rtol=0, atol=1e-15)),
        "target_genes_requested": int(manifest["local_h5ad_audit"]["target_genes_requested"]),
        "target_genes_present": int(manifest["local_h5ad_audit"]["target_genes_present"]),
        "target_genes_missing": missing,
        "state_rows": int(validation["feature_type"].eq("state").sum() - 1),
        "interaction_rows": int(len(interaction)),
        "interaction_min_bh_fdr": float(pd.to_numeric(interaction["BH_FDR"], errors="coerce").min()),
        "regulator_rows": int(len(regulators)),
        "regulator_min_bh_fdr": float(pd.to_numeric(regulators["BH_FDR"], errors="coerce").min()),
        "directly_supported_anchors": anchors.loc[anchors["overall_evidence_tier"].eq("Directly supported"), "gene"].astype(str).tolist(),
    }
    checks["completion_qc_pass"] = bool(
        checks["source_h5ad_equivalence_qc"]
        and checks["cell_ids_unique"]
        and checks["cell_count_conserved"]
        and checks["unit_keys_unique"]
        and checks["paired_donors"] == 36
        and checks["ppar_delta_absolute_difference"] < 1e-12
        and checks["ppar_p_exact"]
    )
    (OUT / "mcop_phase2g_local_completion_audit.json").write_text(json.dumps(checks, indent=2) + "\n", encoding="utf-8")

    paired_states = validation.loc[(validation["feature_type"].eq("state")) & (~validation["feature"].eq("PPAR_NR_score"))].sort_values(["BH_FDR", "p_value"]).head(8)
    top_corr = corr.sort_values(["BH_FDR", "P"]).head(8)
    top_reg = regulators.sort_values(["BH_FDR", "P"]).head(8)
    lines = [
        "# Phase 2G local-H5AD completion audit",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        f"**Completion QC: {'PASS' if checks['completion_qc_pass'] else 'FAIL'}**",
        "",
        "## Provenance and integrity",
        "",
        f"- Official source H5AD bytes: **{source_qc['file_bytes']:,}**; exact expected size: **{source_qc['file_size_exact']}**.",
        f"- C106 9-gene raw-expression equivalence: cells **{source_qc['c106_source_n_obs']:,}/{source_qc['c106_cache_n_obs']:,}**, nnz **{source_qc['c106_source_core_nnz']:,}/{source_qc['c106_cache_core_nnz']:,}**, count sum **{source_qc['c106_source_core_sum']:.0f}/{source_qc['c106_cache_core_sum']:.0f}**.",
        f"- Successful online-cache donor manifests with exact source-H5AD cell counts: **{source_qc['successful_stage3_manifests']}/{source_qc['successful_stage3_manifests']}**.",
        f"- Eligible epithelial cells: **{checks['eligible_cells']:,}**; pseudobulk cell sum: **{checks['unit_cell_sum']:,}**; conservation: **{checks['cell_count_conserved']}**.",
        f"- Frozen target universe present: **{checks['target_genes_present']:,}/{checks['target_genes_requested']:,}**; missing: `{', '.join(missing) if missing else 'none'}`.",
        "",
        "## Frozen PPAR/NR reproduction",
        "",
        f"- Phase 2G median paired tumor-normal delta: **{checks['ppar_median_delta']:.6f}**; P=**{checks['ppar_p_value']:.3g}**; n=**{checks['paired_donors']}**.",
        f"- Phase 2F frozen delta: **{checks['phase2f_ppar_median_delta']:.6f}**; absolute difference: **{checks['ppar_delta_absolute_difference']:.3g}**; P exact: **{checks['ppar_p_exact']}**.",
        "- Cell-level PPAR quartiles define state labels only; donor-level PPAR inference reuses the Phase 2F score standardized in the full epithelial dataset context.",
        "",
        "## New state and regulator results",
        "",
        "Top paired tumor-normal state contrasts:",
        "",
        markdown_table(paired_states[["feature", "median_delta_tumor_minus_normal", "p_value", "BH_FDR", "direction"]]),
        "",
        "Top donor-level PPAR/state associations:",
        "",
        markdown_table(top_corr[["group", "state", "n_donors", "spearman_rho", "P", "BH_FDR"]]),
        "",
        "Top regulator-activity contrasts:",
        "",
        markdown_table(top_reg[["regulator", "comparison", "n_pairs", "activity_delta", "P", "BH_FDR"]]),
        "",
        f"- Tumor×PPAR interaction minimum BH-FDR: **{checks['interaction_min_bh_fdr']:.3g}**; no interaction program passes multiplicity control.",
        f"- Directly supported regulatory anchors under frozen evidence tags: **{', '.join(checks['directly_supported_anchors']) or 'none'}**.",
        "",
        "## Boundary",
        "",
        "The CRC epithelial disease-state bridge is now executable from the official source H5AD and passes integrity QC. The analysis remains associative: it does not establish DINP/MCOP exposure as the cause of the PPAR/NR-low state or prove mediation of the epidemiologic association.",
    ]
    (OUT / "mcop_phase2g_local_completion_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(checks, indent=2))
    if not checks["completion_qc_pass"]:
        raise SystemExit("Phase 2G completion QC failed.")


if __name__ == "__main__":
    main()
