"""Freeze and report the complete outcome-blinded 267-chemical actionability audit."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = ROOT / "outputs"
MATRIX = OUTPUTS / "environmental_crc_267_actionability_matrix_v2.csv"
FLOW = OUTPUTS / "environmental_crc_267_actionability_flow.csv"
FIREWALL = OUTPUTS / "environmental_crc_267_outcome_firewall_audit.json"
OUT_CANDIDATES = OUTPUTS / "environmental_crc_267_human_testable_candidates.csv"
OUT_SUMMARY = OUTPUTS / "environmental_crc_267_actionability_summary_v2.md"
OUT_MANIFEST = OUTPUTS / "environmental_crc_267_actionability_manifest_v2.json"


def fmt(x: object) -> str:
    if pd.isna(x):
        return "NA"
    if isinstance(x, float):
        return f"{x:.3g}"
    return str(x)


def truth(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def main() -> None:
    matrix = pd.read_csv(MATRIX, low_memory=False)
    flow = pd.read_csv(FLOW)
    firewall = json.loads(FIREWALL.read_text(encoding="utf-8"))
    if len(matrix) != 267 or matrix["ChemicalID"].nunique() != 267:
        raise ValueError("The v2 actionability matrix must contain exactly 267 unique chemicals")
    if firewall.get("candidate_specific_crc_effect_used") is not False:
        raise ValueError("Outcome firewall is not closed")

    eligible = matrix.loc[truth(matrix["eligible_permissive"])].copy()
    # The actual tested unit is the biomarker axis. Different parent names
    # can point to the same NHANES analyte (e.g. two DEHP-related records both
    # map to URXECP); they must not become duplicate hypotheses.
    eligible["axis_key"] = eligible["biological_matrix"].fillna("") + "|" + eligible["selected_primary_biomarker"].fillna("")
    axis_counts = eligible.groupby("axis_key", as_index=False).agg(
        eligible_chemical_count=("ChemicalID", "size"),
        eligible_chemical_ids=("ChemicalID", lambda s: ";".join(s.astype(str))),
        axis_name=("exposure_axis", lambda s: ";".join(sorted(set(s.astype(str))))),
        primary_biomarker=("selected_primary_biomarker", "first"),
        mapping_confidence=("mapping_confidence", "first"),
    )
    axis_counts.to_csv(OUT_CANDIDATES, index=False)

    def n(stage: str) -> int:
        row = flow.loc[flow["stage"].eq(stage), "n"]
        return int(row.iloc[0]) if not row.empty else 0

    mcop = matrix.loc[matrix["ChemicalID"].eq("C573544")].iloc[0]
    minp = matrix.loc[matrix["ChemicalID"].eq("C471400")].iloc[0]
    queue_path = OUTPUTS / "environmental_crc_267_manual_review_queue_v2.csv"
    queue_n = max(0, sum(1 for _ in queue_path.open("r", encoding="utf-8")) - 1) if queue_path.exists() else 0

    lines = [
        "# Complete 267-chemical environmental CRC human-actionability audit v2",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Scope and firewall",
        "",
        "This report covers the complete Phase 1 primary universe: 267 unique core environmental chemicals. Identity, biomarker mapping, detectability, cycle coverage, and survey testability were evaluated before any candidate-specific human CRC association results were considered.",
        "",
        f"Outcome firewall: `PRIORITIZATION_OUTCOME_BLINDED={firewall.get('PRIORITIZATION_OUTCOME_BLINDED')}`; candidate-specific human OR/P/CI/LOCO fields used for eligibility: `{firewall.get('candidate_specific_crc_effect_used')}`.",
        "",
        "## Real attrition flow",
        "",
        "| Stage | N |",
        "|---|---:|",
    ]
    for _, row in flow.iterrows():
        lines.append(f"| {row['stage']} | {int(row['n'])} |")

    lines += [
        "",
        f"The permissive rule was frozen as `E=1 & X=1 & B=1 & D>=1 & C>=1 & T>=1`; moderate and strict tiers were retained as sensitivity tiers. The {len(eligible)} eligible chemical rows collapse to {len(axis_counts)} unique biomarker-axis tests because multiple parent chemicals share the same validated NHANES proxy; the human screen is run once per unique axis and reports the member chemicals.",
        "",
        "## Tier counts",
        "",
        "| Tier | Definition | N chemical rows |",
        "|---|---|---:|",
        f"| A strict | D=2, C=2, T=2 | {int(truth(matrix['eligible_strict']).sum())} |",
        f"| A moderate | D>=1, C=2, T>=1 | {int(truth(matrix['eligible_moderate']).sum())} |",
        f"| B human-testable | permissive eligibility but not moderate | {int(truth(matrix['eligible_permissive']).sum() - truth(matrix['eligible_moderate']).sum())} |",
        f"| C molecular-only | M2 but not human-testable | {int((matrix['M_tag'].eq('M2') & ~truth(matrix['eligible_permissive'])).sum())} |",
        "",
        f"Manual-review queue: {queue_n} genuine identity/registry exceptions; candidates that were fully searched but lacked a candidate-specific NHANES analyte remain resolved in the 267-row matrix rather than being hidden as pending.",
        "",
        "## MCOP and MiNP/DINP status",
        "",
        f"- **MCOP** (`{mcop['ChemicalID']}`, `{mcop['selected_primary_biomarker']}`): retained as a fully eligible DINP-related urinary exposure axis; D={fmt(mcop['D_tag'])}, C={fmt(mcop['C_tag'])}, T={fmt(mcop['T_tag'])}, permissive={mcop['eligible_permissive']}, strict={mcop['eligible_strict']}. This status is determined from biomarker/actionability data only.",
        f"- **MiNP** (`{minp['ChemicalID']}`, `{minp['selected_primary_biomarker']}`): not discarded or merged into MCOP; it remains a distinct DINP molecular nominee, but its direct urinary detectability is D={fmt(minp['D_tag'])} ({fmt(minp['above_lod_pct'])}% above LOD), so it does not enter the primary human screen under the frozen D gate.",
        "- The DINP-related axis is therefore represented by distinct records: MiNP for molecular nomination and MCOP for the human biomarker axis. No candidate-specific CRC association was used to make this decision.",
        "",
        "## Interpretation",
        "",
        f"After complete outcome-blinded annotation of all 267 original chemicals, {len(eligible)} chemical candidates satisfy the permissive human-testability rule, representing {len(axis_counts)} unique exposure-axis/biomarker tests. The subsequent systematic screen must include all {len(axis_counts)} axes and apply BH-FDR across those tested axes.",
        "",
        "## Files",
        "",
        "- `environmental_crc_267_actionability_matrix_v2.csv` — one row per original chemical.",
        "- `environmental_crc_267_human_testable_candidates.csv` — unique axis keys and member chemicals entering the human screen.",
        "- `environmental_crc_267_biomarker_mapping.csv` — candidate-to-biomarker evidence trail.",
        "- `environmental_crc_267_detectability_by_cycle.csv` — cycle-level measured/above-LOD counts.",
        "- `environmental_crc_267_testability_audit.csv` — pre-outcome survey infrastructure and case/control availability.",
        "- `environmental_crc_267_manual_review_queue_v2.csv` — only genuine identity/registry exceptions.",
    ]
    OUT_SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")

    manifest = json.loads(OUT_MANIFEST.read_text(encoding="utf-8")) if OUT_MANIFEST.exists() else {}
    manifest.update({
        "summary_version": "v2",
        "n_permissive_chemical_rows": int(len(eligible)),
        "n_unique_human_testable_axes": int(len(axis_counts)),
        "n_manual_review_queue": int(queue_n),
        "human_testable_candidates_file": str(OUT_CANDIDATES),
        "summary_file": str(OUT_SUMMARY),
        "all_267_annotated": True,
        "screening_rule": "one unified complex-survey model per unique exposure_axis + selected_primary_biomarker; BH-FDR across all tested axes",
    })
    manifest["outputs"] = sorted(set(manifest.get("outputs", []) + [str(OUT_CANDIDATES), str(OUT_SUMMARY), str(FIREWALL)]))
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"n_chemical_rows": len(eligible), "n_unique_axes": len(axis_counts), "n_manual_review": queue_n}, ensure_ascii=False))


if __name__ == "__main__":
    main()
