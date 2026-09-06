from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


REPO = Path(r"C:/Users/21634/Documents/Codex/2026-08-22/non/work/whynot17")
M4_DIR = REPO / "analysis/disease_agnostic_environmental_framework/urxp02_m4_sex_tissue_mapping"
M3_DIR = REPO / "analysis/disease_agnostic_environmental_framework/urxp02_m3_disease_branch_analysis"
OUT = REPO / "analysis/disease_agnostic_environmental_framework/urxp02_m4b_hypertension_branch_readout"

FOCUS_TISSUES = [
    "Artery - Aorta",
    "Artery - Tibial",
    "Artery - Coronary",
    "Heart - Left Ventricle",
    "Heart - Atrial Appendage",
    "Kidney - Cortex",
    "Adrenal Gland",
]
GROUPS = {
    "Artery - Aorta": "artery",
    "Artery - Tibial": "artery",
    "Artery - Coronary": "artery",
    "Heart - Left Ventricle": "heart",
    "Heart - Atrial Appendage": "heart",
    "Kidney - Cortex": "kidney",
    "Adrenal Gland": "adrenal",
    "Thyroid": "thyroid_context",
    "Liver": "liver_context",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def yes_count(series: pd.Series) -> int:
    return int(series.astype(str).str.lower().isin(["true", "1", "yes"]).sum())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    m4 = pd.read_csv(M4_DIR / "04_hypertension_branch_sex_expression_by_tissue.csv")
    modules = pd.read_csv(M3_DIR / "05_ppi_modules.csv")
    modules = modules.loc[modules["branch_class"].eq("hypertension-branch")].copy()
    module_meta = modules.drop_duplicates("module_id")[["module_id", "module_size", "module_internal_edges", "module_density"]]
    gene_module = modules[["gene_symbol", "module_id", "module_size"]].drop_duplicates("gene_symbol")

    m4["focus_group"] = m4["tissue"].map(GROUPS)
    m4["focus_relevant"] = m4["tissue"].isin(FOCUS_TISSUES)
    ok = m4["model_status"].eq("OK")
    relevant = m4.loc[ok & m4["focus_relevant"]].copy()
    hits = relevant.loc[relevant["FDR"] < 0.05].copy()
    hits = hits.merge(gene_module, on="gene_symbol", how="left", suffixes=("", "_module"))

    hit_cols = {
        "gene_symbol": "gene_symbol", "gene_id": "gene_id", "branch_membership": "branch_membership",
        "tissue": "tissue", "focus_group": "focus_group", "tissue_role": "tissue_role",
        "n_male": "n_male", "n_female": "n_female", "female_minus_male_beta": "beta_female_minus_male",
        "female_minus_male_se": "se_female_minus_male", "t_statistic": "t_statistic", "raw_p": "raw_p",
        "FDR": "fdr_fixed_3960", "direction": "direction", "meaningful_abs_effect_ge_0_25": "meaningful_abs_effect_ge_0_25",
        "is_expressed_median_tpm_ge_0_1": "is_expressed_median_tpm_ge_0_1",
        "exact_2NAP_human_support": "exact_2NAP_human_support", "exact_2NAP_experimental_support": "exact_2NAP_experimental_support",
        "parent_naphthalene_support": "parent_naphthalene_support", "number_of_sources": "number_of_sources",
        "multi_source_ge2": "multi_source_ge2", "network_central_top10pct_any": "network_central_top10pct_any",
        "m3_priority_hub_candidate": "m3_priority_hub_candidate", "module_id": "module_id", "module_size": "module_size",
    }
    hit_table = hits[list(hit_cols)].rename(columns=hit_cols).sort_values(["tissue", "fdr_fixed_3960", "gene_symbol"])
    hit_table.to_csv(OUT / "01_hypertension_branch_fdr_hits.csv", index=False)

    summary_rows: list[dict[str, object]] = []
    for tissue, sub in m4.groupby("tissue", sort=False):
        tissue_ok = sub.loc[sub["model_status"].eq("OK")]
        tissue_hits = tissue_ok.loc[tissue_ok["FDR"] < 0.05]
        relevant_flag = tissue in FOCUS_TISSUES
        summary_rows.append({
            "summary_level": "tissue",
            "focus_group": GROUPS.get(tissue, "other"),
            "focus_relevant": relevant_flag,
            "tissue": tissue,
            "tissue_role": sub["tissue_role"].iloc[0],
            "n_planned_tests": 440,
            "n_modelled": len(tissue_ok),
            "n_expressed": yes_count(tissue_ok["is_expressed_median_tpm_ge_0_1"]),
            "n_male": int(tissue_ok["n_male"].dropna().iloc[0]) if tissue_ok["n_male"].notna().any() else "NA",
            "n_female": int(tissue_ok["n_female"].dropna().iloc[0]) if tissue_ok["n_female"].notna().any() else "NA",
            "n_fdr_hits": len(tissue_hits),
            "n_female_biased_fdr": int((tissue_hits["direction"] == "female-biased").sum()),
            "n_male_biased_fdr": int((tissue_hits["direction"] == "male-biased").sum()),
            "n_meaningful_abs_effect_ge_0_25": yes_count(tissue_ok["meaningful_abs_effect_ge_0_25"]),
            "median_beta_all_modelled": tissue_ok["female_minus_male_beta"].median(),
            "median_beta_fdr_hits": tissue_hits["female_minus_male_beta"].median() if len(tissue_hits) else "NA",
            "interpretation": "primary hypertension-related tissue" if relevant_flag else "context-only tissue in same fixed family",
        })

    # Fixed descriptive organ grouping; no new hypothesis test or correction.
    for group in ["artery", "heart", "kidney", "adrenal"]:
        sub = relevant.loc[relevant["focus_group"].eq(group)]
        group_hits = sub.loc[sub["FDR"] < 0.05]
        summary_rows.append({
            "summary_level": "organ_group",
            "focus_group": group,
            "focus_relevant": True,
            "tissue": "; ".join(FOCUS_TISSUES if group == "all" else [t for t in FOCUS_TISSUES if GROUPS[t] == group]),
            "tissue_role": "aggregated descriptive readout; not a new test",
            "n_planned_tests": int(sub["tissue"].nunique() * 440),
            "n_modelled": len(sub),
            "n_expressed": yes_count(sub["is_expressed_median_tpm_ge_0_1"]),
            "n_male": "; ".join(sorted({str(int(x)) for x in sub["n_male"].dropna().unique()})),
            "n_female": "; ".join(sorted({str(int(x)) for x in sub["n_female"].dropna().unique()})),
            "n_fdr_hits": len(group_hits),
            "n_female_biased_fdr": int((group_hits["direction"] == "female-biased").sum()),
            "n_male_biased_fdr": int((group_hits["direction"] == "male-biased").sum()),
            "n_meaningful_abs_effect_ge_0_25": yes_count(sub["meaningful_abs_effect_ge_0_25"]),
            "median_beta_all_modelled": sub["female_minus_male_beta"].median(),
            "median_beta_fdr_hits": group_hits["female_minus_male_beta"].median() if len(group_hits) else "NA",
            "interpretation": "descriptive aggregation of fixed tissue-level tests; not a new FDR family",
        })
    pd.DataFrame(summary_rows).to_csv(OUT / "02_hypertension_branch_tissue_summary.csv", index=False)

    recurrence_rows: list[dict[str, object]] = []
    for gene, sub in hits.groupby("gene_symbol", sort=False):
        directions = sorted(set(sub["direction"].dropna()))
        recurrence_rows.append({
            "gene_symbol": gene,
            "n_relevant_tissues_with_fdr": int(sub["tissue"].nunique()),
            "tissues_with_fdr": "; ".join(sorted(sub["tissue"].unique())),
            "directions": "; ".join(directions),
            "direction_consistency": "female_only" if directions == ["female-biased"] else "male_only" if directions == ["male-biased"] else "mixed",
            "female_biased_hit_count": int((sub["direction"] == "female-biased").sum()),
            "male_biased_hit_count": int((sub["direction"] == "male-biased").sum()),
            "min_fdr": float(sub["FDR"].min()),
            "min_raw_p": float(sub["raw_p"].min()),
            "max_abs_beta": float(sub["female_minus_male_beta"].abs().max()),
            "mean_beta": float(sub["female_minus_male_beta"].mean()),
            "repeat_flag": "yes" if sub["tissue"].nunique() >= 2 else "no",
            "module_id": "; ".join(sorted(set(sub["module_id"].dropna()))),
            "exact_2NAP_human_support": bool(sub["exact_2NAP_human_support"].astype(bool).any()),
            "multi_source_ge2": bool(sub["multi_source_ge2"].astype(bool).any()),
            "network_central_top10pct_any": bool(sub["network_central_top10pct_any"].astype(bool).any()),
            "m3_priority_hub_candidate": bool(sub["m3_priority_hub_candidate"].astype(bool).any()),
        })
    pd.DataFrame(recurrence_rows).sort_values(["n_relevant_tissues_with_fdr", "min_fdr", "gene_symbol"], ascending=[False, True, True]).to_csv(OUT / "03_hypertension_branch_gene_recurrence.csv", index=False)

    module_rows: list[dict[str, object]] = []
    for module_id, meta in module_meta.set_index("module_id").iterrows():
        sub = hits.loc[hits["module_id"].eq(module_id)]
        module_rows.append({
            "module_id": module_id,
            "module_size": int(meta["module_size"]),
            "module_internal_edges": int(meta["module_internal_edges"]),
            "module_density": float(meta["module_density"]),
            "n_fdr_hit_tests_relevant": len(sub),
            "n_unique_fdr_hit_genes_relevant": int(sub["gene_symbol"].nunique()),
            "n_tissues_with_fdr_hits": int(sub["tissue"].nunique()),
            "female_biased_fdr_hits": int((sub["direction"] == "female-biased").sum()),
            "male_biased_fdr_hits": int((sub["direction"] == "male-biased").sum()),
            "hit_tissues": "; ".join(sorted(sub["tissue"].unique())) if len(sub) else "",
            "hit_genes": "; ".join(sorted(sub["gene_symbol"].unique())) if len(sub) else "",
            "descriptive_module_status": "observed_fdr_hits" if len(sub) else "no_relevant_fdr_hits",
        })
    pd.DataFrame(module_rows).sort_values(["n_unique_fdr_hit_genes_relevant", "n_fdr_hit_tests_relevant", "module_id"], ascending=[False, False, True]).to_csv(OUT / "04_hypertension_branch_module_context.csv", index=False)

    tissue_summary = pd.DataFrame(summary_rows)
    focus_tissue_summary = tissue_summary.loc[(tissue_summary.summary_level == "tissue") & tissue_summary.focus_relevant].copy()
    top_recurrence = pd.DataFrame(recurrence_rows).sort_values(["n_relevant_tissues_with_fdr", "min_fdr"], ascending=[False, True])
    module_summary = pd.DataFrame(module_rows)
    relevant_hit_n = len(hits)
    female_n = int((hits["direction"] == "female-biased").sum())
    male_n = int((hits["direction"] == "male-biased").sum())
    kidney_row = focus_tissue_summary.loc[focus_tissue_summary.tissue.eq("Kidney - Cortex")].iloc[0]
    top_tissues = focus_tissue_summary.sort_values(["n_fdr_hits", "n_female_biased_fdr"], ascending=False).head(4)
    top_tissue_lines = "\n".join(
        f"- {r.tissue}: {int(r.n_fdr_hits)} FDR hits ({int(r.n_female_biased_fdr)} female-biased, {int(r.n_male_biased_fdr)} male-biased)."
        for r in top_tissues.itertuples()
    )
    recurring_lines = "\n".join(
        f"- {r.gene_symbol}: {int(r.n_relevant_tissues_with_fdr)} tissues — {r.tissues_with_fdr} ({r.direction_consistency})."
        for r in top_recurrence.head(9).itertuples()
    )
    module_lines = "\n".join(
        f"- {r.module_id}: {int(r.n_unique_fdr_hit_genes_relevant)} unique hit genes across {int(r.n_tissues_with_fdr_hits)} tissues; {int(r.female_biased_fdr_hits)} female-biased / {int(r.male_biased_fdr_hits)} male-biased."
        for r in module_summary.loc[module_summary.n_unique_fdr_hit_genes_relevant.gt(0)].head(6).itertuples()
    )
    report = f"""# M4b — hypertension-branch tissue readout

## Scope

This is a readout of the frozen M4 hypertension-branch results only. The fixed family remains 440 genes × 9 prespecified GTEx tissues = 3,960 tests; no model was refit, no new multiple-testing correction was introduced, and no single-cell or mechanism analysis was performed.

## Main result

Across the seven hypertension-relevant tissues (three arteries, two heart tissues, kidney cortex, and adrenal gland), there are **{relevant_hit_n} fixed-FDR hits**: **{female_n} female-biased** and **{male_n} male-biased**. This is a mixed localized pattern, not a coherent female-biased hypertension-branch shift.

{top_tissue_lines}

Kidney cortex has **{int(kidney_row.n_fdr_hits)} fixed-FDR hits**, despite {int(kidney_row.n_meaningful_abs_effect_ge_0_25)} tests meeting the prespecified absolute-effect flag. The kidney GTEx sample size is much smaller than the large artery tissues, so this is a power/context limitation rather than evidence of a kidney null mechanism.

## Recurring genes

The repeated-gene pattern is descriptive only; recurrence was not assigned a new P value:

{recurring_lines}

The clearest repeated female-biased gene is CES1 (adrenal, aorta, tibial). PLA2G5 repeats across adrenal and arterial tissues with a male-biased direction. BCL2 is repeated but directionally mixed across tissues. These observations do not establish disease causality.

## Module context

FDR hits occur in **{int((module_summary.n_unique_fdr_hit_genes_relevant > 0).sum())} of {len(module_summary)} fixed hypertension-branch STRING modules**. The largest descriptive concentrations are:

{module_lines}

This is a module-membership overlay on already-estimated gene×tissue results, not a new module enrichment test.

## Relation to the NHANES female-hypertension phenotype

The relevant-tissue readout provides some female-biased hits in aorta, tibial artery, heart, and adrenal tissue, but female-biased hits are fewer than male-biased hits and kidney has no fixed-FDR hits. Therefore M4b does **not** support a simple, globally female-biased molecular hypertension program. It supports localized and directionally mixed sex-by-tissue differences. The kidney result should be treated as unresolved because of limited GTEx female sample size, motivating targeted—but explicitly exploratory—kidney validation in a donor-balanced single-cell dataset.

## Files and guardrails

- `01_hypertension_branch_fdr_hits.csv`: all fixed-FDR hits in the seven relevant tissues.
- `02_hypertension_branch_tissue_summary.csv`: tissue-level and descriptive organ-group counts, including thyroid/liver context rows.
- `03_hypertension_branch_gene_recurrence.csv`: repeated-gene overlay across relevant tissues.
- `04_hypertension_branch_module_context.csv`: descriptive overlay onto fixed M3 STRING modules.

No figures, no new FDR family, no candidate ranking, and no causal or mechanistic interpretation were added.

Generated UTC: {datetime.now(timezone.utc).isoformat()}
"""
    (OUT / "M4B_HYPERTENSION_BRANCH_READOUT_REPORT.md").write_text(report, encoding="utf-8")

    inputs = [
        M4_DIR / "04_hypertension_branch_sex_expression_by_tissue.csv",
        M4_DIR / "manifest.json",
        M3_DIR / "05_ppi_modules.csv",
    ]
    output_paths = [
        OUT / "01_hypertension_branch_fdr_hits.csv",
        OUT / "02_hypertension_branch_tissue_summary.csv",
        OUT / "03_hypertension_branch_gene_recurrence.csv",
        OUT / "04_hypertension_branch_module_context.csv",
        OUT / "M4B_HYPERTENSION_BRANCH_READOUT_REPORT.md",
    ]
    manifest = {
        "analysis": "URXP02 M4b hypertension-branch tissue readout",
        "scope": "read-only reorganization of frozen M4 outputs; no refit, no new correction, no single-cell analysis",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "fixed_family": {"genes": 440, "tissues": 9, "tests": 3960, "fdr_method": "inherited fixed BH FDR from M4"},
        "focus_tissues": FOCUS_TISSUES,
        "inputs_sha256": {str(p.relative_to(REPO)): sha256(p) for p in inputs},
        "outputs_sha256": {p.name: sha256(p) for p in output_paths},
        "outputs": [
            "01_hypertension_branch_fdr_hits.csv",
            "02_hypertension_branch_tissue_summary.csv",
            "03_hypertension_branch_gene_recurrence.csv",
            "04_hypertension_branch_module_context.csv",
            "M4B_HYPERTENSION_BRANCH_READOUT_REPORT.md",
            "manifest.json",
        ],
        "constraints": [
            "No model refitting",
            "No new FDR or scoring",
            "No single-cell analysis",
            "No figures",
            "No causal or mechanism claims",
        ],
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
