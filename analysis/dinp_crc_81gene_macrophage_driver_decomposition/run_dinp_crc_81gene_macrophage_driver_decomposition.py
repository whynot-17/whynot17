#!/usr/bin/env python
"""Decompose the frozen 81-gene DINP--CRC program in macrophages.

The analysis is donor-aware: macrophage-level expression is first averaged
within donor and disease state, then tumor-minus-normal paired contrasts are
tested across donors.  Cell-level observations are used only to quantify
expression prevalence, not as independent inferential replicates.

Candidate drivers are deliberately named ``network-prioritized`` rather than
causal drivers.  The deterministic screening rule requires a positive paired
effect, paired t-test BH-FDR < 0.05, tumor-cell detection prevalence >= 25%,
and membership in at least one previously observed prostaglandin/arachidonic
acid/inflammatory term.  This is a prioritization step for the next PPI
analysis, not a causal claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "analysis" / "dinp_crc_multi_database_target_convergence" / "outputs"
OUT_DIR = Path(__file__).resolve().parent / "outputs"
DEFAULT_H5AD = Path(
    r"D:\mcop and CRC\cellxgene_census\2025-11-08\16023185-de21-4c0d-a9c8-73abdd52d142.h5ad"
)
INPUT_GENES = INPUT_DIR / "dinp_crc_intersection.csv"
CENSUS_VERSION = "2025-11-08"
TUMOR_LABEL = "colon adenocarcinoma"
NORMAL_LABEL = "normal"
MACROPHAGE_CELL_TYPE = "macrophage"

# Exact memberships returned by g:Profiler for the frozen 81-gene query,
# using the same human annotation release as the preceding enrichment run.
# The full response was requested with no_evidences=False so the query-order
# intersection arrays could be mapped back to symbols and retained below.
PATHWAY_MEMBERS = {
    "GO:0006693": {
        "source": "GO:BP",
        "name": "prostaglandin metabolic process",
        "category": "prostaglandin_AA_inflammation",
        "genes": ["AKR1C1", "AKR1C3", "HPGD", "PLA2G4A", "PTGES", "PTGES2", "PTGES3", "PTGS1", "PTGS2", "SIRT1"],
    },
    "GO:0001516": {
        "source": "GO:BP",
        "name": "prostaglandin biosynthetic process",
        "category": "prostaglandin_AA_inflammation",
        "genes": ["PLA2G4A", "PTGES", "PTGES2", "PTGES3", "PTGS1", "PTGS2", "SIRT1"],
    },
    "KEGG:00590": {
        "source": "KEGG",
        "name": "Arachidonic acid metabolism",
        "category": "prostaglandin_AA_inflammation",
        "genes": ["AKR1C3", "HPGD", "PLA2G4A", "PTGES", "PTGES2", "PTGES3", "PTGS1", "PTGS2"],
    },
    "GO:0006954": {
        "source": "GO:BP",
        "name": "inflammatory response",
        "category": "prostaglandin_AA_inflammation",
        "genes": ["AHSG", "CCL20", "CRP", "CXCL13", "CXCR4", "ESR1", "MMP9", "NEAT1", "NR3C1", "PPARA", "PPARD", "PPARG", "PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTGES", "PTGFR", "PTGS2", "PTX3", "RELA", "SIRT2", "STAT3", "TIMP1", "TNC"],
    },
    "GO:0004879": {
        "source": "GO:MF",
        "name": "nuclear receptor activity",
        "category": "nuclear_receptor",
        "genes": ["ESR1", "NR1I2", "NR1I3", "NR3C1", "PGR", "PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "STAT3"],
    },
    "GO:0009410": {
        "source": "GO:BP",
        "name": "response to xenobiotic stimulus",
        "category": "xenobiotic",
        "genes": ["ABCC4", "AKR1C1", "ATG5", "BECN1", "CA9", "CPS1", "CYP2A6", "GSTA2", "MMP2", "NR1I2", "PPARG"],
    },
    "GO:0006805": {
        "source": "GO:BP",
        "name": "xenobiotic metabolic process",
        "category": "xenobiotic",
        "genes": ["ABCC4", "AKR1C1", "CYP2A6", "GSTA2", "NR1I2"],
    },
}
GP_VERSION = "e114_eg62_p19_27110d83"
GP_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


def sha256_file(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def dense_matrix(value: object) -> np.ndarray:
    if sparse.issparse(value):
        return value.toarray().astype(np.float64, copy=False)
    if hasattr(value, "toarray"):
        return value.toarray().astype(np.float64, copy=False)
    return np.asarray(value, dtype=np.float64)


def bh(values: pd.Series) -> pd.Series:
    output = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(values)
    if valid.any():
        output.loc[valid] = multipletests(
            values.loc[valid].to_numpy(dtype=float), method="fdr_bh"
        )[1]
    return output


def pathway_gene_table(genes: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    for term_id, details in PATHWAY_MEMBERS.items():
        members = set(details["genes"])
        for gene in genes:
            if gene in members:
                rows.append(
                    {
                        "gene_symbol": gene,
                        "term_id": term_id,
                        "source": details["source"],
                        "term_name": details["name"],
                        "pathway_category": details["category"],
                    }
                )
    return pd.DataFrame(rows)


def get_pathway_summary(gene: str) -> tuple[str, str, int, int]:
    hits = [
        (term_id, details)
        for term_id, details in PATHWAY_MEMBERS.items()
        if gene in set(details["genes"])
    ]
    categories = sorted({details["category"] for _, details in hits})
    term_ids = sorted(term_id for term_id, _ in hits)
    core_hits = [details for _, details in hits if details["category"] == "prostaglandin_AA_inflammation"]
    return ";".join(categories), ";".join(term_ids), len(hits), len(core_hits)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    args = parser.parse_args()
    h5ad_path = args.h5ad.expanduser().resolve()
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Source H5AD not found: {h5ad_path}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    source = pd.read_csv(INPUT_GENES)
    genes = sorted(set(source["gene_symbol"].dropna().astype(str).str.upper()))
    if len(genes) != 81:
        raise ValueError(f"Expected 81 frozen genes, found {len(genes)}")
    pd.DataFrame({"gene_symbol": genes}).to_csv(OUT_DIR / "input_81_genes.csv", index=False)
    pathway_gene_table(genes).to_csv(OUT_DIR / "pathway_membership_audit.csv", index=False)

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        feature_names = adata.var["feature_name"].astype(str).str.upper()
        missing = sorted(set(genes) - set(feature_names))
        duplicated = feature_names[feature_names.isin(genes)].value_counts()
        if missing or (duplicated > 1).any():
            raise ValueError(
                f"Gene mapping failed; missing={missing}, duplicated={duplicated[duplicated > 1].to_dict()}"
            )
        gene_indices = [int(np.flatnonzero(feature_names.eq(g).to_numpy())[0]) for g in genes]

        obs = adata.obs.copy()
        obs["donor_id"] = obs["donor_id"].astype(str)
        obs["cell_type"] = obs["cell_type"].astype(str)
        obs["disease_label"] = obs["disease"].astype(str)
        eligible = (
            obs["is_primary_data"].eq(True)
            & obs["cell_type"].eq(MACROPHAGE_CELL_TYPE)
            & obs["disease_label"].isin([TUMOR_LABEL, NORMAL_LABEL])
            & ~obs["donor_id"].isin(["", "nan", "None"])
        )
        positions = np.flatnonzero(eligible.to_numpy())
        meta = obs.iloc[positions][["donor_id", "disease_label"]].reset_index(drop=True)
        if len(meta) < 100:
            raise ValueError(f"Too few eligible macrophage cells: {len(meta)}")

        target_values = dense_matrix(adata.X[:, gene_indices]).astype(np.float32, copy=False)
        values = target_values[positions]
    finally:
        source_shape = [int(adata.n_obs), int(adata.n_vars)]
        adata.file.close()

    n_cells = len(meta)
    tumor_mask = meta["disease_label"].eq(TUMOR_LABEL).to_numpy()
    normal_mask = meta["disease_label"].eq(NORMAL_LABEL).to_numpy()
    donor_ids = sorted(set(meta["donor_id"]))
    donor_mean_rows: list[dict] = []
    donor_delta_rows: list[dict] = []
    stat_rows: list[dict] = []

    # Donor-level means and paired deltas for every frozen gene.
    donor_lookup: dict[str, dict[str, np.ndarray]] = {}
    for donor_id in donor_ids:
        dmask = meta["donor_id"].eq(donor_id).to_numpy()
        donor_lookup[donor_id] = {
            "tumor": values[dmask & tumor_mask].mean(axis=0) if np.any(dmask & tumor_mask) else np.full(len(genes), np.nan),
            "normal": values[dmask & normal_mask].mean(axis=0) if np.any(dmask & normal_mask) else np.full(len(genes), np.nan),
        }
        for group, array in donor_lookup[donor_id].items():
            if np.all(np.isnan(array)):
                continue
            record = {"donor_id": donor_id, "group": group}
            record.update({gene: float(value) for gene, value in zip(genes, array)})
            donor_mean_rows.append(record)

    donor_mean_df = pd.DataFrame(donor_mean_rows)
    donor_mean_df.to_csv(OUT_DIR / "macrophage_donor_gene_means.csv", index=False)
    tumor_donor = donor_mean_df.loc[donor_mean_df["group"].eq("tumor")].set_index("donor_id")[genes]
    normal_donor = donor_mean_df.loc[donor_mean_df["group"].eq("normal")].set_index("donor_id")[genes]
    paired_donor_ids = tumor_donor.index.intersection(normal_donor.index)
    paired_donor_ids = paired_donor_ids[
        tumor_donor.loc[paired_donor_ids].notna().any(axis=1)
        & normal_donor.loc[paired_donor_ids].notna().any(axis=1)
    ]

    for gene in genes:
        tvals = tumor_donor.loc[paired_donor_ids, gene].to_numpy(dtype=float)
        nvals = normal_donor.loc[paired_donor_ids, gene].to_numpy(dtype=float)
        deltas = tvals - nvals
        valid = np.isfinite(deltas)
        deltas = deltas[valid]
        tvals = tvals[valid]
        nvals = nvals[valid]
        for donor_id, tumor_value, normal_value in zip(paired_donor_ids[valid], tvals, nvals):
            donor_delta_rows.append(
                {
                    "donor_id": str(donor_id),
                    "gene_symbol": gene,
                    "tumor_mean": float(tumor_value),
                    "normal_mean": float(normal_value),
                    "tumor_minus_normal": float(tumor_value - normal_value),
                }
            )
        n_pair = int(len(deltas))
        mean_delta = float(np.mean(deltas)) if n_pair else np.nan
        median_delta = float(np.median(deltas)) if n_pair else np.nan
        sd_delta = float(np.std(deltas, ddof=1)) if n_pair > 1 else np.nan
        cohen_dz = mean_delta / sd_delta if np.isfinite(sd_delta) and sd_delta > 0 else np.nan
        if n_pair > 1 and np.isfinite(sd_delta) and sd_delta > 0:
            t_stat, p_t = stats.ttest_1samp(deltas, 0.0)
            ci_half = stats.t.ppf(0.975, n_pair - 1) * sd_delta / math.sqrt(n_pair)
            ci_low, ci_high = mean_delta - ci_half, mean_delta + ci_half
        else:
            t_stat, p_t, ci_low, ci_high = np.nan, np.nan, np.nan, np.nan
        if n_pair >= 5 and np.any(np.abs(deltas) > 0):
            try:
                w_stat, p_w = stats.wilcoxon(deltas, alternative="two-sided", method="auto")
            except ValueError:
                w_stat, p_w = np.nan, np.nan
        else:
            w_stat, p_w = np.nan, np.nan

        tumor_gene = values[tumor_mask, genes.index(gene)]
        normal_gene = values[normal_mask, genes.index(gene)]
        tumor_prev = float(np.mean(tumor_gene != 0)) if len(tumor_gene) else np.nan
        normal_prev = float(np.mean(normal_gene != 0)) if len(normal_gene) else np.nan
        categories, term_ids, term_count, core_term_count = get_pathway_summary(gene)
        stat_rows.append(
            {
                "gene_symbol": gene,
                "n_paired_donors": n_pair,
                "mean_tumor_expression": float(np.mean(tvals)) if len(tvals) else np.nan,
                "mean_normal_expression": float(np.mean(nvals)) if len(nvals) else np.nan,
                "mean_delta_tumor_minus_normal": mean_delta,
                "median_delta_tumor_minus_normal": median_delta,
                "sd_delta": sd_delta,
                "cohen_dz_paired": cohen_dz,
                "mean_delta_95ci_low": ci_low,
                "mean_delta_95ci_high": ci_high,
                "paired_t_statistic": float(t_stat) if np.isfinite(t_stat) else np.nan,
                "paired_t_p": float(p_t) if np.isfinite(p_t) else np.nan,
                "paired_wilcoxon_statistic": float(w_stat) if np.isfinite(w_stat) else np.nan,
                "paired_wilcoxon_p": float(p_w) if np.isfinite(p_w) else np.nan,
                "tumor_cell_detection_fraction": tumor_prev,
                "normal_cell_detection_fraction": normal_prev,
                "detection_fraction_delta": tumor_prev - normal_prev,
                "pathway_categories": categories,
                "pathway_term_ids": term_ids,
                "pathway_term_count": term_count,
                "core_prostaglandin_AA_inflammatory_term_count": core_term_count,
                "pathway_supported": bool(core_term_count > 0),
            }
        )

    donor_delta_df = pd.DataFrame(donor_delta_rows)
    donor_delta_df.to_csv(OUT_DIR / "macrophage_donor_gene_deltas.csv", index=False)
    stats_df = pd.DataFrame(stat_rows)
    stats_df["paired_t_BH_FDR"] = bh(stats_df["paired_t_p"])
    stats_df["paired_wilcoxon_BH_FDR"] = bh(stats_df["paired_wilcoxon_p"])
    stats_df["primary_statistical_eligibility"] = (
        stats_df["mean_delta_tumor_minus_normal"].gt(0)
        & stats_df["paired_t_BH_FDR"].lt(0.05)
        & stats_df["tumor_cell_detection_fraction"].ge(0.25)
    )
    stats_df["network_priority_eligibility"] = (
        stats_df["primary_statistical_eligibility"] & stats_df["pathway_supported"]
    )
    stats_df = stats_df.sort_values(
        ["network_priority_eligibility", "primary_statistical_eligibility", "cohen_dz_paired", "paired_t_BH_FDR"],
        ascending=[False, False, False, True],
        na_position="last",
    ).reset_index(drop=True)
    stats_df["driver_priority_rank"] = np.nan
    eligible_idx = stats_df.index[stats_df["network_priority_eligibility"]]
    stats_df.loc[eligible_idx, "driver_priority_rank"] = np.arange(1, len(eligible_idx) + 1)
    stats_df.to_csv(OUT_DIR / "macrophage_gene_driver_stats.csv", index=False)

    candidates = stats_df.loc[stats_df["network_priority_eligibility"]].copy()
    candidates = candidates.sort_values("driver_priority_rank")
    candidates["candidate_set"] = np.where(candidates["driver_priority_rank"].le(10), "top10_network_priority", "additional_eligible")
    candidates.to_csv(OUT_DIR / "macrophage_driver_candidates.csv", index=False)
    candidates.loc[candidates["driver_priority_rank"].le(10), ["gene_symbol", "driver_priority_rank", "cohen_dz_paired", "paired_t_BH_FDR", "tumor_cell_detection_fraction", "pathway_categories", "pathway_term_ids"]].to_csv(
        OUT_DIR / "macrophage_top10_network_priority_genes.csv", index=False
    )

    source_gene_hash = hashlib.sha256("\n".join(genes).encode()).hexdigest()
    membership_manifest = {
        "source": "g:Profiler human GO/KEGG annotation",
        "endpoint": GP_ENDPOINT,
        "release": GP_VERSION,
        "query_gene_count": len(genes),
        "query_gene_sha256": source_gene_hash,
        "membership_rule": "gene is mapped to the query-order intersection for the exact term ID; no term was added because of macrophage results",
        "terms": [
            {"term_id": term_id, "source": details["source"], "term_name": details["name"], "category": details["category"], "genes": details["genes"]}
            for term_id, details in PATHWAY_MEMBERS.items()
        ],
    }
    (OUT_DIR / "pathway_membership_provenance.json").write_text(
        json.dumps(membership_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    manifest = {
        "analysis": "Frozen 81-gene DINP-CRC macrophage driver decomposition",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_gene_file": str(INPUT_GENES),
        "input_gene_count": len(genes),
        "input_gene_sha256": source_gene_hash,
        "source_h5ad": str(h5ad_path),
        "source_h5ad_sha256": sha256_file(h5ad_path),
        "source_h5ad_bytes": int(h5ad_path.stat().st_size),
        "census_version": CENSUS_VERSION,
        "source_shape": source_shape,
        "expression_matrix": "adata.X (non-raw source matrix); adata.raw was not used",
        "primary_filter": "is_primary_data == True; cell_type == macrophage; disease in {colon adenocarcinoma, normal}; donor_id present",
        "eligible_macrophage_cells": n_cells,
        "paired_donor_count": int(len(paired_donor_ids)),
        "unit_of_inference": "donor-level macrophage gene means; paired tumor-minus-normal contrasts",
        "cell_level_use": "cell-level values used only for detection prevalence; not used as independent inferential replicates",
        "multiple_testing": "BH-FDR across 81 paired gene-level t-tests; Wilcoxon FDR reported separately",
        "primary_driver_rule": "positive paired mean delta AND paired_t_BH_FDR < 0.05 AND tumor detection fraction >= 0.25",
        "pathway_gate": "at least one exact prior g:Profiler term in prostaglandin/arachidonic-acid/inflammatory category",
        "candidate_label_boundary": "network-prioritized candidate; not a causal driver",
        "pathway_membership_source": membership_manifest,
        "outputs": [str(p) for p in sorted(OUT_DIR.glob("*.csv"))],
        "counts": {
            "primary_statistical_eligible": int(stats_df["primary_statistical_eligibility"].sum()),
            "network_priority_eligible": int(stats_df["network_priority_eligibility"].sum()),
            "top10_candidates": int(candidates["driver_priority_rank"].le(10).sum()),
        },
    }
    (OUT_DIR / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    top = candidates.loc[candidates["driver_priority_rank"].le(10)]
    report_lines = [
        "# Frozen 81-gene DINP–CRC program: macrophage driver decomposition",
        "",
        f"Generated: {manifest['run_timestamp_utc']}",
        "",
        "## Analysis boundary",
        "",
        f"- Source: Census release `{CENSUS_VERSION}` H5AD; source file is not copied into the repository.",
        f"- Macrophage cells: `{n_cells:,}`; paired donors: `{len(paired_donor_ids)}`; frozen genes: `{len(genes)}`.",
        "- Expression was summarized within donor and disease state before inference. Cells were not treated as independent replicates.",
        "- `adata.X` was used; `adata.raw` was not used.",
        "",
        "## Frozen prioritization rule",
        "",
        "A gene is primary-statistically eligible when the paired tumor-minus-normal mean is positive, paired t-test BH-FDR is <0.05 across the 81-gene family, and tumor-cell detection fraction is >=25%. A gene is network-priority eligible only if it also belongs to at least one exact, previously observed prostaglandin/arachidonic-acid/inflammatory term. The top 10 are ranked by paired Cohen’s dz, then FDR; these are candidates for network analysis, not causal drivers.",
        "",
        f"- Primary-statistical eligible: `{int(stats_df['primary_statistical_eligibility'].sum())}` genes.",
        f"- Network-priority eligible: `{int(stats_df['network_priority_eligibility'].sum())}` genes.",
        f"- Top network-priority set: `{len(top)}` genes.",
        "",
        "## Top network-prioritized macrophage genes",
        "",
        "| Rank | Gene | Mean Δ | Cohen dz | t-test BH-FDR | Tumor detection | Prior pathway lens |",
        "|---:|---|---:|---:|---:|---:|---|",
    ]
    for _, row in top.iterrows():
        report_lines.append(
            f"| {int(row['driver_priority_rank'])} | **{row['gene_symbol']}** | {row['mean_delta_tumor_minus_normal']:.3f} | {row['cohen_dz_paired']:.2f} | {row['paired_t_BH_FDR']:.3g} | {row['tumor_cell_detection_fraction']:.1%} | {row['pathway_categories']} |"
        )
    report_lines += [
        "",
        "## Interpretation",
        "",
        "This decomposition indicates which members of the frozen 81-gene program are most suitable for the next macrophage-focused PPI/network step. It does not establish that any gene is caused by DINP exposure, that it mediates the epidemiologic association, or that expression change is specific to malignant biology.",
        "",
        "The full 81-gene statistics, donor-level means/deltas, exact pathway-membership audit, source hash, and deterministic candidate set are retained in `outputs/`.",
    ]
    (OUT_DIR / "report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "eligible_macrophage_cells": n_cells,
                "paired_donors": int(len(paired_donor_ids)),
                "primary_statistical_eligible": int(stats_df["primary_statistical_eligibility"].sum()),
                "network_priority_eligible": int(stats_df["network_priority_eligibility"].sum()),
                "top10": top[["gene_symbol", "driver_priority_rank", "mean_delta_tumor_minus_normal", "cohen_dz_paired", "paired_t_BH_FDR"]].to_dict("records"),
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
