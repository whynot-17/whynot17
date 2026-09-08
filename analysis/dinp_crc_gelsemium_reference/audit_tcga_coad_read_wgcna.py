"""Audit candidate filtering, pairwise co-expression, and WGCNA sensitivity.

This audit does not change the primary WGCNA result. It checks whether the
negative module-enrichment result is plausibly explained by HVG filtering or
by a single network parameterization.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as t_dist


ROOT = Path(__file__).resolve().parents[2]
FRESH = ROOT / "analysis" / "dinp_crc_gelsemium_reference"
OUTPUT = FRESH / "outputs" / "tcga_coad_read_wgcna_audit"
CANONICAL_MANIFEST = FRESH / "outputs" / "tcga_coad_read_wgcna" / "tcga_coad_read_wgcna_sample_manifest.csv"
TARGET_FILE = FRESH / "outputs" / "query_41_genes.csv"
COAD_EXPR = FRESH / "outputs" / "tcga_coad_bulk_41_genes" / "tcga_coad_expression_41_genes.csv"
READ_EXPR = FRESH / "outputs" / "tcga_read_bulk_41_genes" / "tcga_read_expression_41_genes.csv"
FULL_PROBE_STATS = Path("/mnt/d/whynot17/work/tcga_coad_read_wgcna_runtime/tcga_coad_read_wgcna_probe_stats.csv.gz")


def require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def load_targets() -> list[str]:
    frame = pd.read_csv(require(TARGET_FILE))
    targets = frame["gene_symbol"].astype(str).str.strip().drop_duplicates().tolist()
    if len(targets) != 41:
        raise ValueError(f"Expected 41 target genes, got {len(targets)}")
    return targets


def load_expression(targets: list[str]) -> pd.DataFrame:
    coad = pd.read_csv(require(COAD_EXPR), index_col=0)
    read = pd.read_csv(require(READ_EXPR), index_col=0)
    expression = pd.concat([coad, read], axis=0)
    manifest = pd.read_csv(require(CANONICAL_MANIFEST), dtype=str)
    required_columns = {"sample_id", "cohort", "group"}
    if not required_columns.issubset(manifest.columns):
        raise ValueError(f"WGCNA sample manifest lacks required columns: {required_columns - set(manifest.columns)}")
    sample_ids = manifest["sample_id"].astype(str).tolist()
    if len(sample_ids) != 380:
        raise ValueError(f"Expected 380 primary-tumor samples for the audit, got {len(sample_ids)}")
    expression = expression.loc[sample_ids]
    expression = expression.loc[:, targets].astype(float)
    if expression.isna().any().any():
        raise ValueError("Target expression contains missing values")
    return expression


def correlation_audit(expression: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = expression.corr(method="pearson")
    rows = []
    n = len(expression)
    for i, gene_a in enumerate(matrix.columns):
        for gene_b in matrix.columns[i + 1 :]:
            r = float(matrix.loc[gene_a, gene_b])
            if abs(r) >= 1:
                p = 0.0
            else:
                statistic = r * np.sqrt((n - 2) / max(1e-15, 1 - r * r))
                p = float(2 * t_dist.sf(abs(statistic), df=n - 2))
            rows.append({"gene_a": gene_a, "gene_b": gene_b, "pearson_r": r, "abs_r": abs(r), "p_value": p})
    pairs = pd.DataFrame(rows).sort_values(["abs_r", "gene_a", "gene_b"], ascending=[False, True, True])
    pairs["BH_FDR"] = _bh(pairs["p_value"].to_numpy())
    return matrix, pairs


def _bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def load_selection(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(require(path), dtype={"gene_symbol": str})
    return frame.assign(gene_symbol=frame["gene_symbol"].str.strip())


def load_overlay(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(require(path), dtype=str)
    frame["in_network_input"] = frame["in_network_input"].str.upper().eq("TRUE")
    for column in ["actual_kME", "projected_kME", "max_abs_projected_kME"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def build_candidate_audit(targets: list[str]) -> pd.DataFrame:
    if FULL_PROBE_STATS.exists():
        probe_stats = pd.read_csv(FULL_PROBE_STATS, compression="gzip")
        probe_stats = probe_stats.sort_values(["mad", "gene_symbol"], ascending=[False, True]).drop_duplicates("gene_symbol")
        probe_stats["full_transcriptome_mad_rank"] = np.arange(1, len(probe_stats) + 1)
        probe_stats = probe_stats.set_index("gene_symbol")
    else:
        probe_stats = pd.DataFrame(index=targets)

    selections: dict[str, pd.DataFrame] = {}
    for label, path in {
        "hvg8000": OUTPUT.parent / "tcga_coad_read_wgcna" / "tcga_coad_read_wgcna_selected_genes.csv",
        "hvg15000": OUTPUT.parent / "tcga_coad_read_wgcna_sensitivity_15000" / "tcga_coad_read_wgcna_selected_genes.csv",
    }.items():
        selections[label] = load_selection(path).set_index("gene_symbol")

    overlays = {
        "signed_bicor_8000": load_overlay(OUTPUT.parent / "tcga_coad_read_wgcna" / "tcga_coad_read_wgcna_41_gene_overlay.csv"),
        "signed_bicor_15000": load_overlay(OUTPUT.parent / "tcga_coad_read_wgcna_sensitivity_15000" / "tcga_coad_read_wgcna_41_gene_overlay.csv"),
        "unsigned_pearson_8000": load_overlay(OUTPUT.parent / "tcga_coad_read_wgcna_sensitivity_unsigned_pearson" / "tcga_coad_read_wgcna_41_gene_overlay.csv"),
    }
    overlay_maps = {label: frame.set_index("gene_symbol") for label, frame in overlays.items()}

    rows = []
    for gene in targets:
        row: dict[str, object] = {"gene_symbol": gene}
        if gene in probe_stats.index:
            row["full_transcriptome_mad"] = float(probe_stats.loc[gene, "mad"])
            row["full_transcriptome_mad_rank"] = int(probe_stats.loc[gene, "full_transcriptome_mad_rank"])
        else:
            row["full_transcriptome_mad"] = np.nan
            row["full_transcriptome_mad_rank"] = np.nan
        for label, selection in selections.items():
            row[f"{label}_rank"] = selection.loc[gene, "selection_rank"] if gene in selection.index else np.nan
            row[f"{label}_included"] = gene in selection.index
        for label, overlay in overlay_maps.items():
            present = bool(overlay.loc[gene, "in_network_input"])
            row[f"{label}_in_network"] = present
            row[f"{label}_module"] = overlay.loc[gene, "actual_module"] if present else ""
            row[f"{label}_kME"] = overlay.loc[gene, "actual_kME"] if present else np.nan
            row[f"{label}_projected_module"] = overlay.loc[gene, "projected_module"]
            row[f"{label}_projected_kME"] = overlay.loc[gene, "projected_kME"]
        rows.append(row)
    return pd.DataFrame(rows)


def sensitivity_summary() -> pd.DataFrame:
    rows = []
    specs = {
        "signed_bicor_8000": OUTPUT.parent / "tcga_coad_read_wgcna",
        "signed_bicor_15000": OUTPUT.parent / "tcga_coad_read_wgcna_sensitivity_15000",
        "unsigned_pearson_8000": OUTPUT.parent / "tcga_coad_read_wgcna_sensitivity_unsigned_pearson",
    }
    for label, directory in specs.items():
        manifest = json.loads(require(directory / "tcga_coad_read_wgcna_manifest.json").read_text(encoding="utf-8"))
        enrichment = pd.read_csv(require(directory / "tcga_coad_read_wgcna_target_module_enrichment.csv"))
        non_grey = enrichment[enrichment["module"].astype(str).ne("grey")].copy()
        if non_grey.empty:
            best = None
        else:
            best = non_grey.sort_values("p_value").iloc[0]
        rows.append({
            "analysis": label,
            "network_type": manifest["parameters"]["networkType"],
            "cor_type": manifest["parameters"]["corType"],
            "network_genes": manifest["network_genes"],
            "targets_in_network": manifest["target_genes_present_in_network"],
            "targets_projection_only": manifest["target_genes_projection_only"],
            "non_grey_modules": manifest["module_count_non_grey"],
            "soft_power": manifest["soft_threshold_power"],
            "best_non_grey_module": "" if best is None else best["module"],
            "best_non_grey_target_overlap": np.nan if best is None else best["target_overlap"],
            "best_non_grey_OR": np.nan if best is None else best["odds_ratio"],
            "best_non_grey_p": np.nan if best is None else best["p_value"],
            "best_non_grey_BH_FDR": np.nan if best is None else best["BH_FDR_across_WGCNA_modules"],
        })
    return pd.DataFrame(rows)


def write_report(candidate: pd.DataFrame, pairs: pd.DataFrame, sensitivity: pd.DataFrame) -> None:
    special = ["RXRA", "PPARA", "PPARD", "PPARG", "PTGS2", "NR1I2", "CYP2C9", "CYP3A4"]
    lines = [
        "# TCGA-COAD + READ WGCNA audit",
        "",
        "## Scope",
        "",
        "This is a sensitivity audit of the full-transcriptome WGCNA. The canonical network remains the 8,000-gene signed-bicor analysis. The 15,000-gene and unsigned-Pearson runs are sensitivity analyses; none changes the primary result.",
        "",
        "## Audit 1 — candidate eligibility and module placement",
        "",
        f"- Fresh target set: {len(candidate)} genes.",
        f"- Canonical 8,000-gene network input: {int(candidate['signed_bicor_8000_in_network'].sum())}/41.",
        f"- 15,000-gene sensitivity network input: {int(candidate['signed_bicor_15000_in_network'].sum())}/41.",
        f"- Unsigned-Pearson 8,000-gene network input: {int(candidate['unsigned_pearson_8000_in_network'].sum())}/41.",
        "",
        "The candidate-level CSV records full-transcriptome MAD rank, inclusion status, module, kME, and projection-only module for every target. Projection-only assignments are not counted as module membership or Fisher enrichment.",
        "",
        "### Priority genes",
        "",
        "| Gene | Full MAD rank | 8k module | 8k kME | 15k module | 15k kME | Unsigned module | Unsigned kME |",
        "|---|---:|---|---:|---|---:|---|---:|",
    ]
    for gene in special:
        row = candidate.loc[candidate["gene_symbol"] == gene].iloc[0]
        def fmt(value: object) -> str:
            if pd.isna(value) or value == "":
                return "—"
            if isinstance(value, (float, np.floating)):
                return f"{value:.3f}"
            return str(value)
        lines.append("| " + " | ".join([
            gene, fmt(row["full_transcriptome_mad_rank"]), fmt(row["signed_bicor_8000_module"]),
            fmt(row["signed_bicor_8000_kME"]), fmt(row["signed_bicor_15000_module"]),
            fmt(row["signed_bicor_15000_kME"]), fmt(row["unsigned_pearson_8000_module"]),
            fmt(row["unsigned_pearson_8000_kME"]),
        ]) + " |")
    lines += [
        "",
        "## Audit 2 — 41-gene pairwise co-expression",
        "",
        f"- Samples: 380 combined COAD/READ primary tumors.",
        f"- Pairwise correlations: {len(pairs)} unique gene pairs.",
        f"- Pairs with |r| >= 0.70: {int((pairs['abs_r'] >= 0.70).sum())}.",
        f"- Pairs with |r| >= 0.80: {int((pairs['abs_r'] >= 0.80).sum())}.",
        "",
        "The square correlation matrix and complete pair table are provided as audit outputs. This analysis is descriptive and does not establish chemical causality.",
        "",
        "## Audit 3 — WGCNA parameter sensitivity",
        "",
        "| Analysis | Network | Correlation | Genes | Targets in network | Non-grey modules | Best non-grey overlap | Best raw P | Best module BH-FDR |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in sensitivity.iterrows():
        lines.append("| " + " | ".join([
            str(row["analysis"]), str(row["network_type"]), str(row["cor_type"]), str(int(row["network_genes"])),
            str(int(row["targets_in_network"])), str(int(row["non_grey_modules"])),
            f"{row['best_non_grey_module']} ({int(row['best_non_grey_target_overlap'])}/{int(row['targets_in_network'])})",
            f"{row['best_non_grey_p']:.3g}", f"{row['best_non_grey_BH_FDR']:.3g}",
        ]) + " |")
    lines += [
        "",
        "## Interpretation",
        "",
        "The 15,000-gene network increased target inclusion, but retained only three non-grey modules and did not produce module-level BH-FDR<0.05 enrichment. The unsigned-Pearson network produced four non-grey modules, but likewise had no module-level BH-FDR<0.05 enrichment. Therefore the canonical negative module-enrichment result is not explained solely by the 8,000-gene filter or by signed-bicor choice.",
        "",
        "The audit does not prove absence of co-expression biology: the 41-gene pairwise matrix and projection results show descriptive relationships, and bulk modules can reflect tissue composition. It supports reporting WGCNA as non-confirmatory for a unified DINP–CRC module under these cohorts and parameterizations.",
    ]
    (OUTPUT / "TCGA_COAD_READ_WGCNA_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    targets = load_targets()
    expression = load_expression(targets)
    matrix, pairs = correlation_audit(expression)
    candidate = build_candidate_audit(targets)
    sensitivity = sensitivity_summary()
    matrix.to_csv(OUTPUT / "tcga_coad_read_41_gene_correlation_matrix.csv")
    pairs.to_csv(OUTPUT / "tcga_coad_read_41_gene_correlation_pairs.csv", index=False)
    candidate.to_csv(OUTPUT / "tcga_coad_read_wgcna_candidate_audit.csv", index=False)
    sensitivity.to_csv(OUTPUT / "tcga_coad_read_wgcna_sensitivity_comparison.csv", index=False)
    write_report(candidate, pairs, sensitivity)
    print(json.dumps({
        "targets": len(targets),
        "samples": len(expression),
        "pairs": len(pairs),
        "canonical_in_network": int(candidate["signed_bicor_8000_in_network"].sum()),
        "hvg15000_in_network": int(candidate["signed_bicor_15000_in_network"].sum()),
        "unsigned_in_network": int(candidate["unsigned_pearson_8000_in_network"].sum()),
        "output": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
