from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sex_oxa_resistance_transition_analysis import (
    GEO,
    PLATFORM,
    annotated_gene_expression,
    parse_samples,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


MODULES = {
    # Core pathway from carbamoyl-phosphate synthesis through dihydroorotate
    # oxidation and UMP formation; kept separate from nucleotide demand genes.
    "de_novo_pyrimidine_core": ["CAD", "DHODH", "UMPS"],
    # Broader pyrimidine/nucleotide supply module used as a sensitivity analysis.
    "pyrimidine_nucleotide_supply": [
        "CAD", "DHODH", "UMPS", "CTPS1", "CTPS2", "PPAT", "PAICS",
        "RRM1", "RRM2", "DCTD", "TYMS", "DHFR", "TK1", "DUT", "CMPK1", "UCK2",
    ],
}
GENES = list(dict.fromkeys(gene for genes in MODULES.values() for gene in genes))
FOCAL_GENES = ["DHODH", "RRM2", "DCTD", "CAD", "UMPS", "TYMS"]


def zscore_columns(expr: pd.DataFrame) -> pd.DataFrame:
    return (expr - expr.mean(axis=0)) / expr.std(axis=0, ddof=1).replace(0, np.nan)


def score_modules(expr: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    z = zscore_columns(expr)
    scores = pd.DataFrame(index=expr.index)
    present: dict[str, list[str]] = {}
    for module, genes in MODULES.items():
        present[module] = [gene for gene in genes if gene in z.columns]
        scores[module] = z[present[module]].mean(axis=1, skipna=True) if present[module] else np.nan
    return scores, present


def load_dataset(dataset: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    expr, expr_manifest = annotated_gene_expression(dataset, {"pyrimidine_targets": GENES})
    meta = parse_samples(dataset).set_index("gsm")
    keep = meta.index.intersection(expr.index)
    expr = expr.reindex(keep)
    meta = meta.reindex(keep)
    scores, present = score_modules(expr)
    sample = meta.join(scores)
    manifest = {
        **expr_manifest,
        "dataset": dataset,
        "target_genes": GENES,
        "genes_present": sorted(expr.columns.tolist()),
        "module_genes_present": present,
        "n_selected_samples": int(len(sample)),
    }
    return expr, sample, manifest


def make_pair_tables(expr: pd.DataFrame, sample: pd.DataFrame, dataset: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Each parental→OxR transition is represented once after averaging technical
    # replicate arrays within each state.
    pair_rows = []
    gene_rows = []
    for (cell_line, sex), group in sample.groupby(["cell_line", "sex"], sort=False):
        parental = group[group.state.eq("parental")]
        oxr = group[group.state.eq("OxR")]
        if parental.empty or oxr.empty:
            continue
        p_expr = expr.reindex(parental.index).mean(axis=0)
        r_expr = expr.reindex(oxr.index).mean(axis=0)
        row = {
            "dataset": dataset,
            "cell_line": cell_line,
            "sex": sex,
            "n_parental_arrays": int(len(parental)),
            "n_oxr_arrays": int(len(oxr)),
        }
        p_scores = sample.reindex(parental.index)[list(MODULES)].mean(axis=0)
        r_scores = sample.reindex(oxr.index)[list(MODULES)].mean(axis=0)
        for module in MODULES:
            row[f"parental_{module}"] = p_scores[module]
            row[f"oxr_{module}"] = r_scores[module]
            row[f"delta_{module}"] = r_scores[module] - p_scores[module]
        for gene in FOCAL_GENES:
            row[f"parental_{gene}"] = p_expr.get(gene, np.nan)
            row[f"oxr_{gene}"] = r_expr.get(gene, np.nan)
            row[f"delta_{gene}"] = r_expr.get(gene, np.nan) - p_expr.get(gene, np.nan)
        pair_rows.append(row)
        for gene in sorted(set(GENES) | set(FOCAL_GENES)):
            gene_rows.append(
                {
                    "dataset": dataset,
                    "cell_line": cell_line,
                    "sex": sex,
                    "gene": gene,
                    "parental_mean": p_expr.get(gene, np.nan),
                    "oxr_mean": r_expr.get(gene, np.nan),
                    "delta_OxR_minus_parental": r_expr.get(gene, np.nan) - p_expr.get(gene, np.nan),
                }
            )
    return pd.DataFrame(pair_rows), pd.DataFrame(gene_rows)


def summarize(pair_deltas: pd.DataFrame, gene_deltas: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    module_rows = []
    for module in MODULES:
        column = f"delta_{module}"
        for sex in ["Female", "Male", "All"]:
            d = pair_deltas if sex == "All" else pair_deltas[pair_deltas.sex.eq(sex)]
            values = d[column].dropna()
            module_rows.append(
                {
                    "module": module,
                    "sex": sex,
                    "n_transitions": int(len(values)),
                    "mean_delta": values.mean(),
                    "median_delta": values.median(),
                    "n_down": int((values < 0).sum()),
                    "n_up": int((values > 0).sum()),
                    "fraction_down": (values < 0).mean() if len(values) else np.nan,
                }
            )
    gene_rows = []
    for gene in sorted(set(GENES) | set(FOCAL_GENES)):
        for sex in ["Female", "Male", "All"]:
            d = gene_deltas if sex == "All" else gene_deltas[gene_deltas.sex.eq(sex)]
            values = d.loc[d.gene.eq(gene), "delta_OxR_minus_parental"].dropna()
            gene_rows.append(
                {
                    "gene": gene,
                    "sex": sex,
                    "n_transitions": int(len(values)),
                    "mean_delta": values.mean(),
                    "median_delta": values.median(),
                    "n_down": int((values < 0).sum()),
                    "n_up": int((values > 0).sum()),
                    "fraction_down": (values < 0).mean() if len(values) else np.nan,
                }
            )
    return pd.DataFrame(module_rows), pd.DataFrame(gene_rows)


def make_figure(pair_deltas: pd.DataFrame, gene_deltas: pd.DataFrame) -> None:
    colors = {"Female": "#c45a72", "Male": "#3b6fb6"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=220)
    modules = list(MODULES)
    x = np.arange(len(modules))
    for sex in ["Female", "Male"]:
        d = pair_deltas[pair_deltas.sex.eq(sex)]
        means = [d[f"delta_{module}"].mean() for module in modules]
        axes[0].plot(x, means, marker="o", lw=2, color=colors[sex], label=f"{sex} (n={len(d)})")
        for _, row in d.iterrows():
            axes[0].plot(x, [row[f"delta_{module}"] for module in modules], color=colors[sex], alpha=0.25, lw=0.8)
    axes[0].axhline(0, color="#555", lw=0.8)
    axes[0].set_xticks(x, ["de novo core", "nucleotide supply"], rotation=15, ha="right")
    axes[0].set_ylabel("Δ score: OxR − parental")
    axes[0].set_title("Pyrimidine pathway rewiring")
    axes[0].legend(frameon=False)

    genes = FOCAL_GENES
    for i, gene in enumerate(genes):
        for sex, offset in [("Female", -0.16), ("Male", 0.16)]:
            values = gene_deltas[(gene_deltas.gene.eq(gene)) & (gene_deltas.sex.eq(sex))]["delta_OxR_minus_parental"].dropna()
            if len(values):
                axes[1].scatter(np.repeat(i + offset, len(values)), values, color=colors[sex], alpha=0.8, s=42, label=sex if i == 0 else None)
                axes[1].plot([i + offset - 0.08, i + offset + 0.08], [values.mean(), values.mean()], color=colors[sex], lw=2)
    axes[1].axhline(0, color="#555", lw=0.8)
    axes[1].set_xticks(np.arange(len(genes)), genes, rotation=30, ha="right")
    axes[1].set_ylabel("Δ log2 expression: OxR − parental")
    axes[1].set_title("Focal gene transitions")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=0.2)
    fig.suptitle("CRC: acquired oxaliplatin resistance and pyrimidine rewiring", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "oxa_r_pyrimidine_rewiring.png", bbox_inches="tight")
    plt.close(fig)


def write_report(pair_deltas: pd.DataFrame, gene_deltas: pd.DataFrame, module_summary: pd.DataFrame, gene_summary: pd.DataFrame, manifest: dict[str, object]) -> None:
    lines = [
        "# OXA-R pyrimidine rewiring",
        "",
        "## Question",
        "",
        "Does acquired oxaliplatin resistance in CRC show a reproducible de novo pyrimidine/nucleotide-supply transition, providing a rational entry point for DHODH dependency analysis?",
        "",
        "## Design",
        "",
        "- Datasets: GSE42387 (HCT116, HT29, LoVo parental/OxPt-resistant) and GSE76092 (HT29/HTOXAR3 untreated basal comparison).",
        "- Expression: platform-annotated microarray expression, summarized to gene level by median across probes.",
        "- Unit of analysis: one parental→OxR transition per dataset×cell-line pair; replicate arrays are averaged within state before calculating Δ=OxR−parental.",
        "- Modules: a prespecified three-gene de novo core (CAD/DHODH/UMPS) and a broader pyrimidine/nucleotide-supply panel.",
        "",
        "## Available transitions",
        "",
        f"The analysis retained {len(pair_deltas)} transitions across {len(manifest['datasets'])} datasets. This is a small transition panel and is intended as a screen for a stable direction, not as a powered meta-analysis.",
        "",
        pair_deltas[["dataset", "cell_line", "sex", "delta_de_novo_pyrimidine_core", "delta_pyrimidine_nucleotide_supply", "delta_DHODH", "delta_RRM2", "delta_DCTD", "delta_CAD", "delta_UMPS", "delta_TYMS"]].to_string(index=False),
        "",
        "## Module summary",
        "",
        module_summary.to_string(index=False),
        "",
        "## Focal-gene summary",
        "",
        gene_summary[gene_summary.gene.isin(FOCAL_GENES)].to_string(index=False),
        "",
        "## Interpretation",
        "",
        "The key criterion for moving to dependency analysis is cross-transition directionality of the de novo core and/or a coherent DHODH-centered response. A mixed module score with isolated gene changes should be treated as rewiring rather than proof of pathway activation.",
        "",
        "No CRISPR dependency or pharmacologic vulnerability analysis is included in this first pass. Those will be run only if the rewiring pattern is sufficiently coherent.",
        "",
        "## Provenance",
        "",
        json.dumps(manifest, ensure_ascii=False, indent=2),
    ]
    (OUT / "oxa_r_pyrimidine_rewiring_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    all_pairs = []
    all_genes = []
    manifest: dict[str, object] = {"datasets": {}, "modules": MODULES, "focal_genes": FOCAL_GENES}
    for dataset in ["GSE42387", "GSE76092"]:
        expr, sample, info = load_dataset(dataset)
        pairs, genes = make_pair_tables(expr, sample, dataset)
        all_pairs.append(pairs)
        all_genes.append(genes)
        manifest["datasets"][dataset] = info
    pair_deltas = pd.concat(all_pairs, ignore_index=True)
    gene_deltas = pd.concat(all_genes, ignore_index=True)
    module_summary, gene_summary = summarize(pair_deltas, gene_deltas)
    pair_deltas.to_csv(OUT / "oxa_r_pyrimidine_pair_deltas.csv", index=False)
    gene_deltas.to_csv(OUT / "oxa_r_pyrimidine_gene_deltas.csv", index=False)
    module_summary.to_csv(OUT / "oxa_r_pyrimidine_module_summary.csv", index=False)
    gene_summary.to_csv(OUT / "oxa_r_pyrimidine_gene_summary.csv", index=False)
    (OUT / "oxa_r_pyrimidine_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    make_figure(pair_deltas, gene_deltas)
    write_report(pair_deltas, gene_deltas, module_summary, gene_summary, manifest)
    print("Pair transitions")
    print(pair_deltas[["dataset", "cell_line", "sex", "delta_de_novo_pyrimidine_core", "delta_pyrimidine_nucleotide_supply", "delta_DHODH", "delta_RRM2", "delta_DCTD", "delta_CAD", "delta_UMPS", "delta_TYMS"]].to_string(index=False))
    print("Module summary")
    print(module_summary.to_string(index=False))
    print("Focal gene summary")
    print(gene_summary[gene_summary.gene.isin(FOCAL_GENES)].to_string(index=False))


if __name__ == "__main__":
    main()
