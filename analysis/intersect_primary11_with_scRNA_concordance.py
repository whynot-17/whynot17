"""Intersect primary BH-FDR genes with previously reported cross-cohort single-cell patterns.

The concordance set is the eight unique genes explicitly summarized in the current
single-cell report as having a matched major cell compartment and same tumor-normal
expression direction in GSE132465 and GSE200997. Stromal cells map to inferred
fibroblasts; epithelial and immune labels are matched by corresponding annotation.
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/2026-10-11/single_cell_16genes_GSE132465_GSE200997"
INPUT = RESULTS / "sample_level_paired_tests.tsv"
OUTPUT = RESULTS / "primary11_scRNA_concordance_intersection.tsv"

CONCORDANT_PATTERNS = {
    "TIMP1": [("Epithelial cells", "Epithelial", "up"), ("Stromal cells", "Fibroblast", "up")],
    "SCD": [("Epithelial cells", "Epithelial", "up")],
    "UBE2C": [("Epithelial cells", "Epithelial", "up")],
    "PCK1": [("Epithelial cells", "Epithelial", "down")],
    "INHBA": [("Stromal cells", "Fibroblast", "up")],
    "TPM1": [("Stromal cells", "Fibroblast", "up")],
    "FHL1": [("Stromal cells", "Fibroblast", "down")],
    "PID1": [("Stromal cells", "Fibroblast", "down")],
}


def main() -> None:
    tests = pd.read_csv(INPUT, sep="\t")
    primary = tests[(tests.cohort == "GSE132465") & (tests.q_within_cohort_lt_0_05)]
    primary_genes = set(primary.gene.unique())
    rows = []
    for gene, patterns in CONCORDANT_PATTERNS.items():
        if gene not in primary_genes:
            continue
        for cell_132, cell_200, direction in patterns:
            a = tests[(tests.cohort == "GSE132465") & (tests.gene == gene) & (tests.cell_type == cell_132)]
            b = tests[(tests.cohort == "GSE200997") & (tests.gene == gene) & (tests.cell_type == cell_200)]
            if len(a) != 1 or len(b) != 1:
                raise ValueError(f"Missing or duplicate matched comparison: {gene}, {cell_132}/{cell_200}")
            x, y = a.iloc[0], b.iloc[0]
            rows.append({
                "gene": gene,
                "GSE132465_cell_type": cell_132,
                "GSE200997_cell_type": cell_200,
                "concordant_tumor_normal_direction": direction,
                "GSE132465_paired_n": int(x.n_complete_pairs),
                "GSE132465_tumor_minus_normal_mean_log1p_CP10K": x.mean_paired_difference_tumor_minus_normal,
                "GSE132465_primary_BH_q": x.q_value_BH_within_cohort,
                "GSE200997_paired_n": int(y.n_complete_pairs),
                "GSE200997_tumor_minus_normal_mean_log1p_CP10K": y.mean_paired_difference_tumor_minus_normal,
                "GSE200997_nominal_p": y.p_value,
                "GSE200997_within_cohort_BH_q": y.q_value_BH_within_cohort,
            })
    out = pd.DataFrame(rows)
    out.to_csv(OUTPUT, sep="\t", index=False, float_format="%.10g")
    print("Primary BH-FDR unique genes:", len(primary_genes), sorted(primary_genes))
    print("Cross-cohort concordance unique genes:", len(CONCORDANT_PATTERNS), sorted(CONCORDANT_PATTERNS))
    print("Intersection unique genes:", out.gene.nunique(), sorted(out.gene.unique()))
    print("Cell-type rows (TIMP1 has two):", len(out))
    print("Wrote:", OUTPUT)


if __name__ == "__main__":
    main()
