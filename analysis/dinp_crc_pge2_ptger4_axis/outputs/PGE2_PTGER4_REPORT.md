# Focused PGE2 -> PTGER4 axis analysis

## Question

Which cell group has the strongest **PGE2 synthesis program**, and which macrophage subtype has the strongest **PTGER4 receiving signal**? This is a descriptive molecular axis analysis and does not prove extracellular PGE2 production, binding, or physical cell-cell communication.

## Data and scoring

- Dataset: local GSE144735 CRC matrix (27,414 annotated cells; 6 donors: KUL01, KUL19, KUL21, KUL28, KUL30, KUL31).
- Synthesis genes: PLA2G4A, PTGS1, PTGS2, PTGES, PTGES2, PTGES3.
- Raw counts were library-size normalized to 10,000 and log1p transformed; each synthesis gene was z-scored across all cells and averaged into `PGE2_synthesis_score`.
- `PGE2_synthesis_high` is the global cell-level upper quartile (threshold 0.127).
- Macrophage PTGER4-high is the macrophage-compartment upper quartile of normalized log1p PTGER4; here the 75th percentile is zero, so the prespecified fallback is **detected-only** (threshold 0.000 (detected-only fallback because the macrophage 75th percentile is zero)).
- Donor-level tables are descriptive aggregates; no cell-level significance is used.

## Main readouts

- Highest broad-group PGE2 synthesis mean: **mast_cells**. Among the prespecified microenvironment comparison groups, the highest mean is **SPP1_like_TAM**.
- Highest macrophage PTGER4 mean: **C1QC_like_TAM**.
- Highest macrophage PTGER4 detection fraction: **C1QC_like_TAM**.
- Highest macrophage PTGER4-high fraction: **C1QC_like_TAM**.
- At the donor level, C1QC-like TAM is the top PTGER4 subtype by mean in **2/6** tumor donors and by detection fraction in **3/6**; it is therefore not uniformly highest across all six donors.
- The complete group and subtype rankings are in `pge2_synthesis_by_group.csv` and `ptger4_macrophage_subtype_statistics.csv`.

## Donor-level source -> receiver check

Tumor donor rows compare each candidate source with each macrophage subtype. A donor is directionally consistent for a candidate source -> C1QC-like TAM pair when that source is the top PGE2 synthesis group among mast cells plus the prespecified comparison groups **and** C1QC-like TAM is the top PTGER4 macrophage subtype for that donor. The product in `pge2_ptger4_donor_axis.csv` is a descriptive ranking score, not a communication probability.

For the most directionally consistent non-C1QC source pair, **mast_cells -> C1QC_like_TAM**, 1/6 tumor donors meet both top-ranking conditions.

Interpretation should use **supports/suggests a potential PGE2-producing source -> PTGER4-high macrophage axis**. It must not be written as demonstrated PGE2 communication without spatial, biochemical, or functional validation.

## Files

- `pge2_synthesis_by_group.csv`: group-level synthesis ranking.
- `ptger4_macrophage_subtype_statistics.csv`: PTGER4 mean, detection fraction, and high-cell fraction.
- `pge2_ptger4_axis_donor_summary.csv`: six-donor consistency summary for each source -> C1QC pair.
- `pge2_ptger4_top_calls_by_donor.csv`: per-donor top source and receiver calls; sparse groups remain visible through their cell counts.
- `figures/`: publication-ready PNG, PDF, and SVG summaries.
