# URXP02 M4 sex-biased tissue mapping

Generated 2026-08-30T08:15:51.395029+00:00. This is a GTEx bulk-tissue expression-context audit, not a causal or cell-type analysis.

## Data and model

- **GTEx release:** GTEx V11; GENCODE v47 / GRCh38.
- **Tissue panel:** Thyroid; Aorta, Tibial, and Coronary arteries; left ventricle and atrial appendage; kidney cortex; adrenal gland; and liver as a toxicokinetic reference only.
- **Expression:** GTEx RNASeQCv2.4.3 gene TPM; analysed as `log2(TPM + 1)`.
- **Per-gene model:** female-minus-male coefficient adjusted for age-bracket midpoint, RIN, ischemic time, sequencing center, and RNA-extraction protocol. Fine expression batch (`SMGEBTCH`) is retained in the sample audit but not fit because its high cardinality is not stable in all prespecified tissues.
- **FDR:** fixed branch families: shared core 189 × 9 = 1701 tests; thyroid branch 219 × 9 = 1971; hypertension branch 440 × 9 = 3960. Failed/missing planned fits enter the corresponding BH family as p=1.

## Sample availability

- **Thyroid**: 459 male and 225 female complete-case samples (662 target genes found).
- **Artery - Aorta**: 308 male and 164 female complete-case samples (662 target genes found).
- **Artery - Tibial**: 471 male and 220 female complete-case samples (662 target genes found).
- **Artery - Coronary**: 168 male and 100 female complete-case samples (662 target genes found).
- **Heart - Left Ventricle**: 308 male and 143 female complete-case samples (662 target genes found).
- **Heart - Atrial Appendage**: 318 male and 143 female complete-case samples (662 target genes found).
- **Kidney - Cortex**: 80 male and 24 female complete-case samples (662 target genes found).
- **Adrenal Gland**: 184 male and 111 female complete-case samples (662 target genes found).
- **Liver**: 184 male and 78 female complete-case samples (662 target genes found).

## Shared-core summary

- Shared-core gene-by-tissue FDR-significant tests: **64**.
- Male-biased shared-core FDR hits in thyroid: **16**.
- Female-biased shared-core FDR hits in the non-thyroid disease-relevant panel: **9**.
- Rank-based shared-core gene-set shifts that survived its nine-tissue FDR family: **0/9**.
- Of the 828 frozen molecular-universe symbols, **662** mapped to GTEx V11 gene symbols in every panel tissue. The remaining planned symbols are retained as unavailable rows and enter the fixed FDR families as p=1.

These counts test the proposed directional pattern; they do not establish that the expression differences explain the NHANES associations.

No tissue showed a shared-core-wide FDR-significant shift. Accordingly, the gene-level hits are localized M5 handoff candidates rather than evidence for a global male-thyroid/female-vascular shared-core program.


## Evidence caveat carried forward

The thyroid-specific 30-gene set has no exact 2-NAP human or experimental support in M1b and is parent-naphthalene-only. It is retained for branch comparison, but is explicitly not interpreted as an exact 2-NAP mechanism.

## M5 handoff

The transparent M5 gate selected **6** shared-core candidates: JUN, TP53, AHR, CYP19A1, NFE2L2, NFKB1.

No figures, single-cell analyses, new NHANES models, pathway/PPI reruns, or causal claims were produced.
