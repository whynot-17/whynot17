# M5b — URXP02 single-cell rescue audit

## Scope

This audit only verifies dataset availability, donor sex metadata, healthy/control composition, processed/raw data access, annotations, and donor-level eligibility. No candidate-gene expression, differential expression, module analysis, clustering, reannotation, figures, or new epidemiologic analysis was run.

## Decision summary

- **Primary recommendation: GSE202109 kidney.** It contains 19 healthy living kidney donors (9 male, 10 female) and 27,677 cells, satisfying the prespecified minimum of at least three donors per sex.
- **Secondary kidney validation: GSE183276 and GSE183277.** These are mixed kidney atlases, but their deposited `Normal Reference` subsets are sex-balanced (scCv3: 20 donors, 7M/13F, 21,650 cells; snCv3: 18 donors, 8M/10F, 88,460 cells). Disease samples must remain separate.
- **Vascular backup: GSE207784 controls.** The series has 13 donors and 71,689 nuclei overall; the seven non-aneurysm controls are 3M/4F. This meets the minimum numerically but is low-powered.
- **Thyroid rescue failed for sex comparison in GSE182416.** All seven samples are female, and the sample characteristics contain PTC or follicular adenoma pathology. It can provide a cell-type reference, not donor-sex DE.
- **GSE189795 is not eligible for sex comparison.** All four normal controls are male (and all nine series samples are male).
- **SCP1265/GSE165824 is a separate thoracic-aorta study, not a duplicate of GSE207784.** It has three paired normal donors and 54,092 nuclei. The deposited h5ad `obs` metadata identifies two female donors (Ao4, Ao8) and one male donor (Ao12), but this is below the preferred donor threshold and is descriptive only.

## Eligibility rule

The preferred criterion is at least 3 male and 3 female donors in the eligible healthy/control subset. A dataset can be retained below that threshold for descriptive or exploratory purposes, but it is flagged low-powered and is not treated as a confirmatory donor-sex resource.

## Audit provenance

The audit used official GEO family SOFT records for GSE182416, GSE202109, GSE207784, GSE189795, and GSE165824, deposited kidney atlas metadata for GSE183276/GSE183277, and metadata-only inspection of the deposited GSE165824 h5ad. The kidney metadata parser explicitly accounts for the deposited extra cell-barcode field before the published header. The incomplete local GSE207784 h5ad download was not used to infer counts; the 71,689-nucleus total is taken from the official GEO series record. For GSE165824, only h5ad dimensions and donor/cell-type/disease metadata were read; expression values were not accessed.

Generated: 2026-08-30T12:48:54.371977+00:00
