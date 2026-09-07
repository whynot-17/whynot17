# DINP exposure-side source integration

Generated: `2026-09-07T13:27:48.611337+00:00`

## Scope

The parent chemical was fixed to diisononyl phthalate (DINP), CASRN 28553-12-0, CTD C012125, DTXSID4022521, PubChem CID 590836, and T3DB T3D3648. CTD, EPA CompTox/ToxCast/Tox21, and T3DB evidence are preserved as separate layers; support counts are descriptive and are not a biological truth score.

## Integrated counts

| Layer | Result |
|---|---:|
| CTD human interaction rows | 166 |
| CTD unique human genes | 86 |
| CompTox dashboard bioactivity endpoints | 445 |
| CompTox active endpoints | 11 |
| Active human ToxCast target genes | 5 |
| Active human Tox21 target genes | 3 |
| T3DB target rows | 1 |
| Expanded DINP exposure union | 93 |
| Genes supported by ≥2 source layers | 1 |

## CompTox active human targets

`BACE1, CYP1A2, CYP2C9, CYP3A4, F3, IL1A, NR1I2, SELE`

The CompTox extraction retained all 445 endpoint rows in the assay-level CSV. The gene layer uses only human-labelled active endpoints with an official gene symbol. Target-free active endpoints remain in the assay export but do not contribute a gene. Inactive or missing records are not treated as negative evidence.

## T3DB

The public DINP T3DB record is T3D3648 and reports one human target, NR1I2/PXR (UniProt O75469). The raw PubChem PUG-View response is retained because it cross-references the T3DB source record; the curated T3DB target row is kept separately.

## CRC overlay

The optional overlay against the existing GeneCards relevance-score ≥10 intersection contains 16 genes. It is an overlay only; it does not alter the prior disease-side threshold or multiplicity analysis.

## Outputs

- `dinp_toxcast_tox21_bioactivity.csv`: source-level CompTox assay/end-point evidence.
- `dinp_t3db_target_records.csv`: curated T3DB DINP target row.
- `dinp_exposure_gene_matrix_ctd_toxcast_tox21_t3db.csv`: source-preserving gene matrix.
- `dinp_crc_genecards_relevance10_expanded_exposure_overlay.csv`: optional exposure expansion over the existing CRC GeneCards threshold intersection.
- `dinp_exposure_source_integration_manifest.json`: counts, rules, provenance, and hashes.
