# DINP target mining search log

Run date: 2026-09-08

## Scope and non-goals

Scope was limited to DINP (diisononyl phthalate) chemical identity resolution, public database target/interaction records, public DINP exposure transcriptomics, and primary literature. No CRC overlap, MR, TCGA, single-cell, PPI, GO/KEGG, docking, or disease-specific integration was performed.

## Chemical entity control

Primary entity: PubChem CID 590836, ChEMBL1905899, ChEBI:35459, formula C26H42O4, InChIKey HBGGXOJOCNVPFY-UHFFFAOYSA-N. DINP-associated identifiers include CAS 28553-12-0 and technical-mixture CAS 68515-48-0; CTD chemical ID is C012125 and CTD maps it to DTXSID4022521. Search aliases included diisononyl phthalate, di-isononyl phthalate, diisononylphthalate, DINP/DiNP, ENJ 2065, Jayflex DINP, Vestinol NN, Palatinol DINP, and phthalic acid diisononyl ester.

Exclusion control: raw acronym DINP was not used alone because it is ambiguous; dinonyl phthalate (PubChem CID 6787), DEHP, DIDP, DBP, MEHP, DINCH and DEHT were not treated as DINP parent entities. Metabolite-only CTD rows and co-treated/mixed-phthalate CTD rows were excluded from the DINP-only count.

## Evidence grading used

A = primary wet experimental receptor activity, expression/protein response, phosphorylation/activation, or functional perturbation explicitly involving DINP; A does not mean purified-protein binding. B = CTD curated interaction or public transcriptomic response, plus ChEMBL-indexed PubChem functional records when the activity comment is inconclusive. C = explicitly in-silico binding/prediction only. C evidence is marked prediction-only in the evidence partition and is never silently upgraded by database wording.

## Retrieval and filtering audit

| Source | Query/record | Result and filter |
|---|---|---|
| CTD | C012125; CTD_chem_gene_ixns.tsv | 1268 raw DINP rows → 296 standalone DINP rows → 251 raw gene symbols; 251 raw symbols retained before human canonical validation; 232 human canonical genes in output. Co-treated/mixed rows and 4 parent-metabolite rows excluded. |
| T3DB | T3D3648 | 7 target associations. NR1I2/NR1I3 linked to PMID 21227907 and retained as A via primary experimental study; PPARA, PPARD, PPARG, RXRA, RXRB linked to PMID 23843199 and retained as C prediction-only evidence. |
| ChEMBL/PubChem BioAssay | CHEMBL1905899 exact molecule | 2 functional activities: THRB / AID 588545 / 398.1 nM and NFE2L2 / AID 720636 / 3.0638 µM. Both ChEMBL records have activity_comment=inconclusive; no direct-binding A label. |
| BindingDB | exact DINP SMILES; cutoff=1.0 | 0 exact hits; similar-compound targets were not imported. |
| STITCH | PubChem/CID 590836 resolution | Current STITCH resolve returned CID identifiers, but current network/interactions endpoints returned 404; 0 STITCH target rows imported. |
| EPA CompTox/ToxCast | DTXSID4022521 / DTXSID60860420 | Dashboard/API and current v4.3 source were identified. API access requires a free API key and no key was available in this run; no ToxCast assay rows were imported, so no EPA A claims are made. |
| GEO GSE40342 | DINP rat primary hepatocytes | 30 samples; 4 DINP-vs-vehicle contrasts; 21,419 annotated probes used. Conservative output filter q<0.05 and |log2FC|≥1.0 yielded 64 top raw candidate symbols before human-symbol validation. |
| GEO GSE158473 | adult mouse ovary, DiNP and DEHP arms | 30 samples; DiNP arms were analyzed against time-matched vehicle only, never against DEHP. With n=3 per group and q<0.05/|log2FC|≥0.5, no candidate passed; dataset was logged but not added to master. |
| PubMed | exact title query: diisononyl phthalate[Title] | 67 title hits; raw DINP acronym searches were not used as a primary query. Primary records were screened for DINP-specific arms and retractions. |

## Public source URLs

- PubChem: https://pubchem.ncbi.nlm.nih.gov/compound/590836
- CTD chemical–gene interaction download: https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz
- T3DB DINP card: https://www.t3db.ca/toxins/T3D3648
- T3DB MOA/target download: https://www.t3db.org/system/downloads/current/moas.csv.zip
- ChEMBL molecule: https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL1905899.json
- ChEMBL activities: https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=CHEMBL1905899&limit=1000
- BindingDB exact-SMILES endpoint: https://bindingdb.org/rest/getTargetByCompound
- STITCH: https://stitch-db.org/
- EPA CompTox API information: https://www.epa.gov/comptox-tools/computational-toxicology-and-exposure-apis-about
- EPA ToxCast data exploration: https://www.epa.gov/comptox-tools/exploring-toxcast-data
- GEO GSE40342: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40342&targ=self&form=text&view=brief
- GEO GSE158473: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158473&targ=self&form=text&view=brief
- PubMed exact-title query: https://pubmed.ncbi.nlm.nih.gov/?term=%28diisononyl+phthalate%5BTitle%5D%29

## Output summary

- Unique human canonical targets: **286**. Counts are based on the highest grade seen per unique target: A=27, B=258, C=1.
- Evidence-item counts before target-level aggregation: A=34, B=291, C=10.
- Prediction-only targets: **1**; targets with any prediction evidence: 5.
- Unique targets supported by ≥2 source families: **32**. Source-family definition: CTD, T3DB, ChEMBL/PubChem BioAssay, GEO, and distinct primary-literature PMIDs; CTD/T3DB may trace back to overlapping primary papers, so this is a provenance count, not proof of biological independence.

Targets with ≥2 source families: NR1I2, NR1I3, PPARA, PPARG, CEBPA, CEBPB, CYP2B6, CYP3A4, DDIT4, FABP4, FABP5, FASN, GATA2, NFE2L2, NR3C1, RELA, SIRT1, SIRT2, SIRT3, SIRT5, STAT3, STAT6, TSLP, ACAA2, ACOX1, ACSL1, AQP7, PPARD, RXRA, STEAP4, THRB, RXRB.

| Source family | Evidence items | Unique genes |
|---|---:|---:|
| CTD | 238 | 232 |
| ChEMBL/PubChem BioAssay | 2 | 2 |
| GEO | 56 | 54 |
| Literature | 32 | 27 |
| T3DB | 7 | 7 |

## Important audit flags

1. The master table is a union of distinct evidence classes, but every row carries grade, source family, species, assay/context and caveat fields. The evidence summary is the compact per-gene audit view; prediction-only targets are explicitly labeled.
2. Non-human observations are supporting evidence only. Their final `gene_symbol` is retained only when the symbol could be directly validated against current human NCBI Gene_info; no unsupported species-specific symbol was promoted.
3. GSE40342 dose labels are recorded as GEO reports them (25 mM and 100 mM); units were not silently changed.
4. The GSE re-analysis is intentionally conservative and does not import all 1,000+ nominal DE genes. It uses the top tranche per contrast and |log2FC|≥1.0 to keep the target table auditable.
5. A retracted DINP brain-energy-metabolism parent study (PMID 39224312; retraction notice PMID 40201523) was excluded from evidence.
