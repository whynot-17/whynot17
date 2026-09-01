# DHODH expression versus dependency in CRC

## Question

Do CRC models with higher DHODH expression show a more negative DHODH CRISPR/Chronos gene effect? This is the first, constitutive dependency check and does not claim an acquired OXA-R dependency.

## Design

- CRC restriction: DepMap `OncotreePrimaryDisease = Colorectal Adenocarcinoma` with known biological sex.
- Expression: DepMap 24Q4 expression matrix, DHODH column.
- Dependency: DepMap 26Q1 `CRISPRGeneEffect`/Chronos gene effect, DHODH column.
- Statistical unit: one DepMap model; more-negative gene effect means greater fitness loss after knockout.
- Prespecified association: Spearman/Pearson correlation and HC3-robust linear regression; no genome-wide multiple-testing correction because only one gene pair was tested.

## Matched panel

n=45 CRC models; sex counts={'Male': 27, 'Female': 18}.

## Result

Spearman rho=-0.1783, P=0.2414; Pearson r=-0.1445, P=0.3435.
HC3-robust OLS beta for expression=-0.06883, SE=0.05367, P=0.1997, R²=0.02089.
Mean gene effect=-0.5612; models with effect ≤−0.5: 25; ≤−1: 4.

## Interpretation

A negative expression–effect association would be compatible with higher DHODH expression marking greater constitutive DHODH dependency. A null result would mean that DHODH upregulation in OXA-R cannot be justified as a general CRC DHODH dependency from this analysis alone.

Even a positive association is not evidence that OXA resistance caused the dependency. The next step, only if warranted, is to project the OXA-R DHODH-centered transition state into the CRC panel and test whether that context strengthens the association.

## Provenance

{
  "gene": "DHODH",
  "model_file": "C:\\Users\\21634\\Documents\\Codex\\2026-08-31\\1-janssen-kp-et-al-extrinsic\\work\\depmap_Model.csv",
  "expression_file": "C:\\Users\\21634\\Documents\\Codex\\2026-08-31\\1-janssen-kp-et-al-extrinsic\\work\\depmap_expression_24Q4.csv",
  "dependency_file": "C:\\Users\\21634\\Documents\\Codex\\2026-08-31\\1-janssen-kp-et-al-extrinsic\\work\\data\\CRISPRGeneEffect_26Q1.csv",
  "expression_column": "DHODH (1723)",
  "dependency_column": "DHODH (1723)",
  "disease_filter": "OncotreePrimaryDisease == Colorectal Adenocarcinoma",
  "sex_filter": [
    "Male",
    "Female"
  ],
  "n_model_metadata": 86,
  "n_expression": 1480,
  "n_dependency": 1208,
  "n_matched": 45,
  "sex_counts": {
    "Male": 27,
    "Female": 18
  },
  "dependency_release": "DepMap 26Q1 CRISPRGeneEffect / Chronos gene effect",
  "expression_release": "DepMap 24Q4 expression"
}
