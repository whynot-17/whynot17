# SLC7A11 virtual knockout redo

## Question

Does SLC7A11 CRISPR/Chronos gene effect differ between manually curated right- and left-sided COAD models? More-negative gene effect means a larger fitness loss after SLC7A11 knockout.

## Design

- Primary unit: one DepMap COAD cell model.
- Primary sidedness analysis: explicit high-confidence anatomical primary-site labels only.
- Sensitivity analysis: high-confidence plus three medium-confidence published assignments.
- Rectal and problematic/misclassified ovarian-site records are not included in the Right-versus-Left comparison.
- Secondary models test `SLC7A11_gene_effect ~ side + state + side:state` for both PUFA-pressure and AA-routing states; all state-by-side terms are exploratory and do not establish causality.

## Primary result

- Right: n=12; mean gene effect=0.10049788928824104
- Left: n=5; mean gene effect=0.14628157355646282
- Right-minus-Left mean difference=-0.04578368426822178
- Exact one-sided permutation p (Right more negative)=0.3015026660203587 (exact_6188_partitions)
- One-sided Mann–Whitney p=0.19133807369101485
- Cliff's delta (Right versus Left)=-0.3

## Sensitivity result

- Right: n=13; Left: n=7
- Right-minus-Left mean difference=-0.0348975543385068
- Exact one-sided permutation p=0.3138891397169799 (exact_77520_partitions)

## Pressure interaction

- Matched labeled models=13; Right-specific slope=-0.005098159331579949; linear-contrast P=0.9852369936915752
- Right-minus-Left slope interaction beta=0.14628977454477107; interaction P=0.6511624944161923
- The PUFA-pressure proxy is the mean within-COAD z-score of the supplied DepMap expression values for ACSL4, LPCAT3, ALOX5, ALOX12 and ALOX15, retained only when at least 4/5 genes are available. It is not a measurement of AA concentration or flux.

## AA-routing interaction

- Matched labeled models=13; Right-specific slope=-0.11875655775478966; linear-contrast P=0.31630379099101774
- Right-minus-Left slope interaction beta=0.086605291998123; interaction P=0.6233194105059274
- The AA-routing proxy is the mean within-COAD z-score of PLA2G4A, PTGS2, ALOX5 and ALOX5AP, retained when at least 3/4 genes are available. It is a transcriptional proxy, not a direct eicosanoid or AA-flux measurement.

## Interpretation guardrails

A negative Right-minus-Left effect estimate is compatible with stronger SLC7A11 knockout sensitivity in right-sided models. A null result does not prove equal biology because the curated cell-line sample is small and sidedness labels are incomplete. A positive result or a pressure interaction remains an association, not proof that AA caused the dependency.

## Provenance

{
  "model_file": "C:\\Users\\21634\\Documents\\Codex\\2026-08-31\\1-janssen-kp-et-al-extrinsic\\work\\depmap_Model.csv",
  "dependency_file": "C:\\Users\\21634\\Documents\\Codex\\2026-08-31\\1-janssen-kp-et-al-extrinsic\\work\\data\\CRISPRGeneEffect_26Q1.csv",
  "expression_file": "C:\\Users\\21634\\Documents\\Codex\\2026-08-31\\1-janssen-kp-et-al-extrinsic\\work\\depmap_expression_24Q4.csv",
  "dependency_column": "SLC7A11 (23657)",
  "pressure_expression_columns": {
    "ACSL4": "ACSL4 (2182)",
    "LPCAT3": "LPCAT3 (10162)",
    "ALOX5": "ALOX5 (240)",
    "ALOX12": "ALOX12 (239)",
    "ALOX15": "ALOX15 (246)"
  },
  "aa_routing_expression_columns": {
    "PLA2G4A": "PLA2G4A (5321)",
    "PTGS2": "PTGS2 (5743)",
    "ALOX5": "ALOX5 (240)",
    "ALOX5AP": "ALOX5AP (241)"
  },
  "pufa_min_available_genes": 4,
  "aa_routing_min_available_genes": 3,
  "dependency_release": "DepMap Public 26Q1 CRISPRGeneEffect / Chronos gene effect",
  "expression_release": "DepMap 24Q4 expression",
  "manual_curation_entries": 22,
  "seed": 20260901,
  "no_external_download": true
}
