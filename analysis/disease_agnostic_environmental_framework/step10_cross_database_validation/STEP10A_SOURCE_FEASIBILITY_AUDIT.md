# Step 10A — Source-replacement feasibility audit

## Locked status

`COMPLETE — source registry and frozen source set locked; Step 10B not run`

Step 10A audited source feasibility, evidence semantics, public/controlled
access, identifier systems, release handling, and independence from the current
CTD/GeneCards/NHANES architecture. No candidate-level result, MCOP result, T2D
association, CRC association, pathway result, or ranking was loaded or used to
choose a source.

The complete registry is in `STEP10_SOURCE_REPLACEMENT_REGISTRY.csv`, and the
predeclared source set is in `STEP10_FROZEN_SOURCE_SET.json`.

## Registry summary

| Layer | Registry entries | Primary replacements | Conditional / other |
|---|---:|---:|---:|
| Environmental | 6 | 3 | 3 |
| Disease | 5 | 2 | 3 |
| Epidemiology | 5 | 0 | 5 |
| Disease-context | 1 | 0 | supplementary-only |
| **Total** | **17** | **5** | **12** |

### Frozen primary replacement set

- Environmental: **E1 ChEMBL**, **E2 BindingDB**, **E3 PubChem BioAssay**.
- Disease: **D1 Open Targets Platform**, **D2 NHGRI-EBI GWAS Catalog**.
- Epidemiology: **none yet**.

The empty epidemiology primary set is intentional. A dataset cannot enter the
primary replacement set merely because it mentions environmental exposure and
diabetes. It must provide, in an approved individual-level extract, a
harmonizable exposure, a prespecified T2D outcome, core covariates, design
variables/weights where appropriate, and an auditable access route.

## Environmental layer decisions

ChEMBL, BindingDB, and PubChem BioAssay are frozen as three different but
operationally usable evidence layers. They must not be merged into a single
undifferentiated chemical–gene edge list:

- ChEMBL supplies curated compound bioactivity and target relationships.
- BindingDB supplies measured ligand–protein affinity records.
- PubChem BioAssay supplies deposited assay activity, with gene/protein targets
  only when an assay defines them.

EPA CompTox/ToxCast and Tox21 remain conditional because they are toxicology
assay activity resources rather than directly comparable chemical–gene sources.
LINCS/L1000 remains conditional as a perturbational transcriptomic-response
layer, not as direct target evidence. Promotion requires endpoint semantics,
human/cell-context rules, identifier mapping, and a release snapshot to be
locked before inspection of candidate coverage.

Official access references: [EPA CompTox downloadable data](https://www.epa.gov/comptox-tools/downloadable-computational-toxicology-data),
[BindingDB REST services](https://www.bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp),
[PubChem BioAssay PUG-REST tutorial](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial), and
[LINCS/CMap data access](https://clue.io/connectopedia/pdf/lincs_cmap_data).

## Disease layer decisions

Open Targets is frozen as an integrated target–disease evidence source, while
the GWAS Catalog is frozen as a human-genetics source. This deliberately
contrasts knowledge-integrated evidence with trait-locus evidence. They will be
queried using source-native disease ontology identifiers and separately
reported evidence fields.

DisGeNET is conditional because academic access is plan-dependent and free
academic access does not provide the full database. OMIM is conditional because
authorized access and release/licensing details must be documented. ClinVar is
conditional because it is a variant-interpretation resource and requires a
prespecified variant-to-gene disease rule before it can generate a comparable
disease-gene set.

Official access references: [Open Targets API](https://api.platform.opentargets.org/),
[GWAS Catalog downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads),
[DisGeNET academic download policy](https://support.disgenet.com/support/solutions/articles/202000089413-can-i-download-the-disgenet-database),
[OMIM](https://www.omim.org/), and [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/).

## Epidemiology layer decisions

CHMS, KNHANES, HBM4EU/IPCHEM, and MESA remain conditional feasibility sources.
Their official documentation establishes relevant health, biomonitoring, or
cohort infrastructure, but Step 10A did not establish the complete joint
availability of the exact exposure, T2D outcome, covariates, weights/design
variables, and usable individual-level access. They must therefore undergo a
source-specific data dictionary and access audit before any can replace NHANES.

HBM4EU/IPCHEM is particularly useful for biomonitoring metadata and exposure
data, but individual-level access is controlled and participant-level T2D
linkage is not assumed. CHMS has exposure and diabetes-related health measures,
but individual microdata access is through Statistics Canada Research Data
Centres. KNHANES provides chronic disease and biomonitoring infrastructure, but
cycle-specific laboratory panels and harmonization remain to be checked. MESA
is a controlled-access cohort and its relevant assay/sample availability must be
confirmed rather than inferred from the general study description.

Official feasibility references: [CHMS access and cycles](https://www.statcan.gc.ca/en/statistical-programs/document/5071_D5_V2),
[KNHANES official overview](https://knhanes.kdca.go.kr/knhanes/eng/main.do),
[HBM4EU/IPCHEM overview](https://ipchem.jrc.ec.europa.eu/hbm4eu_overview.html), and
[MESA at NHLBI](https://www.nhlbi.nih.gov/science/multi-ethnic-study-atherosclerosis-mesa).

NHANES remains `reference_only` because it is the current anchor population.
The previous expanded multi-disease Step10R panel remains
`supplementary_only`; it is a disease-context stress test, not an independent
replacement of the foundational population dataset.

## Independence interpretation

“Independent from CTD/GeneCards/NHANES” in the registry means independent at
the product or data-architecture level. It does **not** guarantee no overlap in
primary publications, source laboratories, chemical structures, or upstream
annotations. Step 10B must quantify overlap and source-component dependence
before treating concordance as robustness.

## Frozen Step 10B gate

Step 10B may begin only after:

1. source-specific release/version and input hashes are captured;
2. chemical, gene, disease, and participant entity-resolution tables are frozen;
3. missingness and unresolved-identifier rules are frozen;
4. comparable evidence units are defined without pooling different evidence
   semantics;
5. at least two environmental replacement sources and both disease replacement
   sources pass their source-specific data audits;
6. an epidemiologic source is promoted from conditional only after its full
   exposure + outcome + covariate + design + access requirements are verified.

Success will be reported using coverage, fixed-top-k retention, overlap,
rank/edge concordance, direction concordance where defined, missingness, and
source-dropout sensitivity. No candidate will be promoted or demoted solely
because a replacement source agrees or disagrees with an existing result.
