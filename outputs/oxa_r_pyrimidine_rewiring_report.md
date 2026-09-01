# OXA-R pyrimidine rewiring

## Question

Does acquired oxaliplatin resistance in CRC show a reproducible de novo pyrimidine/nucleotide-supply transition, providing a rational entry point for DHODH dependency analysis?

## Design

- Datasets: GSE42387 (HCT116, HT29, LoVo parental/OxPt-resistant) and GSE76092 (HT29/HTOXAR3 untreated basal comparison).
- Expression: platform-annotated microarray expression, summarized to gene level by median across probes.
- Unit of analysis: one parental→OxR transition per dataset×cell-line pair; replicate arrays are averaged within state before calculating Δ=OxR−parental.
- Modules: a prespecified three-gene de novo core (CAD/DHODH/UMPS) and a broader pyrimidine/nucleotide-supply panel.

## Available transitions

The analysis retained 4 transitions across 2 datasets. This is a small transition panel and is intended as a screen for a stable direction, not as a powered meta-analysis.

 dataset cell_line    sex  delta_de_novo_pyrimidine_core  delta_pyrimidine_nucleotide_supply  delta_DHODH  delta_RRM2  delta_DCTD  delta_CAD  delta_UMPS  delta_TYMS
GSE42387    HCT116   Male                      -0.039512                           -0.421898     0.227259   -0.412726   -0.259566  -0.258939    0.027634   -0.294863
GSE42387      HT29 Female                       0.130097                            0.652573     0.181119    0.627245    0.075208   0.503969   -0.383150    0.422512
GSE42387      LoVo   Male                       0.685008                            0.541309     0.527958    0.104476    0.278496   0.244181    0.151669   -0.254156
GSE76092      HT29 Female                      -0.085331                            0.292205     0.126667    0.023333    0.096667  -0.120000   -0.040000    0.240000

## Module summary

                      module    sex  n_transitions  mean_delta  median_delta  n_down  n_up  fraction_down
     de_novo_pyrimidine_core Female              2    0.022383      0.022383       1     1           0.50
     de_novo_pyrimidine_core   Male              2    0.322748      0.322748       1     1           0.50
     de_novo_pyrimidine_core    All              4    0.172566      0.045293       2     2           0.50
pyrimidine_nucleotide_supply Female              2    0.472389      0.472389       0     2           0.00
pyrimidine_nucleotide_supply   Male              2    0.059705      0.059705       1     1           0.50
pyrimidine_nucleotide_supply    All              4    0.266047      0.416757       1     3           0.25

## Focal-gene summary

 gene    sex  n_transitions  mean_delta  median_delta  n_down  n_up  fraction_down
  CAD Female              2    0.191985      0.191985       1     1           0.50
  CAD   Male              2   -0.007379     -0.007379       1     1           0.50
  CAD    All              4    0.092303      0.062091       2     2           0.50
 DCTD Female              2    0.085937      0.085937       0     2           0.00
 DCTD   Male              2    0.009465      0.009465       1     1           0.50
 DCTD    All              4    0.047701      0.085937       1     3           0.25
DHODH Female              2    0.153893      0.153893       0     2           0.00
DHODH   Male              2    0.377608      0.377608       0     2           0.00
DHODH    All              4    0.265751      0.204189       0     4           0.00
 RRM2 Female              2    0.325289      0.325289       0     2           0.00
 RRM2   Male              2   -0.154125     -0.154125       1     1           0.50
 RRM2    All              4    0.085582      0.063905       1     3           0.25
 TYMS Female              2    0.331256      0.331256       0     2           0.00
 TYMS   Male              2   -0.274510     -0.274510       2     0           1.00
 TYMS    All              4    0.028373     -0.007078       2     2           0.50
 UMPS Female              2   -0.211575     -0.211575       2     0           1.00
 UMPS   Male              2    0.089651      0.089651       0     2           0.00
 UMPS    All              4   -0.060962     -0.006183       2     2           0.50

## Interpretation

The key criterion for moving to dependency analysis is cross-transition directionality of the de novo core and/or a coherent DHODH-centered response. A mixed module score with isolated gene changes should be treated as rewiring rather than proof of pathway activation.

No CRISPR dependency or pharmacologic vulnerability analysis is included in this first pass. Those will be run only if the rewiring pattern is sufficiently coherent.

## Provenance

{
  "datasets": {
    "GSE42387": {
      "dataset": "GSE42387",
      "n_samples": 27,
      "n_probes": 15,
      "genes_present": [
        "CAD",
        "CMPK1",
        "CTPS2",
        "DCTD",
        "DHFR",
        "DHODH",
        "DUT",
        "PAICS",
        "PPAT",
        "RRM1",
        "RRM2",
        "TK1",
        "TYMS",
        "UCK2",
        "UMPS"
      ],
      "target_genes": [
        "CAD",
        "DHODH",
        "UMPS",
        "CTPS1",
        "CTPS2",
        "PPAT",
        "PAICS",
        "RRM1",
        "RRM2",
        "DCTD",
        "TYMS",
        "DHFR",
        "TK1",
        "DUT",
        "CMPK1",
        "UCK2"
      ],
      "module_genes_present": {
        "de_novo_pyrimidine_core": [
          "CAD",
          "DHODH",
          "UMPS"
        ],
        "pyrimidine_nucleotide_supply": [
          "CAD",
          "DHODH",
          "UMPS",
          "CTPS2",
          "PPAT",
          "PAICS",
          "RRM1",
          "RRM2",
          "DCTD",
          "TYMS",
          "DHFR",
          "TK1",
          "DUT",
          "CMPK1",
          "UCK2"
        ]
      },
      "n_selected_samples": 18
    },
    "GSE76092": {
      "dataset": "GSE76092",
      "n_samples": 18,
      "n_probes": 16,
      "genes_present": [
        "CAD",
        "CMPK1",
        "CTPS1",
        "CTPS2",
        "DCTD",
        "DHFR",
        "DHODH",
        "DUT",
        "PAICS",
        "PPAT",
        "RRM1",
        "RRM2",
        "TK1",
        "TYMS",
        "UCK2",
        "UMPS"
      ],
      "target_genes": [
        "CAD",
        "DHODH",
        "UMPS",
        "CTPS1",
        "CTPS2",
        "PPAT",
        "PAICS",
        "RRM1",
        "RRM2",
        "DCTD",
        "TYMS",
        "DHFR",
        "TK1",
        "DUT",
        "CMPK1",
        "UCK2"
      ],
      "module_genes_present": {
        "de_novo_pyrimidine_core": [
          "CAD",
          "DHODH",
          "UMPS"
        ],
        "pyrimidine_nucleotide_supply": [
          "CAD",
          "DHODH",
          "UMPS",
          "CTPS1",
          "CTPS2",
          "PPAT",
          "PAICS",
          "RRM1",
          "RRM2",
          "DCTD",
          "TYMS",
          "DHFR",
          "TK1",
          "DUT",
          "CMPK1",
          "UCK2"
        ]
      },
      "n_selected_samples": 6
    }
  },
  "modules": {
    "de_novo_pyrimidine_core": [
      "CAD",
      "DHODH",
      "UMPS"
    ],
    "pyrimidine_nucleotide_supply": [
      "CAD",
      "DHODH",
      "UMPS",
      "CTPS1",
      "CTPS2",
      "PPAT",
      "PAICS",
      "RRM1",
      "RRM2",
      "DCTD",
      "TYMS",
      "DHFR",
      "TK1",
      "DUT",
      "CMPK1",
      "UCK2"
    ]
  },
  "focal_genes": [
    "DHODH",
    "RRM2",
    "DCTD",
    "CAD",
    "UMPS",
    "TYMS"
  ]
}
