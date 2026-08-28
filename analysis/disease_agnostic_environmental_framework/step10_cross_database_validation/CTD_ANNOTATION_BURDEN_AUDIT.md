# CTD annotation-burden audit

## Status

`PASS — descriptive audit complete; external source replacement pending`

This audit tests whether the post-firewall T2D candidate set is simply enriched
for highly annotated CTD chemicals. It does not re-rank chemicals and does not
alter any frozen epidemiologic, FDR, robustness, or pathway result.

## Inputs and scope

The analysis used the formal actionability ledger as the reference population:

| Population | Chemicals |
|---|---:|
| Formal actionability universe | 2,042 |
| Actionable chemicals | 409 |
| Post-firewall T2D candidate chemicals | 134 |
| T2D candidates that are actionable in the ledger | 134 |

CTD chemical–gene records were restricted to human records (`OrganismID=9606`
or `Organism=Homo sapiens`). Counts were retained at separate levels:

- raw human interaction rows;
- unique chemical × gene pairs;
- unique human genes;
- unique PMIDs;
- unique interaction actions and interaction terms;
- CTD vocabulary parents, tree numbers, synonyms, and definition length.

The CTD interaction source and vocabulary source are the local files listed in
`CTD_ANNOTATION_BURDEN_AUDIT_SUMMARY.json`; their SHA-256 hashes are recorded in
the summary and manifest.

## Main findings

Only 31/134 (23.1%) T2D candidate chemicals had at least one human CTD gene
annotation in the current file, compared with 124/409 (30.3%) actionable
chemicals. Thus, the candidate set was not uniformly annotation-rich.

The actionable-background 90th percentile was 56 unique human genes and 11.2
unique PMIDs. Ten of 134 T2D candidates fell in the top gene-burden decile and
12/134 fell in the top PMID-burden decile. These counts are descriptive and are
not used as filters.

The candidate and actionable noncandidate distributions were highly
zero-inflated and long-tailed: both groups had a median of zero for the main
burden measures. The recorded two-sided Mann–Whitney tests were nominally
different for unique genes (`P=0.0299`) and PMIDs (`P=0.0448`), but these values
do not establish enrichment or explain candidate selection. The one-sided
Fisher comparison for candidate membership in the actionable top gene-burden
decile gave `OR=0.612`, `P=0.9345`, which is not evidence of enrichment.

Across the full actionability universe, unique human gene count and unique PMID
count were almost perfectly correlated (`Spearman rho=0.9992`). This is an
important audit result: CTD interaction breadth is strongly coupled to study
volume, so raw CTD degree cannot be interpreted as independent biological
support.

The largest burdens in the candidate set were concentrated in a small number
of highly studied entries, including lead, diethylhexyl phthalate, silver,
perfluorohexanesulfonic acid, and uranium. DINP-related entries were not the
largest CTD annotation outliers: the parent DINP entry had 86 unique human genes
and 16 PMIDs, whereas the mono(carboxy-isooctyl)phthalate entry had 19 genes and
2 PMIDs. These values are audit descriptors only.

## Interpretation

The current data do **not** support the simple explanation that the T2D
candidate set consists mainly of the most heavily annotated CTD chemicals. They
also show why annotation burden must be measured: a minority of chemicals carry
very large CTD interaction counts, and those counts are largely explained by
PMID/study volume.

This is not a clean bill of health for CTD. The audit does not test every form of
CTD disease/pathway annotation bias, and a zero CTD count means “not represented
in this source,” not “no biological relevance.” The next test is therefore
source replacement, not a new CTD-derived score.

## Reproducible outputs

- `compute_ctd_annotation_burden_audit.py`
- `ctd_annotation_burden_by_chemical.csv`
- `ctd_annotation_burden_candidate_summary.csv`
- `CTD_ANNOTATION_BURDEN_AUDIT_SUMMARY.json`
- `CTD_ANNOTATION_BURDEN_AUDIT_MANIFEST.json`

