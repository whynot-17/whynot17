# Phase 2 — PubMed literature collision audit

## Scope

- **Search date (UTC):** `2026-08-28`
- **Candidates:** 15 chemical IDs in the Phase 1 advance pool, grouped into 14 search groups because the two DINP parent IDs share the same alias set.
- **Database:** PubMed via NCBI E-utilities `esearch.fcgi`.
- **Search categories:** diabetes-related total, human epidemiology, prospective/longitudinal, mechanism, animal/cell, target/pathway, and network/bioinformatics.
- **Interpretation boundary:** retrieval counts are screening signals, not counts of eligible studies. Every candidate remains subject to title/abstract and, where needed, full-text adjudication.

## Search design

Each query was constructed as `(candidate alias set) AND (category terms)`. The complete query, result count, top relevance-sorted PMIDs, status, and source URL are retained in `pubmed_collision_counts.csv`. DINP parent IDs `C012125` and `C019174` deliberately share the `dinp_parent` search group to avoid duplicate literature counts while preserving both chemical identities.

## Outputs

- `03_literature_collision/pubmed_collision_counts.csv`: one row per search group × category, 98 queries in total.
- `03_literature_collision/literature_counts.csv`: one row per chemical ID with category-level retrieval counts and conservative screening flags.
- `03_literature_collision/candidate_search_vocabulary.csv`: chemical IDs, names, search groups, and aliases.

## Important limitation

The category filters are intentionally sensitive and can retrieve papers that mention an analyte, a mixture, a therapeutic formulation, or a non-T2D metabolic endpoint. Therefore `human_epidemiology_pubmed_hits`, `prospective_pubmed_hits`, and `mechanism_pubmed_hits` are not yet adjudicated evidence counts. A later screening pass must verify whether the exposure is the candidate chemical, whether the outcome is T2D/diabetes, study design, and whether a complete chemical→T2D mechanism exists.

## Phase 2 decision status

This audit establishes the reproducible collision-search layer. It does **not** yet assign final novelty grades, docking feasibility, target status, or a definitive Top 5. Those require review of the retrieved records rather than ranking candidates from raw search counts alone.
