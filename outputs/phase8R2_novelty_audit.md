# Phase 8-R2：post-ranking novelty audit

Audit date: 2026-08-22

## Candidate entering the audit

Only one compound passed the phenotype and identity filters as an approved non-oncology candidate:

| Drug | External replication | Platform identity | Initial disposition |
|---|---|---|---|
| Ciclopirox | PRISM + CTRPv2, Level1 independent external platforms | Launched; infectious disease; onychomycosis | Novelty audit required |

## Targeted literature audit

Queries were run after the biological ranking was frozen:

- `ciclopirox colorectal cancer`
- `ciclopirox oxaliplatin colorectal cancer resistance`
- `ciclopirox OXA-resistant colorectal`
- `ciclopirox FOLFOX colorectal`

The key direct hit is Qi et al., *Cell Death & Disease* 2020, PMID [32719342](https://pubmed.ncbi.nlm.nih.gov/32719342/). It studies ciclopirox in CRC models including HCT-8, DLD-1 and HCT-8/5-FU, and reports ROS-mediated PERK-dependent ER stress and cell death in both chemosensitive and 5-FU-resistant CRC cells. The associated open-access record is [PMC7385140](https://pmc.ncbi.nlm.nih.gov/articles/PMC7385140/).

Earlier work also reported antitumor activity of ciclopirox in colon adenocarcinoma cells ([PMC2888914](https://pmc.ncbi.nlm.nih.gov/articles/PMC2888914/)), and a later review explicitly catalogs CRC and colorectal-tumor evidence ([PMID 33573561](https://pubmed.ncbi.nlm.nih.gov/33573561/)).

The targeted exact-name searches did not identify a primary paper specifically testing ciclopirox as an oxaliplatin-resistance reversal agent in CRC. That is an OXA-specific gap, but it does not restore the broader Drug X novelty: ciclopirox has already been directly studied in CRC and in a chemotherapy-resistant CRC model, with a mechanism centered on ROS/ER stress/proteostasis that overlaps the present state architecture.

## Classification

- CRC evidence: **directly established**.
- Chemotherapy-resistant CRC evidence: **directly established for 5-FU resistance**.
- OXA-resistant CRC evidence: **not identified in the targeted audit**.
- OXA-R Drug X novelty: **Class B / downgrade** — CRC and chemotherapy-resistance work already exists; an OXA-specific extension would be incremental rather than a clean new repurposing discovery.

## Decision

`Ciclopirox` is retained as an external pharmacogenomic positive comparator, but is **not promoted as the novel Drug X**. The Phase 8-R2 cross-sectional pharmacogenomic search therefore yields:

> 0 novel approved non-oncology Drug X candidates after independent external replication and novelty audit.

The next computational step should not be another broad baseline pharmacogenomic search. The project should pivot to matched parental/OXA-R screens or use this result to design a focused wet-lab mini-screen. Any ciclopirox follow-up would need to be framed explicitly as an OXA-R-specific extension of an existing CRC/5-FU-resistance mechanism, not as a first report of ciclopirox in CRC.
