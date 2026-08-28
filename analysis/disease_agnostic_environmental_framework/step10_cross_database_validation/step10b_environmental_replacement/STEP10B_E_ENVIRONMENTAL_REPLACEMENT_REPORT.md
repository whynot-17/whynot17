# Step 10B-E — Environmental knowledge-source replacement audit

Generated: `2026-08-28T09:41:05.542742+00:00`
Candidates analyzed: **134** (inherited frozen Step 4/7 universe)

ChEMBL measured bioactivity, BindingDB affinity, and PubChem BioAssay compound-to-assay membership were queried as separate evidence layers. They were not merged into a chemical-gene edge list.

Missing or unresolved records are reported as not observed/unresolved and are not interpreted as biological negatives. PubChem CID-to-AID membership does not establish a human protein target, so a PubChem target count is not fabricated.

Exact source metadata, query rules, API response hashes, input/output hashes, and call records are in `STEP10B_E_SOURCE_SNAPSHOT.json` and `STEP10B_E_API_CALL_MANIFEST.json`.

Status: complete as a source-specific coverage/replacement audit; mechanistic target inference requires additional source-native assay parsing.