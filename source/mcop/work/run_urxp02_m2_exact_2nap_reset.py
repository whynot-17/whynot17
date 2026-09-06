from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(r"C:/Users/21634/Documents/Codex/2026-08-22/non/work/whynot17")
M1 = ROOT / "analysis/disease_agnostic_environmental_framework/urxp02_expanded_molecular_universe"
MECH = ROOT / "analysis/disease_agnostic_environmental_framework/urxp02_mechanism_gene_discovery"
OLD_M2 = ROOT / "analysis/disease_agnostic_environmental_framework/urxp02_m2_disease_branch"
OUT = ROOT / "analysis/disease_agnostic_environmental_framework/urxp02_m2_exact_2nap_disease_branch"

EXACT_IDENTITY = "CTD:C028405; PubChem CID:8663; ChEMBL:CHEMBL14126"
EXACT_CHEMICAL = "2-naphthol"
THYROID_ID = "MONDO_0003240"
HYPERTENSION_ID = "HP_0000822"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def yes(v: object) -> bool:
    return str(v).strip().lower() in {"true", "1", "yes"}


def role_for_row(row: pd.Series) -> set[str]:
    note = str(row.get("notes", "")).lower()
    source = str(row.get("source_database", ""))
    evidence_type = str(row.get("evidence_type", "")).lower()
    roles: set[str] = set()
    # A gene/protein explicitly bound, inhibited, or assayed as a single protein.
    if re.search(r"\bbind|binding|inhibit|inhibitory concentration|target_type=single protein|target activity", note):
        roles.add("direct_binding")
    # Preserve the direction of chemical clearance/formation separately from response.
    if re.search(r"metabol|sulfat|glucuron|abundance of 2-naphthol|sulfonation", note):
        roles.add("gene_to_chemical_metabolism")
    if source == "PubChem BioAssay" or "bioassay" in evidence_type:
        roles.add("bioassay_target")
    if source == "ChEMBL" and "target_type=SINGLE PROTEIN" in note and not re.search(r"km|vmax|sulfonation", note):
        roles.add("direct_binding")
    if not roles:
        roles.add("other_exact_interaction")
    direction = str(row.get("direction", "")).strip().lower()
    if direction in {"", "nan", "none", "unspecified", "affects"}:
        roles.add("ambiguous_direction")
    return roles


def aggregate_exact(exact_long: pd.DataFrame, m1_summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    summary = m1_summary.set_index("gene_symbol")
    for gene, sub in exact_long.groupby("gene_symbol", sort=True):
        roles = set().union(*(role_for_row(r) for _, r in sub.iterrows()))
        species = sorted({str(x) for x in sub["species"].dropna() if str(x).strip()})
        refs = sorted({str(x) for x in sub["dataset_or_reference"].dropna() if str(x).strip()})
        sources = sorted({str(x) for x in sub["source_database"].dropna() if str(x).strip()})
        directions = sorted({str(x) for x in sub["direction"].dropna() if str(x).strip() and str(x).lower() != "nan"})
        tiers = sorted({str(x) for x in sub["evidence_tier"].dropna() if str(x).strip()})
        evidence_types = sorted({str(x) for x in sub["evidence_type"].dropna() if str(x).strip()})
        provenance = []
        for _, ev in sub.iterrows():
            provenance.append(
                f"{ev.source_database}|{ev.dataset_or_reference}|{ev.species}|{ev.direction}|{ev.evidence_type}"
            )
        # These two fields must reproduce M1b's exact-human / exact-experimental definition.
        m1 = summary.loc[gene]
        exact_human = bool(m1["exact_2NAP_human_support"])
        exact_experimental = bool(m1["exact_2NAP_experimental_support"])
        calculated_human = bool((sub["species"] == "Homo sapiens").any())
        calculated_experimental = bool((sub["species"] != "Homo sapiens").any())
        if exact_human != calculated_human or exact_experimental != calculated_experimental:
            raise AssertionError(f"M1b evidence flag mismatch for {gene}")
        rows.append({
            "gene_symbol": gene,
            "gene_id": ";".join(sorted({str(x) for x in sub["gene_id"].dropna() if str(x).strip()})),
            "chemical": EXACT_CHEMICAL,
            "chemical_identity": EXACT_IDENTITY,
            "exact_evidence_rows": len(sub),
            "exact_human_supported": exact_human,
            "exact_experimental_supported": exact_experimental,
            "CTD_support": bool((sub["source_database"] == "CTD").any()),
            "PubChem_support": bool((sub["source_database"] == "PubChem BioAssay").any()),
            "ChEMBL_support": bool((sub["source_database"] == "ChEMBL").any()),
            "number_independent_sources": len(sources),
            "number_independent_references": len(refs),
            "number_species": len(species),
            "evidence_species": ";".join(species),
            "evidence_tiers": ";".join(tiers),
            "evidence_types": ";".join(evidence_types),
            "evidence_directions": ";".join(directions),
            "evidence_role": ";".join(sorted(roles)),
            "direct_target_binding_flag": "direct_binding" in roles,
            "metabolism_related_flag": "gene_to_chemical_metabolism" in roles,
            "ambiguous_direction_flag": "ambiguous_direction" in roles,
            "evidence_source_databases": ";".join(sources),
            "evidence_provenance": ";".join(provenance),
            "parent_naphthalene_rows_in_primary": 0,
            "primary_universe_rule": "exact 2-NAP evidence required; parent-only evidence excluded",
        })
    out = pd.DataFrame(rows)
    if len(out) != 110 or out.gene_symbol.nunique() != 110:
        raise AssertionError(f"Expected exactly 110 unique exact genes, got {len(out)}")
    return out


def add_disease_fields(exact: pd.DataFrame, thyroid: pd.DataFrame, hyper: pd.DataFrame) -> pd.DataFrame:
    tmap = thyroid.set_index("gene_symbol").to_dict("index")
    hmap = hyper.set_index("gene_symbol").to_dict("index")
    rows = []
    for _, r in exact.iterrows():
        gene = r.gene_symbol
        t = tmap.get(gene)
        h = hmap.get(gene)
        in_t = t is not None
        in_h = h is not None
        branch = "shared" if in_t and in_h else "thyroid-specific" if in_t else "hypertension-specific" if in_h else "neither"
        row = r.to_dict()
        row.update({
            "in_thyroid_disease_set": in_t,
            "in_hypertension_disease_set": in_h,
            "branch_class": branch,
            "thyroid_disease_id": t.get("disease_id", THYROID_ID) if t else "",
            "thyroid_disease_name": t.get("disease_name", "thyroid gland disorder") if t else "",
            "thyroid_disease_relevance_score": t.get("disease_relevance_score", "") if t else "",
            "thyroid_disease_rank": t.get("rank", "") if t else "",
            "thyroid_source": t.get("source", "") if t else "",
            "thyroid_association_scope": t.get("association_scope", "") if t else "",
            "hypertension_disease_id": h.get("disease_id", HYPERTENSION_ID) if h else "",
            "hypertension_disease_name": h.get("disease_name", "Hypertension") if h else "",
            "hypertension_disease_relevance_score": h.get("disease_relevance_score", "") if h else "",
            "hypertension_disease_rank": h.get("rank", "") if h else "",
            "hypertension_source": h.get("source", "") if h else "",
            "hypertension_association_scope": h.get("association_scope", "") if h else "",
            "branch_evidence_note": "Disease membership is direct Open Targets symbol membership; it does not imply causality, target status, or sex specificity. Parent naphthalene-only evidence is excluded from the primary universe.",
        })
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    unified = pd.read_csv(M1 / "04_unified_gene_evidence_long.csv")
    m1_summary = pd.read_csv(M1 / "05_gene_evidence_summary.csv")
    exact_long = unified.loc[
        unified["chemical"].eq(EXACT_CHEMICAL) & unified["chemical_identity"].eq(EXACT_IDENTITY)
    ].copy()
    if exact_long.gene_symbol.nunique() != 110:
        raise AssertionError("Exact identity filter did not produce 110 unique genes")
    if len(exact_long) != 524:
        raise AssertionError(f"Expected 524 exact evidence rows, got {len(exact_long)}")
    exact = aggregate_exact(exact_long, m1_summary)

    thyroid = pd.read_csv(MECH / "02_thyroid_disease_genes.csv")
    hypertension = pd.read_csv(MECH / "03_hypertension_disease_genes.csv")
    classified = add_disease_fields(exact, thyroid, hypertension)
    counts = classified["branch_class"].value_counts().to_dict()
    expected = {"thyroid-specific": 0, "hypertension-specific": 29, "shared": 34, "neither": 47}
    observed = {branch: int(counts.get(branch, 0)) for branch in expected}
    if observed != expected:
        raise AssertionError(f"Unexpected exact branch counts: {observed}")
    if sum(observed.values()) != 110:
        raise AssertionError("Branch identity check failed")
    counts = observed

    base_cols = list(exact.columns)
    exact.to_csv(OUT / "01_exact_2NAP_110_gene_universe.csv", index=False)
    intersection_cols = base_cols + [
        "in_thyroid_disease_set", "in_hypertension_disease_set", "branch_class",
        "thyroid_disease_id", "thyroid_disease_name", "thyroid_disease_relevance_score", "thyroid_disease_rank", "thyroid_source", "thyroid_association_scope",
        "hypertension_disease_id", "hypertension_disease_name", "hypertension_disease_relevance_score", "hypertension_disease_rank", "hypertension_source", "hypertension_association_scope", "branch_evidence_note",
    ]
    classified.loc[classified.in_thyroid_disease_set, intersection_cols].sort_values(["branch_class", "gene_symbol"]).to_csv(OUT / "02_exact_2NAP_thyroid_intersection.csv", index=False)
    classified.loc[classified.in_hypertension_disease_set, intersection_cols].sort_values(["branch_class", "gene_symbol"]).to_csv(OUT / "03_exact_2NAP_hypertension_intersection.csv", index=False)
    classified.sort_values("gene_symbol")[intersection_cols].to_csv(OUT / "04_exact_2NAP_branch_classification.csv", index=False)

    summary_rows = []
    for branch in ["thyroid-specific", "hypertension-specific", "shared", "neither"]:
        sub = classified.loc[classified.branch_class.eq(branch)]
        summary_rows.append({
            "branch_class": branch,
            "n_genes": len(sub),
            "n_exact_human_supported": int(sub.exact_human_supported.sum()),
            "n_exact_experimental_supported": int(sub.exact_experimental_supported.sum()),
            "n_CTD_supported": int(sub.CTD_support.sum()),
            "n_PubChem_supported": int(sub.PubChem_support.sum()),
            "n_ChEMBL_supported": int(sub.ChEMBL_support.sum()),
            "n_ge2_independent_sources": int((sub.number_independent_sources >= 2).sum()),
            "n_direct_target_binding": int(sub.direct_target_binding_flag.sum()),
            "n_metabolism_related": int(sub.metabolism_related_flag.sum()),
            "n_chemical_to_gene_response_or_interaction": int((sub.evidence_role.str.contains("other_exact_interaction|direct_binding|bioassay_target")).sum()),
            "n_gene_to_chemical_metabolism_role": int(sub.evidence_role.str.contains("gene_to_chemical_metabolism").sum()),
            "n_ambiguous_direction": int(sub.ambiguous_direction_flag.sum()),
            "evidence_note": "Exact 2-NAP evidence only; parent naphthalene-only rows are excluded. Metabolism-related genes are not treated as equivalent to genes whose activity/expression changes in response to 2-NAP.",
        })
    pd.DataFrame(summary_rows).to_csv(OUT / "05_exact_2NAP_branch_evidence_summary.csv", index=False)

    old_class = pd.read_csv(OLD_M2 / "03_branch_gene_classification.csv")
    old_counts = old_class.branch_class.value_counts().to_dict()
    if old_counts != {"neither": 358, "shared": 189, "hypertension-specific": 251, "thyroid-specific": 30}:
        raise AssertionError(f"Old M2 branch counts changed: {old_counts}")
    exact_set = set(exact.gene_symbol)
    old_sets = {b: set(old_class.loc[old_class.branch_class.eq(b), "gene_symbol"]) for b in old_counts}
    exact_sets = {b: set(classified.loc[classified.branch_class.eq(b), "gene_symbol"]) for b in expected}
    comparisons = []
    def add(name, old_category, old_n, exact_category, exact_n, overlap, notes, parent_only=0):
        comparisons.append({
            "comparison": name,
            "old_category": old_category,
            "old_n": old_n,
            "exact_category": exact_category,
            "exact_n": exact_n,
            "overlap_n": overlap,
            "removed_n": old_n - overlap,
            "survival_percent": round(100 * overlap / old_n, 4) if old_n else 0.0,
            "removed_percent": round(100 * (old_n - overlap) / old_n, 4) if old_n else 0.0,
            "old_shared_parent_only_n": parent_only,
            "notes": notes,
        })
    add("overall universe", "old expanded universe", 828, "exact 2-NAP universe", 110, len(exact_set), "Exact universe is a hard identity filter over M1b; parent-only genes are excluded.")
    add("thyroid intersection", "old 828 ∩ thyroid", 219, "exact 110 ∩ thyroid", len(exact_set & (old_sets["shared"] | old_sets["thyroid-specific"])), len(exact_set & old_sets["shared"] | exact_set & old_sets["thyroid-specific"]), "Survival of the old thyroid disease intersection under exact restriction.")
    add("hypertension intersection", "old 828 ∩ hypertension", 440, "exact 110 ∩ hypertension", len(exact_set & (old_sets["shared"] | old_sets["hypertension-specific"])), len(exact_set & old_sets["shared"] | exact_set & old_sets["hypertension-specific"]), "Survival of the old hypertension disease intersection under exact restriction.")
    add("thyroid-specific branch", "old thyroid-specific", 30, "exact thyroid-specific", 0, len(exact_sets["thyroid-specific"]), "No exact 2-NAP gene remains thyroid-specific under the frozen disease definitions.")
    add("hypertension-specific branch", "old hypertension-specific", 251, "exact hypertension-specific", 29, len(exact_sets["hypertension-specific"]), "Exact restriction retains a small hypertension-specific branch.")
    add("shared core", "old shared", 189, "exact shared", 34, len(exact_set & old_sets["shared"]), "The old shared core is reduced to 34 genes; 155/189 are not exact 2-NAP genes and are parent-expanded removals.", parent_only=len(old_sets["shared"] - exact_set))
    add("neither branch", "old neither", 358, "exact neither", 47, len(exact_set & old_sets["neither"]), "Exact genes outside both frozen disease lists.")
    # Include an explicit identity audit row for the required 110-gene partition.
    add("exact partition identity", "not applicable", 110, "thyroid-specific + hypertension-specific + shared + neither", 110, 110, "0 + 29 + 34 + 47 = 110; mandatory identity check passed.")
    pd.DataFrame(comparisons).to_csv(OUT / "06_old_vs_exact_M2_comparison.csv", index=False)

    report = f"""# M2-RESET — exact 2-NAP disease branching

Generated UTC: {datetime.now(timezone.utc).isoformat()}

## Primary rule

> **Parent naphthalene-only evidence is excluded from the primary molecular universe.**

The primary universe is constructed only from M1b rows with chemical `2-naphthol` and the exact identity `CTD:C028405; PubChem CID:8663; ChEMBL:CHEMBL14126` (2-naphthol / 2-hydroxynaphthalene / 2-NAP). The filter yields 524 exact evidence rows and exactly **110 unique gene symbols**. No gene enters because of parent naphthalene evidence alone.

The 110 symbols are retained as supplied by M1b, including source labels that may not be canonical HGNC symbols. They are evidence-linked entities, not 110 proven causal targets of 2-NAP.

## Disease definitions

The original M2 definitions are frozen and reused without broadening:

- Thyroid disease: `{THYROID_ID}` (`enableIndirect=false` Open Targets direct associations).
- Hypertension: `{HYPERTENSION_ID}` (`enableIndirect=false` Open Targets direct associations).

No thyroid cancer, thyroid hormone traits, blood-pressure GWAS traits, cardiovascular surrogates, or kidney surrogates were introduced.

## Exact-2-NAP branch counts

| Branch | Genes |
|---|---:|
| Thyroid-specific (A − B) | **{counts['thyroid-specific']}** |
| Hypertension-specific (B − A) | **{counts['hypertension-specific']}** |
| Shared (A ∩ B) | **{counts['shared']}** |
| Neither | **{counts['neither']}** |
| **Total** | **110** |

The mandatory partition identity passes: **0 + 29 + 34 + 47 = 110**.

There is no exact-2-NAP thyroid-specific gene under these frozen disease lists. The exact universe therefore has a hypertension-specific component (29), a shared component (34), and a sizeable neither component (47), but no exact thyroid-only branch.

## Evidence characterization

Evidence flags in the branch summary are recomputed from exact rows only. `exact_experimental_supported` follows the M1b convention of any exact row whose species is not `Homo sapiens` (including source rows with unknown species labels). Evidence roles distinguish `gene_to_chemical_metabolism` from `direct_binding`, `bioassay_target`, and other exact interactions. A metabolism enzyme that changes 2-NAP abundance is not treated as a gene whose expression/activity is changed by 2-NAP.

These records preserve provenance and direction where supplied; they do not establish causality. In particular, this analysis does not call all 110 genes direct targets.

## Comparison with the old parent-expanded 828-gene M2

- Old shared core: **189 → 34 exact shared genes**; **155/189 (82.01%)** are removed under the exact identity restriction.
- Old thyroid intersection: 219 → 34 exact genes (15.53% retained).
- Old hypertension intersection: 440 → 63 exact genes (14.32% retained).
- Old thyroid-specific branch: 30 → **0** exact genes.
- Old hypertension-specific branch: 251 → 29 exact genes (11.55% retained).

Therefore the prior broad 828-gene shared-core interpretation **does not remain valid as the primary exact-2-NAP interpretation**. A 34-gene exact shared subset remains, but the broad 189-gene parent-expanded shared core has collapsed and the exact branch structure is not reciprocal thyroid-vs-hypertension.

## Stop rule

This reset stops after exact-universe construction, disease intersections, branch classification, evidence characterization, and old-vs-exact audit. No enrichment, GO/KEGG/Reactome, PPI, STRING, modules, GTEx, sex-DE, single-cell analysis, immune infiltration, figures, or new NHANES analysis was run.
"""
    (OUT / "M2_EXACT_2NAP_DISEASE_BRANCH_REPORT.md").write_text(report, encoding="utf-8")

    inputs = [
        M1 / "04_unified_gene_evidence_long.csv",
        M1 / "05_gene_evidence_summary.csv",
        MECH / "02_thyroid_disease_genes.csv",
        MECH / "03_hypertension_disease_genes.csv",
        OLD_M2 / "03_branch_gene_classification.csv",
        OLD_M2 / "manifest.json",
    ]
    output_paths = [
        OUT / "01_exact_2NAP_110_gene_universe.csv",
        OUT / "02_exact_2NAP_thyroid_intersection.csv",
        OUT / "03_exact_2NAP_hypertension_intersection.csv",
        OUT / "04_exact_2NAP_branch_classification.csv",
        OUT / "05_exact_2NAP_branch_evidence_summary.csv",
        OUT / "06_old_vs_exact_M2_comparison.csv",
        OUT / "M2_EXACT_2NAP_DISEASE_BRANCH_REPORT.md",
    ]
    manifest = {
        "analysis": "URXP02 M2-RESET exact 2-NAP disease branching",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "primary_universe": {
            "chemical": EXACT_CHEMICAL,
            "ctd_id": "C028405",
            "pubchem_cid": 8663,
            "chembl_id": "CHEMBL14126",
            "exact_identity": EXACT_IDENTITY,
            "evidence_rows": int(len(exact_long)),
            "unique_gene_symbols": int(exact.gene_symbol.nunique()),
            "parent_naphthalene_only_included": False,
        },
        "disease_definitions": {
            "thyroid": {"id": THYROID_ID, "source": "Open Targets direct associations", "enableIndirect": False},
            "hypertension": {"id": HYPERTENSION_ID, "source": "Open Targets direct associations", "enableIndirect": False},
        },
        "counts": {"exact_universe": 110, **{k.replace("-", "_"): int(v) for k, v in counts.items()}},
        "mandatory_identity_check": "0 + 29 + 34 + 47 = 110",
        "old_m2_counts": {k.replace("-", "_"): int(v) for k, v in old_counts.items()},
        "inputs_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
        "outputs_sha256": {p.name: sha256(p) for p in output_paths},
        "constraints": [
            "Parent naphthalene-only evidence excluded from primary universe",
            "No disease broadening",
            "No enrichment",
            "No PPI/STRING/modules",
            "No GTEx or sex-DE",
            "No single-cell analysis",
            "No figures",
            "No new NHANES analysis",
            "No causal interpretation",
        ],
        "output_files": [p.name for p in output_paths] + ["manifest.json"],
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
