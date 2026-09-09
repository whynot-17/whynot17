from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent


def read_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def first(rows, gene):
    return next((row for row in rows if row.get("gene_symbol") == gene), {})


def val(row, key, default=""):
    return row.get(key, default) if row else default


def write_csv(path: Path, rows):
    fields = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    overlap = read_csv(ROOT / "outputs" / "DINP_CRC_overlap.csv")
    target_master = read_csv(ROOT / "outputs" / "DINP_target_master.csv")
    cross_rank = read_csv(ROOT / "outputs" / "DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv")
    consensus = read_csv(ROOT / "outputs" / "DINP_CRC_STRING_cytohubba_8_consensus_hubs.csv")

    # The shortlist is intentionally small. PPARG is the only proposed full
    # docking+MD target; the other rows are docking support/negative-control
    # candidates and are not substitutes for the mechanistic bridge.
    plans = [
        {
            "selection_rank": 1,
            "gene_symbol": "PPARG",
            "analysis_role": "primary_docking_and_100ns_md",
            "pocket_feasibility": "high",
            "representative_pdb": "4XLD; 2PRG",
            "structure_note": "Human PPARgamma LBD holo structure with rosiglitazone (4XLD) plus apo LBD (2PRG); use the holo pocket for docking and apo/alternative holo structure for pose sensitivity.",
            "dinp_to_cebp_bridge": "Strongest available bridge: DINP-specific preadipocyte exposure activated PPARgamma and increased early Cebpb; PPARG is a 97-gene overlap and cross-supported candidate.",
            "decision_note": "Primary target. Docking is hypothesis-generating because the current DINP record is curated/predicted plus functional exposure evidence, not purified-protein direct binding.",
        },
        {
            "selection_rank": 2,
            "gene_symbol": "PPARA",
            "analysis_role": "support_docking",
            "pocket_feasibility": "high",
            "representative_pdb": "6KB1; 3VI8",
            "structure_note": "Human PPARalpha LBD co-crystal structures with ligands; compact nuclear-receptor orthosteric pocket.",
            "dinp_to_cebp_bridge": "Same PPAR family and DINP-related curated/predicted binding plus DINP exposure response; the CEBPB link is weaker and less direct than PPARG.",
            "decision_note": "Support docking candidate; useful for testing PPAR-family selectivity and pose consistency.",
        },
        {
            "selection_rank": 3,
            "gene_symbol": "RXRA",
            "analysis_role": "support_docking",
            "pocket_feasibility": "high",
            "representative_pdb": "1FBY; 1FM6",
            "structure_note": "Human RXRalpha LBD holo structure (1FBY) and PPARgamma/RXRalpha heterodimer (1FM6); both provide a defined nuclear-receptor ligand pocket/context.",
            "dinp_to_cebp_bridge": "RXRA is the obligate heterodimer partner for PPAR nuclear receptors; DINP–RXRA evidence is currently curated/predicted rather than direct purified-protein binding.",
            "decision_note": "Support docking candidate; interpret as receptor-complex context, not proof that RXRA is the primary DINP sensor.",
        },
        {
            "selection_rank": 4,
            "gene_symbol": "NR1I2",
            "analysis_role": "support_docking",
            "pocket_feasibility": "high",
            "representative_pdb": "8CH8",
            "structure_note": "Human PXR/NR1I2 LBD holo structure with liranaftate; the single-chain receptor LBD provides a defined pocket for support docking.",
            "dinp_to_cebp_bridge": "DINP-specific human receptor transactivation/ex vivo evidence exists, but the connection to CEBPB is indirect and not CRC-specific.",
            "decision_note": "Support docking candidate because the DINP receptor activity evidence is stronger than for most non-PPAR candidates; do not overstate a CEBPB mechanism.",
        },
        {
            "selection_rank": 5,
            "gene_symbol": "ESR1",
            "analysis_role": "support_docking_crc_network_control",
            "pocket_feasibility": "high",
            "representative_pdb": "1XPC",
            "structure_note": "Human estrogen-receptor alpha LBD structures with defined ligand pocket; verify the exact released entry and ligand state before docking.",
            "dinp_to_cebp_bridge": "Strong CRC network support (hub/consensus/cross-rank), but no DINP-specific direct-binding evidence in the current master table and weaker DINP-to-CEBPB mechanistic continuity.",
            "decision_note": "Optional support/control docking only; it should not displace PPARG as the primary target.",
        },
    ]

    rows = []
    for plan in plans:
        gene = plan["gene_symbol"]
        ov = first(overlap, gene)
        tm = first(target_master, gene)
        cr = first(cross_rank, gene)
        cc = first(consensus, gene)
        row = dict(plan)
        row.update(
            {
                "in_97_overlap": "yes" if ov else "no",
                "dinp_evidence_grade_current": val(tm, "evidence_grade"),
                "dinp_database_source": val(tm, "database_source"),
                "dinp_evidence_type": val(tm, "evidence_type"),
                "direct_binding_or_activity_audit": val(tm, "direct_binding_or_activity"),
                "dinp_species": val(tm, "species"),
                "dinp_evidence_note": val(tm, "evidence_note"),
                "cross_rank": val(cr, "cross_rank"),
                "candidate_tier": val(cr, "candidate_tier"),
                "consensus_rank": val(cc, "consensus_rank"),
                "stable_ml_flag": val(cc, "stable_ml_flag"),
                "tier1_flag": val(cc, "tier1"),
                "docking_status": "completed_vina_screen",
                "md_status": "deferred_to_external_machine_by_design",
            }
        )
        rows.append(row)

    write_csv(OUT / "DINP_CRC_docking_target_shortlist.csv", rows)

    exclusions = [
        ("CEBPB", "downstream mediator; intrinsically/disorder-rich transcription factor and no defensible primary ligand pocket selected"),
        ("IL6/TNF/IL1B", "high network centrality but cytokine proteins are not clean small-molecule orthosteric targets for this question"),
        ("CD36", "membrane scavenger receptor with a more complex extracellular/membrane system and weaker direct DINP-pocket rationale"),
        ("AKT1/STAT3", "useful network hubs, but no specific DINP-to-CEBPB direct ligand bridge in the current evidence table"),
        ("NR1I3", "DINP-related receptor evidence exists, but the available CAR/RXR structural context is heterodimeric and the gene is not in the 97 overlap; retain as backup rather than primary support"),
    ]

    report = "# DINP–CRC docking target selection\n\n"
    report += f"Generated: {datetime.now(timezone.utc).isoformat()}\n\n"
    report += "## Decision\n\n"
    report += "**PPARG is the proposed primary docking target.** CEBPB remains a downstream core mediator and is explicitly excluded from docking. The primary calculation should use a human PPARG ligand-binding-domain structure with a validated holo pocket; docking scores are not direct-binding proof.\n\n"
    report += "## Shortlist\n\n"
    report += "| Rank | Target | Role | 97 overlap | Pocket | Representative structure(s) |\n|---:|---|---|---|---|---|\n"
    for row in rows:
        report += f"| {row['selection_rank']} | {row['gene_symbol']} | {row['analysis_role']} | {row['in_97_overlap']} | {row['pocket_feasibility']} | {row['representative_pdb']} |\n"
    report += "\n"
    report += "## Why PPARG\n\n"
    report += "PPARG is the cleanest bridge currently present in the project: it is in the 97 DINP–CRC overlap, is cross-supported by PPI/pathway/transcriptomic prioritization, has a well-characterized nuclear-receptor ligand-binding pocket, and the DINP-specific exposure evidence reports PPARG activity with early CEBPB upregulation. The important boundary is that the current DINP evidence is curated/predicted binding plus functional exposure evidence; a docking pose can support structural plausibility but cannot establish biochemical affinity.\n\n"
    report += "## Calculation plan\n\n"
    report += "1. Prepare PPARG LBD and define the grid from the co-crystallized ligand pocket; use a second PPARG structure for pose sensitivity.\n"
    report += "2. Dock a documented representative DINP structure. DINP is a commercial mixture/branched-isomer class, so the exact PubChem/CAS-defined representative and its 3-D conformer must be recorded; PubChem reports a flexible molecule and no deposited 3-D conformer.\n"
    report += "3. Run the full 100 ns explicit-solvent protein–DINP MD only for PPARG after pose and parameterization QC. The support candidates receive docking only.\n"
    report += "4. Report docking score, pose contacts, pocket residues, ligand RMSD/occupancy, protein RMSD/RMSF, radius of gyration, SASA, hydrogen bonds, and replicate/structure sensitivity.\n\n"
    report += "## Explicit exclusions\n\n"
    for gene, reason in exclusions:
        report += f"- **{gene}**: {reason}.\n"
    report += "\n## Current execution status\n\n"
    report += "Target selection and Vina docking are complete: PPARG was run with three seeds, and PPARA/RXRA/NR1I2/ESR1 were run as support screens. MD is intentionally deferred to the user's external machine and is not part of this docking-only package.\n\n"
    report += "## Source records\n\n"
    report += "- Project evidence: `outputs/DINP_CRC_overlap.csv`, `outputs/DINP_target_master.csv`, `outputs/DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv`, `outputs/DINP_CRC_STRING_cytohubba_8_consensus_hubs.csv`.\n"
    report += "- RCSB PDB: [4XLD](https://www.rcsb.org/structure/4XLD), [2PRG](https://www.rcsb.org/structure/2PRG), [6KB1](https://www.rcsb.org/structure/6KB1), [1FBY](https://www.rcsb.org/structure/1FBY), [1FM6](https://www.rcsb.org/structure/1FM6), [6NX1](https://www.rcsb.org/structure/6NX1), [1XPC](https://www.rcsb.org/structure/1XPC).\n"
    report += "- PubChem DINP: [CID 590836](https://pubchem.ncbi.nlm.nih.gov/compound/Diisononyl-phthalate), CAS 28553-12-0 / 68515-48-0.\n"
    (OUT / "DINP_CRC_docking_target_selection_report.md").write_text(report, encoding="utf-8")

    source_files = [
        ROOT / "outputs" / "DINP_CRC_overlap.csv",
        ROOT / "outputs" / "DINP_target_master.csv",
        ROOT / "outputs" / "DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv",
        ROOT / "outputs" / "DINP_CRC_STRING_cytohubba_8_consensus_hubs.csv",
    ]
    manifest = {
        "analysis": "DINP_CRC_docking_target_selection",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "primary_target": "PPARG",
        "support_targets": ["PPARA", "RXRA", "NR1I2", "ESR1"],
        "excluded_from_docking": ["CEBPB"],
        "docking_status": "completed_vina_screen",
        "md_status": "deferred_to_external_machine_by_design",
        "source_files": [
            {"path": str(path.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(path)}
            for path in source_files
        ],
        "public_structure_records": {
            "PPARG": ["4XLD", "2PRG"],
            "PPARA": ["6KB1", "3VI8"],
            "RXRA": ["1FBY", "1FM6"],
            "NR1I2": ["8CH8"],
            "ESR1": ["1XPC"],
        },
        "limitations": [
            "DINP is a branched-isomer mixture rather than a single stereochemically defined ligand.",
            "Current target-master DINP binding records are not equivalent to purified-protein binding assays.",
            "Structure availability is not a docking or MD result.",
        ],
    }
    (OUT / "DINP_CRC_docking_target_selection_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"Wrote {len(rows)} shortlisted targets to {OUT}")


if __name__ == "__main__":
    main()
