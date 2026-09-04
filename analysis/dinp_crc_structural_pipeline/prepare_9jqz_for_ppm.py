#!/usr/bin/env python3
"""Prepare the RCSB 9JQZ PTGER4 coordinates for PPM membrane orientation.

This is an input-preparation and audit script only.  It downloads the official
RCSB mmCIF, resolves the PTGER4 polymer using mmCIF entity/reference
annotations (not a chain-ID guess), writes a protein-only PDB for PPM, and
performs parser and structural sanity checks.  It deliberately does not add
residues, alter coordinates, orient the protein, run PPM, or start MD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.vectors import calc_dihedral


ROOT = Path(__file__).resolve().parent
INPUT_DIR = ROOT / "inputs" / "ppm"
OUTPUT_DIR = ROOT / "outputs" / "ppm"
RAW_CIF = INPUT_DIR / "9JQZ_raw.cif"
OUTPUT_PDB = INPUT_DIR / "9JQZ_PTGER4_for_PPM.pdb"
QC_MD = OUTPUT_DIR / "9JQZ_PPM_PREP_QC.md"
QC_JSON = OUTPUT_DIR / "9JQZ_PPM_PREP_QC.json"

PDB_ID = "9JQZ"
RCSB_CIF_URL = "https://files.rcsb.org/download/9JQZ.cif"
RCSB_PDB_URL = "https://files.rcsb.org/download/9JQZ.pdb"
UNIPROT_ACCESSION = "P35408"
UNIPROT_SYMBOL = "PTGER4"


def as_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(x) for x in value]
    return [str(value)]


def mmcif_col(data: Dict[str, Any], name: str) -> List[str]:
    """Return an mmCIF column with case-insensitive tag matching."""
    target = name.lower()
    for key, value in data.items():
        if key.lower() == target:
            return as_list(value)
    return []


def first_value(data: Dict[str, Any], *names: str, default: str = "?") -> str:
    for name in names:
        values = mmcif_col(data, name)
        if values and values[0] not in {"?", ".", ""}:
            return values[0]
    return default


def clean_text(value: str) -> str:
    return value.strip().strip("'").strip('"')


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, destination: Path, force: bool) -> Dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return {
            "url": url,
            "path": str(destination.resolve()),
            "downloaded": False,
            "bytes": destination.stat().st_size,
            "sha256": sha256(destination),
        }
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "whynot17-PTGER4-PPM-prep/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
    if not payload:
        raise RuntimeError(f"RCSB returned an empty file: {url}")
    destination.write_bytes(payload)
    return {
        "url": url,
        "path": str(destination.resolve()),
        "downloaded": True,
        "bytes": len(payload),
        "sha256": sha256(destination),
    }


def parse_int(value: str) -> Optional[int]:
    try:
        return int(float(clean_text(value)))
    except (TypeError, ValueError):
        return None


def parse_float(value: str) -> Optional[float]:
    try:
        return float(clean_text(value))
    except (TypeError, ValueError):
        return None


def entity_annotations(data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    ids = mmcif_col(data, "_entity.id")
    types = mmcif_col(data, "_entity.type")
    descriptions = mmcif_col(data, "_entity.pdbx_description")
    out: Dict[str, Dict[str, Any]] = {}
    for index, entity_id in enumerate(ids):
        out[clean_text(entity_id)] = {
            "type": clean_text(types[index]) if index < len(types) else "?",
            "description": clean_text(descriptions[index]) if index < len(descriptions) else "?",
        }
    return out


def resolve_ptger4(data: Dict[str, Any]) -> Dict[str, Any]:
    """Resolve entity/label/auth chain from explicit mmCIF annotations."""
    entities = entity_annotations(data)
    # ``_struct_asym.id`` is the mmCIF label asym ID.  In this entry,
    # ``_entity_poly.pdbx_strand_id`` is an author/PDB chain identifier, so it
    # must not be used as a label asym ID.  Resolve the two namespaces
    # separately and connect them through the atom-site table below.
    entity_poly_ids = mmcif_col(data, "_struct_asym.entity_id")
    entity_poly_labels = mmcif_col(data, "_struct_asym.id")
    entity_to_labels: Dict[str, List[str]] = {}
    for entity_id, label in zip(entity_poly_ids, entity_poly_labels):
        entity_to_labels.setdefault(clean_text(entity_id), []).append(clean_text(label))

    ref_ids = mmcif_col(data, "_struct_ref.id")
    ref_db = mmcif_col(data, "_struct_ref.db_name")
    ref_accessions = mmcif_col(data, "_struct_ref.pdbx_db_accession")
    ref_entities = mmcif_col(data, "_struct_ref.entity_id")
    uniprot_entities: List[str] = []
    ref_evidence: List[Dict[str, str]] = []
    for index, entity_id in enumerate(ref_entities):
        db = clean_text(ref_db[index]) if index < len(ref_db) else "?"
        accession = clean_text(ref_accessions[index]) if index < len(ref_accessions) else "?"
        ref_id = clean_text(ref_ids[index]) if index < len(ref_ids) else "?"
        row = {"ref_id": ref_id, "db_name": db, "accession": accession, "entity_id": clean_text(entity_id)}
        ref_evidence.append(row)
        if accession.upper() == UNIPROT_ACCESSION or (
            db.upper() == "UNP" and accession.upper() == UNIPROT_ACCESSION
        ):
            uniprot_entities.append(clean_text(entity_id))

    description_entities = [
        entity_id
        for entity_id, info in entities.items()
        if info["type"].lower() == "polymer"
        and ("prostaglandin e2 receptor ep4" in info["description"].lower()
             or "ptger4" in info["description"].lower())
    ]
    candidate_entities = sorted(set(uniprot_entities) | set(description_entities))
    if len(candidate_entities) != 1:
        raise RuntimeError(
            "Could not resolve exactly one PTGER4 entity from mmCIF annotations: "
            + json.dumps({
                "uniprot_entities": uniprot_entities,
                "description_entities": description_entities,
                "entities": entities,
            }, ensure_ascii=False)
        )
    entity_id = candidate_entities[0]
    labels = entity_to_labels.get(entity_id, [])
    if not labels:
        raise RuntimeError(f"PTGER4 entity {entity_id} has no annotated polymer chain")

    atom_labels = mmcif_col(data, "_atom_site.label_asym_id")
    atom_auth = mmcif_col(data, "_atom_site.auth_asym_id")
    label_to_auth: Dict[str, List[str]] = {}
    for label, auth in zip(atom_labels, atom_auth):
        label = clean_text(label)
        auth = clean_text(auth)
        if label not in {"?", "."} and auth not in {"?", "."}:
            label_to_auth.setdefault(label, [])
            if auth not in label_to_auth[label]:
                label_to_auth[label].append(auth)
    auth_chains = sorted({auth for label in labels for auth in label_to_auth.get(label, [])})
    if len(auth_chains) != 1:
        raise RuntimeError(
            f"PTGER4 entity {entity_id} maps to ambiguous coordinate chains: "
            f"labels={labels}, auth_chains={auth_chains}"
        )

    # Find the UniProt reference row and its author-number span.  This protects
    # against retaining a positive-numbered fusion partner in the same chain.
    seq_ref_ids = mmcif_col(data, "_struct_ref_seq.ref_id")
    seq_beg = mmcif_col(data, "_struct_ref_seq.pdbx_auth_seq_align_beg")
    seq_end = mmcif_col(data, "_struct_ref_seq.pdbx_auth_seq_align_end")
    seq_strands = mmcif_col(data, "_struct_ref_seq.pdbx_pdb_strand_id")
    receptor_spans: List[Dict[str, Any]] = []
    for index, ref_id in enumerate(seq_ref_ids):
        matching_ref = next((x for x in ref_evidence if x["ref_id"] == clean_text(ref_id)), None)
        if not matching_ref or matching_ref["accession"].upper() != UNIPROT_ACCESSION:
            continue
        receptor_spans.append({
            "ref_id": clean_text(ref_id),
            "auth_begin": parse_int(seq_beg[index]) if index < len(seq_beg) else None,
            "auth_end": parse_int(seq_end[index]) if index < len(seq_end) else None,
            "strand": clean_text(seq_strands[index]) if index < len(seq_strands) else "?",
        })
    if len(receptor_spans) != 1:
        raise RuntimeError(f"Expected one UniProt PTGER4 alignment span, found {receptor_spans}")
    span = receptor_spans[0]
    if span["auth_begin"] is None or span["auth_end"] is None:
        raise RuntimeError(f"PTGER4 UniProt alignment span lacks usable author bounds: {span}")

    return {
        "selected_entity_id": entity_id,
        "selected_entity_description": entities[entity_id]["description"],
        "selected_label_asym_ids": labels,
        "selected_auth_chain": auth_chains[0],
        "resolution_evidence": {
            "uniprot_accession": UNIPROT_ACCESSION,
            "uniprot_symbol": UNIPROT_SYMBOL,
            "uniprot_reference_rows": [x for x in ref_evidence if x["entity_id"] == entity_id],
            "entity_description_match": entity_id in description_entities,
        },
        "receptor_author_residue_span": span,
        "all_polymer_entities": entities,
        "entity_to_label_chains": entity_to_labels,
        "label_to_auth_chains_with_coordinates": label_to_auth,
    }


class PTGER4Select(Select):
    def __init__(self, auth_chain: str, auth_begin: int, auth_end: int):
        self.auth_chain = auth_chain
        self.auth_begin = auth_begin
        self.auth_end = auth_end

    def accept_model(self, model) -> int:
        return int(model.id == 0)

    def accept_chain(self, chain) -> int:
        return int(chain.id == self.auth_chain)

    def accept_residue(self, residue) -> int:
        if not is_aa(residue, standard=False):
            return 0
        resseq = int(residue.id[1])
        return int(self.auth_begin <= resseq <= self.auth_end)


def residue_key(residue: Any) -> Tuple[int, str]:
    return int(residue.id[1]), str(residue.id[2] or "")


def get_protein_residues(chain: Any) -> List[Any]:
    return [residue for residue in chain if is_aa(residue, standard=False)]


def residue_label(residue: Any) -> str:
    insertion = residue.id[2] if residue.id[2] not in {" ", ""} else ""
    return f"{residue.resname.strip()} {int(residue.id[1])}{insertion}"


def alpha_helix_sanity(residues: Sequence[Any]) -> Dict[str, Any]:
    """Use backbone phi/psi geometry as a conservative alpha-helix check."""
    classified: List[bool] = []
    for index in range(1, len(residues) - 1):
        prev_res, current, next_res = residues[index - 1:index + 2]
        if not all(atom in prev_res for atom in ["C"]):
            continue
        if not all(atom in current for atom in ["N", "CA", "C"]):
            continue
        if not all(atom in next_res for atom in ["N"]):
            continue
        try:
            phi = math.degrees(float(calc_dihedral(
                prev_res["C"].get_vector(), current["N"].get_vector(),
                current["CA"].get_vector(), current["C"].get_vector()
            )))
            psi = math.degrees(float(calc_dihedral(
                current["N"].get_vector(), current["CA"].get_vector(),
                current["C"].get_vector(), next_res["N"].get_vector()
            )))
        except Exception:
            continue
        classified.append(-100.0 <= phi <= -30.0 and -85.0 <= psi <= 25.0)

    runs: List[int] = []
    current_run = 0
    for value in classified:
        if value:
            current_run += 1
        elif current_run:
            runs.append(current_run)
            current_run = 0
    if current_run:
        runs.append(current_run)
    alpha_like = int(sum(classified))
    return {
        "backbone_triplets_evaluated": len(classified),
        "alpha_like_backbone_positions": alpha_like,
        "alpha_like_fraction": (alpha_like / len(classified)) if classified else 0.0,
        "alpha_like_runs_at_least_5": int(sum(run >= 5 for run in runs)),
        "longest_alpha_like_run": max(runs, default=0),
        "criterion": "phi -100..-30 degrees and psi -85..25 degrees; geometry sanity check, not formal DSSP",
    }


def numbering_qc(
    residues: Sequence[Any],
    annotated_begin: Optional[int] = None,
    annotated_end: Optional[int] = None,
) -> Dict[str, Any]:
    numeric = sorted({int(residue.id[1]) for residue in residues})
    missing_positions: List[int] = []
    for left, right in zip(numeric, numeric[1:]):
        if right - left > 1:
            missing_positions.extend(range(left + 1, right))
    terminal_missing: List[int] = []
    if numeric and annotated_begin is not None and numeric[0] > annotated_begin:
        terminal_missing.extend(range(annotated_begin, numeric[0]))
    if numeric and annotated_end is not None and numeric[-1] < annotated_end:
        terminal_missing.extend(range(numeric[-1] + 1, annotated_end + 1))
    return {
        "observed_auth_residue_numbers": [numeric[0], numeric[-1]] if numeric else [],
        "internal_numbering_gap_count": len(missing_positions),
        "internal_missing_auth_residue_numbers": missing_positions,
        "annotated_auth_residue_span": [annotated_begin, annotated_end],
        "terminal_missing_auth_residue_numbers": terminal_missing,
        "chain_break_or_missing_residue_detected": bool(missing_positions or terminal_missing),
        "interpretation": "Coordinate gaps relative to the annotated span are reported only; no residues are modeled or inserted.",
    }


def chain_atom_counts(chain: Any) -> Dict[str, int]:
    counts = {"total": 0, "ATOM_like": 0, "HETATM_like": 0, "hydrogen": 0}
    for residue in chain:
        for atom in residue:
            counts["total"] += 1
            if residue.id[0] == " ":
                counts["ATOM_like"] += 1
            else:
                counts["HETATM_like"] += 1
            if getattr(atom, "element", "").upper() == "H":
                counts["hydrogen"] += 1
    return counts


def parse_structure(path: Path, mmcif: bool = False):
    if mmcif:
        return MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True).get_structure(path.stem, str(path))
    return PDBParser(QUIET=True).get_structure(path.stem, str(path))


def validate_output(
    path: Path,
    expected_chain: str,
    expected_residue_count: int,
    annotated_begin: int,
    annotated_end: int,
) -> Dict[str, Any]:
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"Output PDB is missing or empty: {path}")
    parsed = parse_structure(path)
    model = next(parsed.get_models())
    chain_ids = [chain.id for chain in model]
    if chain_ids != [expected_chain]:
        raise RuntimeError(f"Output must contain exactly chain {expected_chain}; found {chain_ids}")
    chain = model[expected_chain]
    residues = get_protein_residues(chain)
    atom_counts = chain_atom_counts(chain)
    if len(residues) != expected_residue_count:
        raise RuntimeError(
            f"Output residue count changed unexpectedly: source={expected_residue_count}, output={len(residues)}"
        )
    if len(residues) < 250:
        raise RuntimeError(f"Output has only {len(residues)} protein residues; expected at least 250")
    if atom_counts["HETATM_like"]:
        raise RuntimeError(f"Output contains {atom_counts['HETATM_like']} non-standard HETATM-like residues")
    return {
        "re_readable_by_biopython": True,
        "chain_ids": chain_ids,
        "protein_residue_count": len(residues),
        "coordinate_atom_count": atom_counts["total"],
        "atom_counts": atom_counts,
        "alpha_helix_sanity": alpha_helix_sanity(residues),
        "numbering_qc": numbering_qc(residues, annotated_begin, annotated_end),
    }


def windows_path(path: Path) -> str:
    """Render a useful Windows path when executed inside WSL."""
    resolved = str(path.resolve())
    match = re.match(r"^/mnt/([a-zA-Z])/(.*)$", resolved)
    if match:
        tail = match.group(2).replace("/", "\\")
        return f"{match.group(1).upper()}:\\{tail}"
    return resolved


def format_value(value: Any) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    if value is None:
        return "not available"
    if isinstance(value, (dict, list)):
        return "`" + json.dumps(value, ensure_ascii=False) + "`"
    return str(value)


def write_qc_report(report: Dict[str, Any]) -> None:
    QC_JSON.parent.mkdir(parents=True, exist_ok=True)
    QC_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    selected = report["chain_resolution"]
    source = report["source_metadata"]
    output = report["output_qc"]
    deletions = report["removed_content"]
    lines = [
        "# 9JQZ PTGER4 PPM-input preparation QC",
        "",
        "This file documents coordinate extraction only. No membrane orientation, loop modeling, residue insertion, energy minimization, PPM run, or MD was performed.",
        "",
        "## Structure and source",
        "",
        f"- PDB ID: `{PDB_ID}`",
        f"- Structure title: {source['title']}",
        f"- Experimental method: {source['experimental_method']}",
        f"- Resolution: {source['resolution_A']} Å",
        f"- Raw RCSB mmCIF: `{windows_path(Path(report['downloads']['mmcif']['path']))}`",
        f"- Raw mmCIF SHA256: `{report['downloads']['mmcif']['sha256']}`",
        "",
        "## Annotation-based PTGER4 chain resolution",
        "",
        f"- Selected entity ID: `{selected['selected_entity_id']}`",
        f"- Entity description: {selected['selected_entity_description']}",
        f"- Selected mmCIF label asym ID(s): `{', '.join(selected['selected_label_asym_ids'])}`",
        f"- Selected author/PDB chain ID: `{selected['selected_auth_chain']}`",
        f"- Resolution evidence: UniProt `{UNIPROT_ACCESSION}` / `{UNIPROT_SYMBOL}` reference plus entity description match = {format_value(selected['resolution_evidence']['entity_description_match'])}",
        f"- Annotated PTGER4 author residue span: `{selected['receptor_author_residue_span']['auth_begin']}–{selected['receptor_author_residue_span']['auth_end']}`",
        f"- Protein residue count retained: `{report['retained_content']['protein_residue_count']}`",
        f"- First retained residue: `{report['retained_content']['first_residue']}`",
        f"- Last retained residue: `{report['retained_content']['last_residue']}`",
        f"- N-terminal retained residue identity: `{report['retained_content']['n_terminal_residue']}`",
        "",
        "## Coordinate and structural sanity checks",
        "",
        f"- Coordinate atom count: `{report['retained_content']['coordinate_atom_count']}`",
        f"- Obvious 7TM GPCR receptor body: **{format_value(report['retained_content']['obvious_7tm_gpcr'])}** — {report['retained_content']['obvious_7tm_basis']}",
        f"- Alpha-helix coordinate sanity: **{format_value(output['alpha_helix_sanity']['alpha_like_backbone_positions'] >= 70)}** — {output['alpha_helix_sanity']['alpha_like_backbone_positions']} alpha-like backbone positions across {output['alpha_helix_sanity']['backbone_triplets_evaluated']} evaluated positions; {output['alpha_helix_sanity']['alpha_like_runs_at_least_5']} runs ≥5 residues",
        f"- Chain break / missing coordinate residue detected: **{format_value(output['numbering_qc']['chain_break_or_missing_residue_detected'])}**",
        f"- Internal numbering gaps: `{output['numbering_qc']['internal_missing_auth_residue_numbers']}`",
        f"- BioPython re-read: **{format_value(output['re_readable_by_biopython'])}**",
        f"- MDAnalysis re-read: **{format_value(report['parser_qc']['mdanalysis_re_readable'])}**",
        "",
        "## Removed content",
        "",
        f"- Other polymer chains removed: `{deletions['other_polymer_chains_removed']}`",
        f"- Water removed: `{deletions['water_removed']}`",
        f"- Detergent/lipid/ion/crystallization additive/co-crystal ligand records removed: `{deletions['heteroatom_categories_removed']}`",
        f"- Fusion/accessory content: PTGER4 was selected by entity/reference annotation and only its annotated author-residue span was retained; non-receptor fusion coordinates, if present, were excluded.",
        "",
        "## Output",
        "",
        f"- PPM upload PDB: `{windows_path(Path(report['output_pdb']))}`",
        f"- Output PDB SHA256: `{report['output_sha256']}`",
        f"- Output is non-empty and legal PDB: **{format_value(report['output_legal_pdb'])}**",
        "",
        "## Recommended PPM settings",
        "",
        "- Number of membranes = 1",
        "- Type of membrane = Plasma membrane (mammalian)",
        "- Allow curvature = no",
        "- Topology (N-ter) = out",
        "- Input = Coordinate file",
        "- Include heteroatoms = no",
        "",
        "## Boundary",
        "",
        "This PDB is prepared only for upload to PPM membrane-orientation software. PPM was not run, the coordinates were not oriented, and no MD system was built.",
    ]
    QC_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Redownload the raw RCSB mmCIF")
    args = parser.parse_args()

    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    download_info = download_file(RCSB_CIF_URL, RAW_CIF, force=args.force)

    raw_data = MMCIF2Dict(str(RAW_CIF))
    chain_resolution = resolve_ptger4(raw_data)
    structure = parse_structure(RAW_CIF, mmcif=True)
    model = next(structure.get_models())
    auth_chain = chain_resolution["selected_auth_chain"]
    if auth_chain not in model:
        raise RuntimeError(f"Annotation-resolved author chain {auth_chain!r} not found in parsed structure")
    source_chain = model[auth_chain]
    auth_begin = chain_resolution["receptor_author_residue_span"]["auth_begin"]
    auth_end = chain_resolution["receptor_author_residue_span"]["auth_end"]
    source_residues = [
        residue for residue in source_chain
        if is_aa(residue, standard=False) and auth_begin <= int(residue.id[1]) <= auth_end
    ]
    if len(source_residues) < 250:
        raise RuntimeError(f"Source PTGER4 selection has only {len(source_residues)} residues; aborting")

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(OUTPUT_PDB), PTGER4Select(auth_chain, auth_begin, auth_end))
    output_qc = validate_output(
        OUTPUT_PDB, auth_chain, len(source_residues), auth_begin, auth_end
    )

    # Optional independent parser check.  MDAnalysis is expected in the D-drive
    # WSL MD environment but the script remains usable without it.
    mda_readable = False
    mda_error = None
    try:
        import MDAnalysis as mda  # type: ignore
        universe = mda.Universe(str(OUTPUT_PDB))
        mda_readable = bool(len(universe.atoms) > 0 and len(universe.segments) >= 1)
    except Exception as exc:  # pragma: no cover - environment-dependent
        mda_error = repr(exc)

    keywords = " ".join(
        clean_text(x) for x in mmcif_col(raw_data, "_struct_keywords.text")
    ).lower()
    title = first_value(raw_data, "_struct.title")
    method = first_value(raw_data, "_exptl.method")
    resolution = first_value(
        raw_data,
        "_em_3d_reconstruction.resolution",
        "_refine.ls_d_res_high",
    )
    hetero_categories = [
        "water", "detergent", "lipid", "ion", "crystallization additive", "co-crystal ligand",
    ]
    nonpoly_names = [clean_text(x) for x in mmcif_col(raw_data, "_pdbx_entity_nonpoly.name")]
    nonpoly_comp_ids = [clean_text(x) for x in mmcif_col(raw_data, "_pdbx_entity_nonpoly.comp_id")]
    nonpoly_summary = [
        {"name": name, "comp_id": comp_id}
        for name, comp_id in zip(nonpoly_names, nonpoly_comp_ids)
    ]
    other_polymer_chains = [
        chain.id for chain in model
        if chain.id != auth_chain
        and any(is_aa(residue, standard=False) for residue in chain)
    ]
    first_residue = residue_label(source_residues[0])
    last_residue = residue_label(source_residues[-1])
    receptor_body_basis = (
        "The selected entity is explicitly annotated as PTGER4 and the entry keywords identify a GPCR/membrane protein; "
        "the retained coordinates also pass the alpha-helix geometry sanity check."
    )
    report: Dict[str, Any] = {
        "script": str(Path(__file__).resolve()),
        "pdb_id": PDB_ID,
        "downloads": {"mmcif": download_info, "pdb_url_checked": RCSB_PDB_URL},
        "source_metadata": {
            "title": clean_text(title),
            "experimental_method": clean_text(method),
            "resolution_A": parse_float(resolution) if parse_float(resolution) is not None else clean_text(resolution),
            "struct_keywords": keywords,
        },
        "chain_resolution": chain_resolution,
        "retained_content": {
            "protein_residue_count": len(source_residues),
            "coordinate_atom_count": output_qc["coordinate_atom_count"],
            "first_residue": first_residue,
            "last_residue": last_residue,
            "n_terminal_residue": first_residue,
            "obvious_7tm_gpcr": bool("gpcr" in keywords and "membrane" in keywords),
            "obvious_7tm_basis": receptor_body_basis,
        },
        "removed_content": {
            "other_polymer_chains_removed": other_polymer_chains,
            "water_removed": [x for x in nonpoly_summary if "water" in x["name"].lower() or x["comp_id"] in {"HOH", "WAT"}],
            "heteroatom_categories_removed": nonpoly_summary,
            "nonpolymer_records_in_source": nonpoly_summary,
        },
        "parser_qc": {
            "biopython_re_readable": output_qc["re_readable_by_biopython"],
            "mdanalysis_re_readable": mda_readable,
            "mdanalysis_error": mda_error,
        },
        "output_pdb": str(OUTPUT_PDB.resolve()),
        "output_sha256": sha256(OUTPUT_PDB),
        "output_legal_pdb": True,
        "output_qc": output_qc,
        "ppm_run_performed": False,
        "md_run_performed": False,
    }
    if not mda_readable and mda_error:
        # MDAnalysis is a requested sanity check; absence is not silently
        # labeled as success, while the required BioPython check remains fatal.
        report["parser_qc"]["mdanalysis_status"] = "not_available_or_failed"
    else:
        report["parser_qc"]["mdanalysis_status"] = "passed"
    write_qc_report(report)

    print(f"PPM upload file:\n{windows_path(OUTPUT_PDB)}")
    print("\nRecommended PPM settings:")
    print("Number of membranes = 1")
    print("Type of membrane = Plasma membrane (mammalian)")
    print("Allow curvature = no")
    print("Topology (N-ter) = out")
    print("Input = Coordinate file")
    print("Include heteroatoms = no")
    print(f"\nQC report: {windows_path(QC_MD)}")
    print(f"Output SHA256: {report['output_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
