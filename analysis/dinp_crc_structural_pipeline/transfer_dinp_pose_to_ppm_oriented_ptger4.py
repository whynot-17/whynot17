#!/usr/bin/env python3
"""Transfer the validated DINP--PTGER4 pose into the local PPM frame.

This is deliberately a coordinate-transfer-only workflow.  It does not redock,
minimize, regenerate, or otherwise alter the ligand.  The receptor-derived
Kabsch transform is the sole transform applied to the DINP pose.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa


ROOT = Path(__file__).resolve().parents[2]
PIPELINE = ROOT / "analysis" / "dinp_crc_structural_pipeline"
TARGET_RECEPTOR = PIPELINE / "outputs" / "ppm_local" / "9JQZ_PTGER4_PPM_oriented.pdb"
OUTPUT_DIR = PIPELINE / "outputs" / "ppm_pose_transfer"
COMPLEX_PDB = OUTPUT_DIR / "PTGER4_DINP_PPM_oriented_complex.pdb"
NATIVE_DUM_COMPLEX_PDB = OUTPUT_DIR / "PTGER4_DINP_PPM_oriented_complex_native_dum.pdb"
LIGAND_PDB = OUTPUT_DIR / "DINP_PPM_oriented.pdb"
TRANSFORM_JSON = OUTPUT_DIR / "original_to_ppm_transform.json"
QC_MD = OUTPUT_DIR / "PTGER4_DINP_PPM_POSE_TRANSFER_QC.md"

EXPECTED_AFFINITY = -7.429
RECEPTOR_RMSD_LIMIT = 0.01
LIGAND_INTERNAL_RMSD_LIMIT = 1e-5
POCKET_RMSD_LIMIT = 1e-4


def fail(message: str) -> "NoReturn":
    raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def recorded_path(value: str | Path) -> Path:
    """Resolve Windows paths recorded by the earlier Windows docking run."""
    raw = str(value)
    direct = Path(raw)
    if direct.exists():
        return direct.resolve()
    normalized = raw.replace("\\", "/")
    match = re.match(r"^([A-Za-z]):/(.*)$", normalized)
    if match:
        candidate = Path("/mnt") / match.group(1).lower() / match.group(2)
        if candidate.exists():
            return candidate.resolve()
    # Keep the search deterministic if provenance contains an old worktree
    # path but the same repository is currently mounted elsewhere.
    marker = "/analysis/dinp_crc_structural_pipeline/"
    if marker in normalized:
        relative = normalized.split(marker, 1)[1]
        candidate = PIPELINE / relative
        if candidate.exists():
            return candidate.resolve()
    fail(f"provenance path does not exist: {value}")


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        fail(f"could not read JSON provenance {path}: {exc}")


def find_final_ptger4_audit() -> tuple[Path, dict, dict]:
    candidates = [PIPELINE / "outputs" / "stage2_refine_and_rescue" / "audit.json"]
    candidates.extend(sorted((PIPELINE / "outputs").rglob("audit.json")))
    seen: set[Path] = set()
    for audit_path in candidates:
        audit_path = audit_path.resolve()
        if audit_path in seen or not audit_path.exists():
            continue
        seen.add(audit_path)
        data = load_json(audit_path)
        info = data.get("ptger4")
        if not isinstance(info, dict):
            continue
        try:
            affinity = float(info.get("affinity_kcal_mol"))
            pose = recorded_path(info["pose"])
            receptor = recorded_path(info["clean_pdb"])
        except (KeyError, TypeError, ValueError, RuntimeError):
            continue
        if abs(affinity - EXPECTED_AFFINITY) > 1e-6:
            continue
        if "rescue_pose" not in pose.name or not pose.exists() or not receptor.exists():
            continue
        return audit_path, data, info
    fail("could not locate the final PTGER4 rescue audit with DINP affinity -7.429 kcal/mol")


def parse_protocol(log_path: Path) -> dict:
    text = log_path.read_text(errors="replace")
    center_match = re.search(
        r"Grid center:\s*X\s*([-+]?\d+(?:\.\d+)?)\s*Y\s*([-+]?\d+(?:\.\d+)?)\s*Z\s*([-+]?\d+(?:\.\d+)?)",
        text,
    )
    size_match = re.search(
        r"Grid size\s*:\s*X\s*([-+]?\d+(?:\.\d+)?)\s*Y\s*([-+]?\d+(?:\.\d+)?)\s*Z\s*([-+]?\d+(?:\.\d+)?)",
        text,
    )
    ex_match = re.search(r"Exhaustiveness:\s*(\d+)", text)
    seed_match = re.search(r"random seed:\s*(\d+)", text)
    if not (center_match and size_match and ex_match and seed_match):
        fail(f"incomplete Vina protocol provenance in {log_path}")
    return {
        "center": [float(value) for value in center_match.groups()],
        "size": [float(value) for value in size_match.groups()],
        "exhaustiveness": int(ex_match.group(1)),
        "seed": int(seed_match.group(1)),
    }


def parse_vina_affinity(path: Path) -> float:
    match = re.search(r"REMARK VINA RESULT:\s*([-+]?\d+(?:\.\d+)?)", path.read_text(errors="replace"))
    if not match:
        fail(f"no Vina affinity found in selected pose: {path}")
    return float(match.group(1))


def atom_element(atom) -> str:
    element = str(getattr(atom, "element", "") or "").strip().upper()
    if element:
        return element
    name = atom.get_name().strip().upper()
    return re.sub(r"^[0-9]+", "", name)[:2]


def protein_atom_rows(path: Path) -> list[dict]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    rows: list[dict] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if not is_aa(residue, standard=False):
                    continue
                hetflag, resseq, icode = residue.id
                for atom in residue:
                    key = (
                        str(chain.id),
                        int(resseq),
                        str(icode).strip(),
                        residue.get_resname().strip(),
                        atom.get_name().strip(),
                    )
                    rows.append(
                        {
                            "key": key,
                            "residue_key": key[:4],
                            "chain": str(chain.id),
                            "resseq": int(resseq),
                            "icode": str(icode).strip(),
                            "resname": residue.get_resname().strip(),
                            "atom_name": atom.get_name().strip(),
                            "element": atom_element(atom),
                            "coord": np.asarray(atom.coord, dtype=float),
                        }
                    )
        break
    return rows


def heavy(row: dict) -> bool:
    return not row["element"].startswith("H") and not row["atom_name"].upper().startswith("H")


def row_map(rows: list[dict]) -> dict[tuple, dict]:
    result: dict[tuple, dict] = {}
    for row in rows:
        if row["key"] in result:
            fail(f"duplicate receptor atom key encountered: {row['key']}")
        result[row["key"]] = row
    return result


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    pre = float(np.sqrt(np.mean(np.sum((source - target) ** 2, axis=1))))
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    source_c = source - source_center
    target_c = target - target_center
    covariance = source_c.T @ target_c
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    translation = target_center - source_center @ rotation
    fitted = source @ rotation + translation
    post = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rotation, translation, pre, post


def pairwise_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) != len(b):
        fail("pairwise RMSD received arrays with different lengths")
    if len(a) < 2:
        return 0.0
    da = np.sqrt(np.maximum(0.0, ((a[:, None, :] - a[None, :, :]) ** 2).sum(axis=2)))
    db = np.sqrt(np.maximum(0.0, ((b[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)))
    tri = np.triu_indices(len(a), k=1)
    return float(np.sqrt(np.mean((da[tri] - db[tri]) ** 2)))


def distance_matrix_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        fail(f"distance matrix shape mismatch: {a.shape} vs {b.shape}")
    da = np.sqrt(np.maximum(0.0, ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)))
    db = np.sqrt(np.maximum(0.0, ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)))
    return float(np.sqrt(np.mean((da - db) ** 2)))


def parse_pdbqt_pose(path: Path) -> list[dict]:
    atoms: list[dict] = []
    in_model = False
    saw_model = False
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("MODEL"):
            if saw_model:
                break
            saw_model = True
            in_model = True
            continue
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if saw_model and not in_model:
            continue
        if len(line) < 54:
            fail(f"malformed PDBQT atom record in {path}: {line!r}")
        try:
            x, y, z = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError as exc:
            fail(f"could not parse PDBQT coordinates in {path}: {line!r} ({exc})")
        atoms.append(
            {
                "serial": int(line[6:11].strip() or len(atoms) + 1),
                "atom_name": line[12:16].strip(),
                "resname": line[17:20].strip() or "UNL",
                "chain": line[21].strip(),
                "resseq": int(line[22:26].strip() or 1),
                "coord": np.array([x, y, z], dtype=float),
                "element": (line[77:79].strip() if len(line) >= 79 else ""),
                "source_line": line,
            }
        )
    if not atoms:
        fail(f"selected DINP pose contains no atom records: {path}")
    for atom in atoms:
        atom_type = atom["element"].upper()
        # PDBQT uses A for aromatic carbon and OA/NA for typed heteroatoms;
        # change only the output element field while retaining atom names.
        if atom_type == "A":
            atom["element"] = "C"
        elif atom_type in {"OA", "OS"}:
            atom["element"] = "O"
        elif atom_type in {"NA", "N"}:
            atom["element"] = "N"
        elif atom_type in {"HD", "H"}:
            atom["element"] = "H"
        elif not atom_type:
            atom["element"] = re.sub(r"^[0-9]+", "", atom["atom_name"]).upper()[:2]
    return atoms


def ligand_heavy_atoms(atoms: list[dict]) -> list[dict]:
    return [atom for atom in atoms if not atom["element"].upper().startswith("H") and not atom["atom_name"].upper().startswith("H")]


def ligand_pdb_lines(
    atoms: list[dict], rotation: np.ndarray, translation: np.ndarray
) -> tuple[list[str], list[dict]]:
    lines = []
    name_mapping = []
    for serial, atom in enumerate(atoms, start=9001):
        x, y, z = atom["coord"] @ rotation + translation
        output_name = f"{atom['element'][:1].upper()}{serial - 9000:02d}"[:4]
        element = atom["element"].upper()[:2].rjust(2)
        name_mapping.append(
            {
                "output_pdb_atom_name": output_name,
                "original_pdbqt_atom_name": atom["atom_name"],
                "serial": serial,
            }
        )
        lines.append(
            f"HETATM{serial:5d} {output_name:<4s} DIN L{1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element}\n"
        )
    return lines, name_mapping


def strip_end_records(lines: list[str]) -> list[str]:
    return [line for line in lines if not line.startswith(("END", "ENDMDL"))]


def try_biopython(path: Path) -> dict:
    try:
        structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
        return {"readable": True, "atoms": int(len(list(structure.get_atoms()))), "error": ""}
    except Exception as exc:
        return {"readable": False, "atoms": None, "error": str(exc)}


def try_mdanalysis(path: Path) -> dict:
    try:
        import MDAnalysis as mda

        universe = mda.Universe(str(path))
        return {"readable": True, "atoms": int(len(universe.atoms)), "residues": int(len(universe.residues)), "error": ""}
    except Exception as exc:
        return {"readable": False, "atoms": None, "residues": None, "error": str(exc)}


def main() -> int:
    if not TARGET_RECEPTOR.exists():
        fail(f"PPM-oriented receptor not found: {TARGET_RECEPTOR}")

    audit_path, stage2_audit, stage2_info = find_final_ptger4_audit()
    source_receptor = recorded_path(stage2_info["clean_pdb"])
    selected_prepared_receptor = recorded_path(stage2_info["receptor_pdbqt"])
    ligand_source = recorded_path(stage2_info["pose"])
    rescue_log = ligand_source.with_name("DINP_PTGER4_9JQZ_rescue_vina.log")
    control_audit_path = PIPELINE / "outputs" / "ptger4_control_docking" / "audit.json"
    if not control_audit_path.exists():
        fail(f"PTGER4 control docking provenance not found: {control_audit_path}")
    control_audit = load_json(control_audit_path)
    control_receptor = recorded_path(control_audit["receptor"])
    control_log = PIPELINE / "outputs" / "ptger4_control_docking" / "DINP_PTGER4_same_protocol.log"
    if not control_log.exists():
        fail(f"PTGER4 control protocol log not found: {control_log}")

    selected_affinity = parse_vina_affinity(ligand_source)
    if abs(selected_affinity - EXPECTED_AFFINITY) > 1e-6:
        fail(f"selected pose affinity is {selected_affinity}, expected {EXPECTED_AFFINITY}")
    if abs(float(control_audit.get("dinp_affinity_kcal_mol")) - EXPECTED_AFFINITY) > 1e-6:
        fail("PTGER4 control docking provenance does not report the expected DINP affinity")
    if selected_prepared_receptor.resolve() != control_receptor.resolve():
        fail("selected DINP pose and PTGER4 control docking do not use the same prepared receptor")

    selected_protocol = parse_protocol(rescue_log)
    control_protocol = parse_protocol(control_log)
    protocol_delta = {
        "center_max_abs_delta": max(abs(a - b) for a, b in zip(selected_protocol["center"], control_protocol["center"])),
        "size_max_abs_delta": max(abs(a - b) for a, b in zip(selected_protocol["size"], control_protocol["size"])),
        "exhaustiveness_same": selected_protocol["exhaustiveness"] == control_protocol["exhaustiveness"],
        "seed_same": selected_protocol["seed"] == control_protocol["seed"],
    }
    same_protocol = (
        protocol_delta["center_max_abs_delta"] < 1e-4
        and protocol_delta["size_max_abs_delta"] < 1e-4
        and protocol_delta["exhaustiveness_same"]
        and protocol_delta["seed_same"]
    )
    if not same_protocol:
        fail(f"selected pose and control protocol differ: {protocol_delta}")

    source_rows = protein_atom_rows(source_receptor)
    target_rows = protein_atom_rows(TARGET_RECEPTOR)
    source_map = row_map(source_rows)
    target_map = row_map(target_rows)
    common_keys = sorted(set(source_map) & set(target_map), key=str)
    unmatched_source = sorted(set(source_map) - set(target_map), key=str)
    unmatched_target = sorted(set(target_map) - set(source_map), key=str)
    if len(common_keys) < 250:
        fail(f"only {len(common_keys)} common receptor atoms matched")

    source_coords = np.array([source_map[key]["coord"] for key in common_keys])
    target_coords = np.array([target_map[key]["coord"] for key in common_keys])
    rotation, translation, receptor_prefit_rmsd, receptor_postfit_rmsd = kabsch(source_coords, target_coords)

    def metric_for(predicate):
        keys = [key for key in common_keys if predicate(source_map[key]) and predicate(target_map[key])]
        src = np.array([source_map[key]["coord"] for key in keys])
        tgt = np.array([target_map[key]["coord"] for key in keys])
        transformed = src @ rotation + translation
        rmsd = float(np.sqrt(np.mean(np.sum((transformed - tgt) ** 2, axis=1))))
        return len(keys), rmsd

    ca_count, ca_rmsd = metric_for(lambda row: row["atom_name"] == "CA")
    backbone_count, backbone_rmsd = metric_for(lambda row: row["atom_name"] in {"N", "CA", "C", "O"})
    heavy_count, heavy_rmsd = metric_for(heavy)

    ligand_atoms = parse_pdbqt_pose(ligand_source)
    ligand_heavy = ligand_heavy_atoms(ligand_atoms)
    ligand_coords = np.array([atom["coord"] for atom in ligand_heavy])
    transformed_ligand_coords = ligand_coords @ rotation + translation
    ligand_internal_rmsd = pairwise_rmsd(ligand_coords, transformed_ligand_coords)

    source_receptor_heavy = [row for row in source_rows if heavy(row)]
    source_receptor_heavy_coords = np.array([row["coord"] for row in source_receptor_heavy])
    ligand_to_receptor_dist = np.sqrt(
        ((ligand_coords[:, None, :] - source_receptor_heavy_coords[None, :, :]) ** 2).sum(axis=2)
    )
    residue_hits: dict[tuple, float] = {}
    for receptor_index, row in enumerate(source_receptor_heavy):
        residue_hits[row["residue_key"]] = min(
            residue_hits.get(row["residue_key"], math.inf), float(ligand_to_receptor_dist[:, receptor_index].min())
        )
    pocket_residues = {key for key, distance in residue_hits.items() if distance <= 5.0}
    if not pocket_residues:
        fail("no receptor residues within 5.0 Å of the selected DINP pose")
    pocket_rows = [row for row in source_receptor_heavy if row["residue_key"] in pocket_residues]
    pocket_coords = np.array([row["coord"] for row in pocket_rows])
    transformed_pocket_coords = pocket_coords @ rotation + translation
    source_relative_distances = np.sqrt(
        ((ligand_coords[:, None, :] - pocket_coords[None, :, :]) ** 2).sum(axis=2)
    )
    transformed_relative_distances = np.sqrt(
        ((transformed_ligand_coords[:, None, :] - transformed_pocket_coords[None, :, :]) ** 2).sum(axis=2)
    )
    pocket_pairwise_rmsd = float(np.sqrt(np.mean((source_relative_distances - transformed_relative_distances) ** 2)))
    source_min_pocket_distance = float(source_relative_distances.min())
    transformed_min_pocket_distance = float(transformed_relative_distances.min())
    source_ligand_com = ligand_coords.mean(axis=0)
    transformed_ligand_com = transformed_ligand_coords.mean(axis=0)
    source_pocket_com = pocket_coords.mean(axis=0)
    transformed_pocket_com = transformed_pocket_coords.mean(axis=0)
    source_com_pocket_distance = float(np.linalg.norm(source_ligand_com - source_pocket_com))
    transformed_com_pocket_distance = float(np.linalg.norm(transformed_ligand_com - transformed_pocket_com))
    nearest_contacts_before = np.sort(source_relative_distances.ravel())[:10].tolist()
    nearest_contacts_after = np.sort(transformed_relative_distances.ravel())[:10].tolist()

    native_lines = TARGET_RECEPTOR.read_text(errors="replace").splitlines(keepends=True)
    dum_z = []
    for line in native_lines:
        if line.startswith("HETATM") and line[17:20].strip() == "DUM":
            try:
                dum_z.append(float(line[46:54]))
            except ValueError:
                pass
    if not dum_z:
        fail("could not parse PPM DUM membrane-boundary z-coordinates")
    membrane_center_z = (min(dum_z) + max(dum_z)) / 2.0
    membrane_lower_z = min(dum_z)
    membrane_upper_z = max(dum_z)
    ligand_com_z = float(transformed_ligand_com[2])
    membrane_half_thickness = max(abs(membrane_lower_z - membrane_center_z), abs(membrane_upper_z - membrane_center_z))
    ligand_near_pocket = transformed_min_pocket_distance <= 5.0

    ligand_lines, ligand_name_mapping = ligand_pdb_lines(ligand_atoms, rotation, translation)
    ligand_name_remarks = [
        "REMARK 950 DINP PDB atom names are unique serialization names; original PDBQT names are mapped below.\n",
        *[
            f"REMARK 950 MAP {entry['output_pdb_atom_name']:<4s} <- {entry['original_pdbqt_atom_name']:<4s} SERIAL {entry['serial']}\n"
            for entry in ligand_name_mapping
        ],
    ]
    target_protein_lines = [line for line in native_lines if line.startswith(("REMARK", "ATOM  ", "TER"))]
    native_complex_lines = strip_end_records(native_lines) + ligand_name_remarks + ligand_lines + ["END\n"]
    protein_complex_lines = target_protein_lines + ligand_name_remarks + ligand_lines + ["END\n"]
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LIGAND_PDB.write_text("".join(ligand_name_remarks + ligand_lines + ["END\n"]), encoding="utf-8")
    COMPLEX_PDB.write_text("".join(protein_complex_lines), encoding="utf-8")
    NATIVE_DUM_COMPLEX_PDB.write_text("".join(native_complex_lines), encoding="utf-8")

    biopython_complex = try_biopython(COMPLEX_PDB)
    mdanalysis_complex = try_mdanalysis(COMPLEX_PDB)
    biopython_native_complex = try_biopython(NATIVE_DUM_COMPLEX_PDB)
    mdanalysis_native_complex = try_mdanalysis(NATIVE_DUM_COMPLEX_PDB)

    qc_pass = (
        ca_rmsd < RECEPTOR_RMSD_LIMIT
        and backbone_rmsd < RECEPTOR_RMSD_LIMIT
        and heavy_rmsd < RECEPTOR_RMSD_LIMIT
        and ligand_internal_rmsd < LIGAND_INTERNAL_RMSD_LIMIT
        and pocket_pairwise_rmsd < POCKET_RMSD_LIMIT
        and biopython_complex["readable"]
        and mdanalysis_complex["readable"]
        and ligand_near_pocket
    )

    timestamp = datetime.now(timezone.utc).isoformat()
    transform_payload = {
        "source_receptor_path": str(source_receptor.resolve()),
        "target_receptor_path": str(TARGET_RECEPTOR.resolve()),
        "ligand_source_path": str(ligand_source.resolve()),
        "selected_docking_affinity_kcal_mol": selected_affinity,
        "ptger4_control_audit": str(control_audit_path.resolve()),
        "matched_atom_count": len(common_keys),
        "matching_key": "chain ID + residue number + insertion code + residue name + atom name",
        "rotation_matrix": rotation.tolist(),
        "translation_vector_angstrom": translation.tolist(),
        "receptor_prefit_rmsd_angstrom": receptor_prefit_rmsd,
        "receptor_postfit_rmsd_angstrom": receptor_postfit_rmsd,
        "git_commit_if_available": None,
        "timestamp_utc": timestamp,
        "selected_protocol": selected_protocol,
        "control_protocol": control_protocol,
        "same_receptor_pocket_protocol_verified": same_protocol,
        "ligand_atom_name_mapping": ligand_name_mapping,
    }
    try:
        import subprocess

        transform_payload["git_commit_if_available"] = subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        transform_payload["git_commit_if_available"] = "unavailable"
    TRANSFORM_JSON.write_text(json.dumps(transform_payload, indent=2) + "\n", encoding="utf-8")

    report = [
        "# PTGER4–DINP PPM Pose-Transfer QC",
        "",
        f"- Overall pose-transfer QC: `{'PASS' if qc_pass else 'FAIL'}`",
        f"- Selected DINP affinity: `{selected_affinity:.3f} kcal/mol`",
        f"- Timestamp (UTC): `{timestamp}`",
        "",
        "## Provenance and protocol identity",
        "",
        f"- Source receptor: `{source_receptor}`",
        f"- PPM-oriented receptor: `{TARGET_RECEPTOR}`",
        f"- Source DINP pose: `{ligand_source}`",
        f"- Stage-2 audit: `{audit_path}`",
        f"- PTGER4 control audit: `{control_audit_path}`",
        f"- Same prepared docking receptor: `{'PASS' if selected_prepared_receptor.resolve() == control_receptor.resolve() else 'FAIL'}`",
        f"- Same pocket/protocol: `{'PASS' if same_protocol else 'FAIL'}`",
        f"- Selected Vina affinity from pose remark: `{selected_affinity:.3f} kcal/mol`",
        "",
        "## A. Receptor mapping and rigid fit",
        "",
        f"- Mapping key: `chain + residue number + insertion code + residue name + atom name`",
        f"- Common receptor atoms: `{len(common_keys)}`",
        f"- Unmatched source atoms: `{len(unmatched_source)}`",
        f"- Unmatched PPM atoms: `{len(unmatched_target)}`",
        f"- Common Cα atoms / RMSD after fit: `{ca_count}` / `{ca_rmsd:.9g} Å` (threshold `<0.01 Å`)",
        f"- Common backbone atoms / RMSD after fit: `{backbone_count}` / `{backbone_rmsd:.9g} Å` (threshold `<0.01 Å`)",
        f"- Common heavy atoms / RMSD after fit: `{heavy_count}` / `{heavy_rmsd:.9g} Å` (threshold `<0.01 Å`)",
        f"- Receptor pre-fit RMSD: `{receptor_prefit_rmsd:.6g} Å`",
        f"- Receptor post-fit all-atom RMSD: `{receptor_postfit_rmsd:.9g} Å`",
        f"- Receptor transform QC: `{'PASS' if ca_rmsd < RECEPTOR_RMSD_LIMIT and backbone_rmsd < RECEPTOR_RMSD_LIMIT and heavy_rmsd < RECEPTOR_RMSD_LIMIT else 'FAIL'}`",
        "",
        "## B. Ligand integrity",
        "",
        f"- DINP atoms before/after: `{len(ligand_atoms)}` / `{len(ligand_atoms)}`",
        f"- DINP heavy atoms before/after: `{len(ligand_heavy)}` / `{len(ligand_heavy)}`",
        f"- Ligand internal pairwise-distance RMSD (pre-serialization rigid transform): `{ligand_internal_rmsd:.9g} Å` (threshold `<1e-5 Å`)",
        f"- Ligand internal geometry changed: `{'YES' if ligand_internal_rmsd >= LIGAND_INTERNAL_RMSD_LIMIT else 'NO'}`",
        "- Ligand processing: no redocking, minimization, embedding, conformer regeneration, coordinate reordering, or bond-order inference.",
        "- PDB serialization note: source PDBQT uses repeated element-only atom names (`C`/`O`); output PDB uses unique names (`C01`, `O08`, ...) so BioPython/MDAnalysis retain every atom. The original names are preserved in the REMARK mapping and transform JSON.",
        "",
        "## C. Receptor–ligand relative pose preservation",
        "",
        f"- Pocket definition: `{len(pocket_residues)} receptor residues with any DINP heavy atom within 5.0 Å in the original pose`",
        f"- Minimum ligand–pocket heavy-atom distance before/after: `{source_min_pocket_distance:.6g}` / `{transformed_min_pocket_distance:.6g} Å`",
        f"- Ligand COM → pocket COM distance before/after: `{source_com_pocket_distance:.6g}` / `{transformed_com_pocket_distance:.6g} Å`",
        f"- Ligand–pocket pairwise distance-matrix RMSD: `{pocket_pairwise_rmsd:.9g} Å` (threshold `<1e-4 Å`)",
        f"- Ten nearest contact distances before: `{[round(x, 4) for x in nearest_contacts_before]}`",
        f"- Ten nearest contact distances after: `{[round(x, 4) for x in nearest_contacts_after]}`",
        f"- Ligand–pocket relative pose preserved: `{'PASS' if pocket_pairwise_rmsd < POCKET_RMSD_LIMIT else 'FAIL'}`",
        "",
        "## D. PPM membrane sanity check",
        "",
        f"- Parsed PPM DUM boundary records: `{len(dum_z)}`",
        f"- Membrane center Z: `{membrane_center_z:.3f} Å`",
        f"- Hydrophobic slab / boundary Z: `{membrane_lower_z:.3f} to {membrane_upper_z:.3f} Å`",
        f"- Half-thickness inferred from DUM records: `{membrane_half_thickness:.3f} Å`",
        f"- Transformed DINP ligand COM Z: `{ligand_com_z:.3f} Å`",
        f"- DINP near receptor pocket: `{'YES' if ligand_near_pocket else 'NO'}`",
        "- Interpretation: geometric sanity check only; no membrane-binding energy or biological conclusion is inferred.",
        "",
        "## E. Output re-read checks",
        "",
        f"- Usable protein-only complex BioPython reread: `{'PASS' if biopython_complex['readable'] else 'FAIL'} ({biopython_complex.get('atoms')})`",
        f"- Usable protein-only complex MDAnalysis reread: `{'PASS' if mdanalysis_complex['readable'] else 'FAIL'} ({mdanalysis_complex.get('atoms')} atoms)`",
        f"- Native-DUM complex BioPython reread: `{'PASS' if biopython_native_complex['readable'] else 'FAIL'}`",
        f"- Native-DUM complex MDAnalysis reread: `{'PASS' if mdanalysis_native_complex['readable'] else 'parser-incompatible with PPM DUM records'}`",
        f"- PPM native DUM information preserved: `YES` in `{TARGET_RECEPTOR.name}` and `{NATIVE_DUM_COMPLEX_PDB.name}`",
        "- The primary complex omits only PPM DUM annotation records so it is readable by standard protein/ligand tools; receptor coordinates are unchanged.",
        "",
        "## Boundary",
        "",
        "This task performed coordinate transfer and QC only. No docking, OpenMM, membrane construction, POPC insertion, or molecular dynamics was run.",
    ]
    QC_MD.write_text("\n".join(report) + "\n", encoding="utf-8")

    if not qc_pass:
        fail(f"pose-transfer QC failed; inspect {QC_MD}")

    print("PTGER4-DINP PPM POSE TRANSFER: PASS")
    print(f"Source receptor: {source_receptor}")
    print(f"PPM receptor: {TARGET_RECEPTOR}")
    print(f"Source DINP pose: {ligand_source}")
    print(f"Matched receptor atoms: {len(common_keys)}")
    print(f"C-alpha RMSD: {ca_rmsd:.9g} Å")
    print(f"Backbone RMSD: {backbone_rmsd:.9g} Å")
    print(f"All-heavy-atom RMSD: {heavy_rmsd:.9g} Å")
    print(f"Ligand internal geometry preserved: {'PASS' if ligand_internal_rmsd < LIGAND_INTERNAL_RMSD_LIMIT else 'FAIL'}")
    print(f"Ligand-pocket relative pose preserved: {'PASS' if pocket_pairwise_rmsd < POCKET_RMSD_LIMIT else 'FAIL'}")
    print(f"Ligand-pocket distance RMSD: {pocket_pairwise_rmsd:.9g} Å")
    print(f"Ligand COM Z: {ligand_com_z:.6f} Å")
    print("Ready for membrane-system build: YES")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"PTGER4-DINP PPM POSE TRANSFER: FAIL\nERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
