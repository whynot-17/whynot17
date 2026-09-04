#!/usr/bin/env bash
set -Eeuo pipefail

# Run the official PPM 3.0 source-code build locally for the frozen PTGER4
# coordinate file.  This script performs membrane orientation only; it does
# not prepare an MD system, add a ligand, or start OpenMM.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
INPUT_PDB="${REPO_ROOT}/analysis/dinp_crc_structural_pipeline/inputs/ppm/9JQZ_PTGER4_for_PPM.pdb"
OUTPUT_ROOT="${REPO_ROOT}/analysis/dinp_crc_structural_pipeline/outputs/ppm_local"
RUN_DIR="${OUTPUT_ROOT}/9JQZ_PTGER4"
BUILD_LOG="${OUTPUT_ROOT}/PPM3_BUILD_LOG.txt"
SOURCE_URL="https://cggit.cc.lehigh.edu/biomembhub/ppm3_server_code.git"
RUNTIME_ROOT="${PPM3_RUNTIME_ROOT:-/mnt/d/whynot17/ppm3_runtime}"
SOURCE_DIR="${PPM3_SOURCE_DIR:-${RUNTIME_ROOT}/ppm3_server_code_wsl}"
PPM_BIN="${SOURCE_DIR}/immers"
QC_PYTHON="${PPM3_QC_PYTHON:-/mnt/d/whynot17/dinp_stage2_md_env_wsl/bin/python}"

EXPECTED_INPUT_SHA256="97459ccad1210e79eff1907742ca5468ee143d143705c7038fc5052efd039ba3"
ORIENTED_PDB="${OUTPUT_ROOT}/9JQZ_PTGER4_PPM_oriented.pdb"
QC_REPORT="${OUTPUT_ROOT}/9JQZ_PPM_LOCAL_QC.md"
QC_JSON="${OUTPUT_ROOT}/9JQZ_PPM_LOCAL_QC.json"
PROVENANCE="${OUTPUT_ROOT}/PPM3_SOURCE_PROVENANCE.txt"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

command -v sha256sum >/dev/null 2>&1 || die "sha256sum is required"
command -v make >/dev/null 2>&1 || die "make is required"
command -v gfortran >/dev/null 2>&1 || die "gfortran is required; install it in Ubuntu-24.04 before running this script"

[[ -f "${INPUT_PDB}" ]] || die "input PDB not found: ${INPUT_PDB}"
[[ -d "${SOURCE_DIR}" ]] || die "official PPM source directory not found: ${SOURCE_DIR}"

mkdir -p "${OUTPUT_ROOT}" "${RUN_DIR}"

input_sha256="$(sha256sum "${INPUT_PDB}" | awk '{print tolower($1)}')"
[[ "${input_sha256}" == "${EXPECTED_INPUT_SHA256}" ]] || die "input PDB SHA256 mismatch: ${input_sha256}"

source_commit="$(git -C "${SOURCE_DIR}" rev-parse HEAD 2>/dev/null || true)"
[[ -n "${source_commit}" ]] || die "could not determine official PPM source commit"
# Compilation creates untracked *.o/immers files in the source directory.  The
# safety check therefore rejects substantive tracked edits, not expected build
# output or the CRLF/LF checkout normalization performed by Windows Git.
source_dirty="$(git -C "${SOURCE_DIR}" diff --ignore-space-at-eol --name-only 2>/dev/null; git -C "${SOURCE_DIR}" diff --cached --ignore-space-at-eol --name-only 2>/dev/null)"
[[ -z "${source_dirty}" ]] || die "official PPM source tree has tracked edits; refusing to run a modified source tree"

clone_date="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
{
  echo "PPM 3.0 source provenance"
  echo "source_url=${SOURCE_URL}"
  echo "source_commit=${source_commit}"
  echo "source_branch=$(git -C "${SOURCE_DIR}" branch --show-current 2>/dev/null || true)"
  echo "recorded_utc=${clone_date}"
  echo "runtime_root=${RUNTIME_ROOT}"
  echo "compiler=$(gfortran --version | head -n 1)"
  echo "make=$(make --version | head -n 1)"
} > "${PROVENANCE}"

# Compile the unmodified official source with its own Makefile.  The build
# log is deliberately replaced on every invocation so it remains auditable.
{
  echo "PPM 3.0 official local build"
  echo "build_started_utc=${clone_date}"
  echo "source_url=${SOURCE_URL}"
  echo "source_commit=${source_commit}"
  echo "compiler=$(gfortran --version | head -n 1)"
  echo "make=$(make --version | head -n 1)"
  echo
  make -C "${SOURCE_DIR}" load
  echo
  echo "build_exit_code=0"
} > "${BUILD_LOG}" 2>&1 || {
  build_rc=$?
  echo "build_exit_code=${build_rc}" >> "${BUILD_LOG}"
  die "official PPM 3.0 build failed; inspect ${BUILD_LOG}"
}

[[ -x "${PPM_BIN}" ]] || die "PPM executable was not produced: ${PPM_BIN}"

# Work in a dedicated output directory so all files written by the official
# Fortran executable (datasub1/datapar*, native *_out.pdb, etc.) are retained.
rm -f "${RUN_DIR}/9JQZ_PTGER4_for_PPMout.pdb" "${ORIENTED_PDB}" \
  "${RUN_DIR}/ppm3_stdout.txt" "${RUN_DIR}/ppm3_stderr.txt"
cp --preserve=mode,timestamps "${INPUT_PDB}" "${RUN_DIR}/9JQZ_PTGER4_for_PPM.pdb"
cp --preserve=mode,timestamps "${SOURCE_DIR}/res.lib" "${RUN_DIR}/res.lib"

PARAM_FILE="${RUN_DIR}/ppm3_input.inp"
printf '%s\n' \
  '2' \
  'no' \
  '9JQZ_PTGER4_for_PPM.pdb' \
  '1' \
  'PMm' \
  'planar' \
  'out' \
  '' > "${PARAM_FILE}"

set +e
(cd "${RUN_DIR}" && "${PPM_BIN}" < "${PARAM_FILE}" > "ppm3_stdout.txt" 2> "ppm3_stderr.txt")
ppm_exit_code=$?
set -e
echo "ppm_exit_code=${ppm_exit_code}" > "${RUN_DIR}/ppm3_exit_code.txt"
[[ "${ppm_exit_code}" -eq 0 ]] || die "PPM 3.0 execution failed (exit ${ppm_exit_code}); inspect ${RUN_DIR}/ppm3_stderr.txt"

mapfile -t generated_outputs < <(find "${RUN_DIR}" -maxdepth 1 -type f -name '*out.pdb' -printf '%f\n' | sort)
[[ "${#generated_outputs[@]}" -eq 1 ]] || die "expected exactly one native PPM oriented *_out.pdb, found ${#generated_outputs[@]}"
cp --preserve=mode,timestamps "${RUN_DIR}/${generated_outputs[0]}" "${ORIENTED_PDB}"

[[ -s "${ORIENTED_PDB}" ]] || die "PPM oriented PDB is empty"

if [[ ! -x "${QC_PYTHON}" ]]; then
  QC_PYTHON="$(command -v python3 || true)"
fi
[[ -n "${QC_PYTHON}" && -x "${QC_PYTHON}" ]] || die "no usable Python interpreter for orientation QC"

"${QC_PYTHON}" - "${INPUT_PDB}" "${ORIENTED_PDB}" "${QC_REPORT}" "${QC_JSON}" \
  "${RUN_DIR}/ppm3_stdout.txt" "${RUN_DIR}/ppm3_stderr.txt" "${PROVENANCE}" "${ppm_exit_code}" <<'PY'
import json
import math
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

input_pdb, oriented_pdb, qc_report, qc_json, stdout_path, stderr_path, provenance_path, exit_code = sys.argv[1:]
exit_code = int(exit_code)

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import is_aa
except Exception as exc:
    raise SystemExit(f"BioPython is required for QC: {exc}")

parser = PDBParser(QUIET=True)
input_structure = parser.get_structure("input", input_pdb)
oriented_structure = parser.get_structure("oriented", oriented_pdb)

def protein_residues(structure):
    rows = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=False):
                    rows.append((chain.id, residue.id, residue))
        break
    return rows

def atom_count(structure):
    return sum(1 for atom in structure.get_atoms())

def protein_atom_count(structure):
    return sum(1 for chain in next(structure.get_models()) for residue in chain for atom in residue if is_aa(residue, standard=False))

def ca_map(structure):
    result = {}
    for chain_id, residue_id, residue in protein_residues(structure):
        if residue.has_id("CA"):
            result[(chain_id, residue_id)] = residue["CA"].coord.astype(float)
    return result

def kabsch_rmsd(a, b):
    # Return optimal rigid-alignment RMSD and aligned coordinates for a -> b.
    import numpy as np
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a_center = a.mean(axis=0)
    b_center = b.mean(axis=0)
    ac = a - a_center
    bc = b - b_center
    covariance = ac.T @ bc
    u, _, vt = np.linalg.svd(covariance)
    d = np.sign(np.linalg.det(u @ vt))
    correction = np.eye(3)
    correction[-1, -1] = d
    rotation = u @ correction @ vt
    aligned = ac @ rotation + b_center
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - b) ** 2, axis=1))))
    return rmsd, aligned

def pairwise_distance_rmsd(a, b):
    import numpy as np
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    da = np.sqrt(np.maximum(0.0, ((a[:, None, :] - a[None, :, :]) ** 2).sum(axis=2)))
    db = np.sqrt(np.maximum(0.0, ((b[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)))
    tri = np.triu_indices(len(a), k=1)
    return float(np.sqrt(np.mean((da[tri] - db[tri]) ** 2)))

input_residues = protein_residues(input_structure)
oriented_residues = protein_residues(oriented_structure)
input_ca = ca_map(input_structure)
oriented_ca = ca_map(oriented_structure)
common = sorted(set(input_ca).intersection(oriented_ca), key=str)
if len(common) < 250:
    raise SystemExit(f"fewer than 250 common protein C-alpha residues for rigid QC: {len(common)}")

input_coords = [input_ca[key] for key in common]
oriented_coords = [oriented_ca[key] for key in common]
rigid_rmsd, aligned_input = kabsch_rmsd(input_coords, oriented_coords)
geometry_rmsd = pairwise_distance_rmsd(input_coords, oriented_coords)

def chains(structure):
    return [chain.id for model in structure for chain in model]

input_chain_ids = chains(input_structure)
oriented_chain_ids = chains(oriented_structure)
oriented_protein_chain_ids = sorted({chain_id for chain_id, _, _ in oriented_residues})
residue_numbers = [residue.id[1] for _, residue_id, residue in input_residues for residue in [residue]]
oriented_numbers = [residue.id[1] for _, residue_id, residue in oriented_residues for residue in [residue]]

def residue_label(row):
    return f"{row[0]}:{row[2].get_resname().strip()} {row[1][1]}"

first_input = input_residues[0]
last_input = input_residues[-1]
first_oriented = oriented_residues[0]
last_oriented = oriented_residues[-1]

stdout = Path(stdout_path).read_text(errors="replace")
stderr = Path(stderr_path).read_text(errors="replace")

number_pattern = r"([-+]?(?:\d+(?:\.\d*)?|\.\d+))"

def find_value(patterns):
    for pattern in patterns:
        match = re.search(pattern, stdout, flags=re.IGNORECASE | re.MULTILINE)
        if match:
            return match.group(1).strip()
    return "not reported; see raw PPM stdout"

membrane_thickness = find_value([
    rf"slab thickn\.\s*=\s*{number_pattern}",
    rf"hydrophobic thickness\s*=\s*{number_pattern}",
    rf"membrane thickness\s*=\s*{number_pattern}",
    rf"\bthickn\s*=\s*{number_pattern}",
])
transfer_energy = find_value([
    rf"transfer energy\s*=\s*{number_pattern}",
    rf"\bemin\s*=\s*{number_pattern}",
    rf"energy\s*=\s*{number_pattern}",
])
tm_segments = re.findall(r"^\s*\d+\s+([A-Za-z])\s+(-?\d+)\s+(-?\d+)\s+\d+\s*$", stdout, flags=re.MULTILINE)
first_oriented_ca_z = float(first_oriented[2]["CA"].coord[2])
last_oriented_ca_z = float(last_oriented[2]["CA"].coord[2])
topology_qc = {
    "requested_n_terminus": "out",
    "observed_n_terminus_ca_z_angstrom": first_oriented_ca_z,
    "observed_c_terminus_ca_z_angstrom": last_oriented_ca_z,
    "n_terminus_on_positive_z_extracellular_side": first_oriented_ca_z > 0.0,
    "c_terminus_on_negative_z_side": last_oriented_ca_z < 0.0,
    "seven_tm_segments_detected": len(tm_segments) == 7,
    "pass": first_oriented_ca_z > 0.0 and last_oriented_ca_z < 0.0 and len(tm_segments) == 7,
}

def try_mdanalysis(path):
    try:
        import MDAnalysis as mda
        universe = mda.Universe(path)
        return {"readable": True, "atoms": int(len(universe.atoms)), "residues": int(len(universe.residues)), "error": ""}
    except Exception as exc:
        return {"readable": False, "atoms": None, "residues": None, "error": str(exc)}

# PPM's native *_out.pdb intentionally appends HETATM DUM records describing
# membrane boundary points.  BioPython reads that native file, but some
# MDAnalysis versions reject the mixed dummy-record layout.  Preserve and QC
# the native output unchanged, while using a temporary ATOM-only view for the
# requested protein-coordinate MDAnalysis re-read.
native_mda_qc = try_mdanalysis(oriented_pdb)
protein_only_qc_pdb = Path(oriented_pdb).with_name(".9JQZ_PTGER4_PPM_oriented_protein_only_for_qc.pdb")
native_lines = Path(oriented_pdb).read_text(errors="replace").splitlines()
protein_only_lines = [line for line in native_lines if line.startswith(("ATOM  ", "TER", "END"))]
protein_only_qc_pdb.write_text("\n".join(protein_only_lines) + "\n")
try:
    mda_qc = try_mdanalysis(protein_only_qc_pdb)
finally:
    protein_only_qc_pdb.unlink(missing_ok=True)
ppm_dummy_atom_count = sum(
    1 for line in native_lines
    if line.startswith("HETATM") and line[17:20].strip() == "DUM"
)
raw_coordinate_record_count = sum(
    1 for line in native_lines if line.startswith(("ATOM  ", "HETATM"))
)
input_raw_coordinate_record_count = sum(
    1 for line in Path(input_pdb).read_text(errors="replace").splitlines()
    if line.startswith(("ATOM  ", "HETATM"))
)
input_sha = __import__("hashlib").sha256(Path(input_pdb).read_bytes()).hexdigest()
oriented_sha = __import__("hashlib").sha256(Path(oriented_pdb).read_bytes()).hexdigest()

payload = {
    "pdb_id": "9JQZ",
    "source_url": "https://files.rcsb.org/download/9JQZ.cif",
    "ppm_source_provenance": Path(provenance_path).read_text(errors="replace").splitlines(),
    "ppm_exit_code": exit_code,
    "input_pdb": str(Path(input_pdb).resolve()),
    "oriented_pdb": str(Path(oriented_pdb).resolve()),
    "input_sha256": input_sha,
    "oriented_sha256": oriented_sha,
    "input_chain_ids": input_chain_ids,
    "oriented_chain_ids": oriented_chain_ids,
    "oriented_protein_chain_ids": oriented_protein_chain_ids,
    "input_protein_residues": len(input_residues),
    "oriented_protein_residues": len(oriented_residues),
    "input_atoms": atom_count(input_structure),
    "oriented_atoms": atom_count(oriented_structure),
    "input_protein_atoms": protein_atom_count(input_structure),
    "oriented_protein_atoms": protein_atom_count(oriented_structure),
    "ppm_dummy_atom_count": ppm_dummy_atom_count,
    "input_raw_coordinate_record_count": input_raw_coordinate_record_count,
    "native_raw_coordinate_record_count": raw_coordinate_record_count,
    "common_ca_residues": len(common),
    "first_input_residue": residue_label(first_input),
    "last_input_residue": residue_label(last_input),
    "first_oriented_residue": residue_label(first_oriented),
    "last_oriented_residue": residue_label(last_oriented),
    "input_ca_rigid_alignment_rmsd_angstrom": rigid_rmsd,
    "internal_pairwise_distance_rmsd_angstrom": geometry_rmsd,
    "membrane_type": "PMm (mammalian plasma membrane)",
    "curvature": "planar",
    "topology_n_terminus": "out",
    "include_heteroatoms": False,
    "membrane_thickness_raw": membrane_thickness,
    "transfer_energy_raw": transfer_energy,
    "tm_segments_raw": tm_segments,
    "topology_qc": topology_qc,
    "biopython_readable": True,
    "mdanalysis": mda_qc,
    "native_mdanalysis": native_mda_qc,
    "mdanalysis_qc_view": "temporary protein-only ATOM/TER/END view; native PPM output preserved",
    "stderr": stderr.strip(),
    "internal_geometry_preserved": geometry_rmsd < 1e-3,
    "generated_utc": datetime.now(timezone.utc).isoformat(),
}

if len(oriented_protein_chain_ids) != 1:
    raise SystemExit(f"oriented output contains {len(oriented_protein_chain_ids)} protein chains, expected one")
if len(oriented_residues) < 250:
    raise SystemExit(f"oriented output contains only {len(oriented_residues)} protein residues")
if len(tm_segments) != 7:
    raise SystemExit(f"PPM did not report the expected 7 transmembrane segments: {len(tm_segments)}")
if not topology_qc["pass"]:
    raise SystemExit(
        "PPM orientation conflicts with requested N-terminal-out topology: "
        f"N-term CA z={first_oriented_ca_z:.3f}, C-term CA z={last_oriented_ca_z:.3f}, "
        f"TM segments={len(tm_segments)}"
    )
if not mda_qc["readable"]:
    raise SystemExit(f"MDAnalysis could not reread protein-only oriented QC view: {mda_qc['error']}")
if not payload["internal_geometry_preserved"]:
    raise SystemExit(f"internal geometry changed beyond tolerance: {geometry_rmsd:.6g} A")

Path(qc_json).write_text(json.dumps(payload, indent=2) + "\n")

def fmt(value):
    return str(value)

lines = [
    "# 9JQZ PPM 3.0 Local Orientation QC",
    "",
    "## Provenance",
    "",
    "- PDB ID: `9JQZ`",
    "- Official source: `https://cggit.cc.lehigh.edu/biomembhub/ppm3_server_code.git`",
    f"- PPM source commit: `{Path(provenance_path).read_text().split('source_commit=', 1)[1].splitlines()[0]}`",
    f"- PPM exit code: `{exit_code}`",
    f"- Input SHA256: `{input_sha}`",
    f"- Oriented output SHA256: `{oriented_sha}`",
    "",
    "## Run settings",
    "",
    "- Membranes: `1`",
    "- Membrane type: `PMm` — mammalian plasma membrane",
    "- Curvature: `planar`",
    "- N-terminal topology: `out`",
    "- Include heteroatoms: `no`",
    "- PPM input mode: official CLI type `2` coordinate-file mode",
    "",
    "## Coordinate QC",
    "",
    f"- Protein chains before/after: `{input_chain_ids}` / `{oriented_protein_chain_ids}`",
    f"- Protein residues before/after: `{len(input_residues)}` / `{len(oriented_residues)}`",
    f"- Parsed protein coordinate atoms before/after: `{protein_atom_count(input_structure)}` / `{protein_atom_count(oriented_structure)}`",
    f"- Raw PDB coordinate records before/native-after: `{input_raw_coordinate_record_count}` / `{raw_coordinate_record_count}`",
    f"- PPM membrane-boundary dummy atoms retained in native output: `{ppm_dummy_atom_count}`",
    f"- First residue before/after: `{residue_label(first_input)}` / `{residue_label(first_oriented)}`",
    f"- Last residue before/after: `{residue_label(last_input)}` / `{residue_label(last_oriented)}`",
    f"- Common C-alpha residues: `{len(common)}`",
    f"- C-alpha RMSD after optimal rigid alignment: `{rigid_rmsd:.6g} Å`",
    f"- Internal pairwise-distance RMSD: `{geometry_rmsd:.6g} Å`",
    f"- Oriented N-terminal CA z-coordinate: `{first_oriented_ca_z:.3f} Å` (positive-z/outside side)",
    f"- Oriented C-terminal CA z-coordinate: `{last_oriented_ca_z:.3f} Å` (negative-z side)",
    f"- PPM-reported TM segments: `{len(tm_segments)}`; 7TM GPCR body: `{'PASS' if len(tm_segments) == 7 else 'FAIL'}`",
    f"- N-terminal-out topology check: `{'PASS' if topology_qc['pass'] else 'FAIL'}`",
    f"- Internal receptor geometry preserved: `{'PASS' if payload['internal_geometry_preserved'] else 'FAIL'}`",
    "- Coordinate manipulation by this workflow: rigid-body orientation only; no loop completion, residue insertion, minimization, or renumbering.",
    "",
    "## PPM-native outputs",
    "",
    f"- Hydrophobic/membrane thickness parsed from stdout: `{membrane_thickness}`",
    f"- Transfer energy parsed from stdout: `{transfer_energy}`",
    f"- TM segment records parsed from stdout: `{tm_segments if tm_segments else 'not parsed; see raw stdout'}`",
    "- Raw stdout: `9JQZ_PTGER4/ppm3_stdout.txt`",
    "- Raw stderr: `9JQZ_PTGER4/ppm3_stderr.txt`",
    "- Native PPM oriented output: `9JQZ_PTGER4/9JQZ_PTGER4_for_PPMout.pdb`",
    "",
    "## Re-read checks",
    "",
    "- BioPython reread: `PASS`",
    f"- MDAnalysis reread of temporary protein-only oriented view: `{'PASS' if mda_qc['readable'] else 'FAIL'}" + (f" ({mda_qc['atoms']} atoms, {mda_qc['residues']} residues)" if mda_qc['readable'] else f" — {mda_qc['error']}") + "`",
    f"- MDAnalysis reread of native mixed PPM file: `{'PASS' if native_mda_qc['readable'] else 'not compatible with this parser'}" + (f" ({native_mda_qc['atoms']} atoms, {native_mda_qc['residues']} residues)" if native_mda_qc['readable'] else f" — {native_mda_qc['error']}") + "`",
    "- Native PPM output retained unchanged; its DUM records are membrane-boundary annotations, not receptor residues.",
    f"- One PTGER4 protein chain: `{'PASS' if len(oriented_protein_chain_ids) == 1 else 'FAIL'}`",
    f"- At least 250 protein residues: `{'PASS' if len(oriented_residues) >= 250 else 'FAIL'}`",
    f"- N-terminal topology constraint submitted as `out`: `PASS`",
    "- Membrane orientation conflict: `not detected by local post-run checks; inspect raw PPM output before downstream MD`",
    "",
    "## Downstream boundary",
    "",
    "This directory contains PPM orientation output only. No DINP pose was added, no POPC system was built, and no OpenMM/MD simulation was run.",
]
Path(qc_report).write_text("\n".join(lines) + "\n")
PY

echo "PPM3 LOCAL DEPLOYMENT: PASS"
echo "PPM version: ${source_commit}"
echo "Input: ${INPUT_PDB}"
echo "Oriented PDB: ${ORIENTED_PDB}"
echo "Membrane settings: PMm / planar / N-ter out / heteroatoms no"
echo "QC report: ${QC_REPORT}"
echo "Input SHA256: ${input_sha256}"
echo "Ready for PTGER4-DINP pose transfer: YES"
