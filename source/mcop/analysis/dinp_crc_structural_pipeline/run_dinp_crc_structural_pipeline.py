#!/usr/bin/env python3
"""DINP–CRC protein-network, docking, and MD preparation pipeline.

Purpose
-------
Take the macrophage-prioritized DINP–CRC candidates and prepare a reproducible
structural follow-up:

1) STRING PPI connectivity for protein-coding candidates;
2) RCSB PDB structure inventory by UniProt accession;
3) DINP ligand retrieval from PubChem;
4) optional AutoDock Vina docking when receptor/ligand PDBQT and docking boxes
   are available;
5) reproducible GROMACS/ACPYPE MD command generation for selected complexes.

The script never substitutes a related phthalate for DINP. Missing structural
or software resources are recorded rather than interpreted as biological
negative evidence.

This is a prioritization/structural-simulation pipeline, not causal proof.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import requests


DEFAULT_TARGETS = {
    # NEAT1 excluded from protein PPI/docking because it is an lncRNA.
    "MMP9": "P14780",
    "TIMP1": "P01033",
    "STAT3": "P40763",
    "PTGER4": "P35408",
    "PTGES3": "Q15185",
    "CXCR4": "P61073",
}

DEFAULT_CONFIG = {
    "species": 9606,
    "string_required_score": 400,
    "ligand_name": "diisononyl phthalate",
    "vina_exhaustiveness": 32,
    "vina_num_modes": 20,
    "md_ns": 100,
    "temperature_k": 310,
    "pressure_bar": 1.0,
    "water_model": "tip3p",
    "force_field": "amber99sb-ildn",
    # Fill only after a biologically justified pocket is selected.
    # Example:
    # "docking_boxes": {
    #   "PTGER4": {"center": [0.0, 0.0, 0.0], "size": [22.0, 22.0, 22.0]},
    # }
    "docking_boxes": {},
    # Optional manual PDB choices after structure review.
    # Example: {"PTGER4": "XXXX"}
    "manual_pdb": {},
}


@dataclass
class StructureHit:
    gene: str
    uniprot: str
    pdb_id: str
    title: str
    method: str
    resolution_angstrom: Optional[float]
    release_date: str


class Audit:
    def __init__(self) -> None:
        self.rows: List[dict] = []

    def add(self, step: str, status: str, detail: str, **kwargs) -> None:
        row = {"step": step, "status": status, "detail": detail, **kwargs}
        self.rows.append(row)
        print(f"[{status}] {step}: {detail}")

    def write(self, path: Path) -> None:
        path.write_text(json.dumps(self.rows, indent=2, ensure_ascii=False), encoding="utf-8")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def load_or_create_config(path: Path) -> dict:
    if not path.exists():
        save_json(path, DEFAULT_CONFIG)
        print(f"Created config template: {path}")
    cfg = json.loads(path.read_text(encoding="utf-8"))
    merged = dict(DEFAULT_CONFIG)
    merged.update(cfg)
    return merged


def _response_from_bytes(content: bytes, url: str) -> requests.Response:
    response = requests.Response()
    response.status_code = 200
    response.url = url
    response._content = content
    response.encoding = "utf-8"
    return response


def _curl_get(url: str, *, params=None, timeout=30) -> requests.Response:
    from urllib.parse import urlencode
    full_url = url
    if params:
        full_url += ("&" if "?" in full_url else "?") + urlencode(params)
    proc = subprocess.run(
        ["curl", "-fsSL", "--connect-timeout", "10", "--max-time", str(timeout), full_url],
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        detail = proc.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"curl GET failed (exit={proc.returncode}): {detail}")
    return _response_from_bytes(proc.stdout, full_url)


def _curl_post(url: str, *, json_body=None, data=None, timeout=30) -> requests.Response:
    body = json.dumps(json_body) if json_body is not None else data
    cmd = [
        "curl", "-fsSL", "--connect-timeout", "10", "--max-time", str(timeout),
        "-H", "Content-Type: application/json", "--data-binary", body or "", url,
    ]
    proc = subprocess.run(cmd, capture_output=True, check=False)
    if proc.returncode != 0:
        detail = proc.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"curl POST failed (exit={proc.returncode}): {detail}")
    return _response_from_bytes(proc.stdout, url)


def http_get(url: str, *, params=None, timeout=30) -> requests.Response:
    last = None
    for attempt in range(4):
        try:
            # curl is more reliable than the WSL Python TLS path for the
            # public RCSB/STRING/PubChem endpoints. Keep requests as a
            # portable fallback for environments without curl.
            if shutil.which("curl"):
                r = _curl_get(url, params=params, timeout=timeout)
            else:
                r = requests.get(url, params=params, timeout=timeout)
            if r.status_code == 429:
                time.sleep(2 ** attempt)
                continue
            r.raise_for_status()
            return r
        except requests.RequestException as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET failed after retries: {url}: {last}")


def http_post(url: str, *, json_body=None, data=None, timeout=30) -> requests.Response:
    last = None
    for attempt in range(4):
        try:
            if shutil.which("curl"):
                r = _curl_post(url, json_body=json_body, data=data, timeout=timeout)
            else:
                r = requests.post(url, json=json_body, data=data, timeout=timeout)
            if r.status_code == 429:
                time.sleep(2 ** attempt)
                continue
            r.raise_for_status()
            return r
        except requests.RequestException as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"POST failed after retries: {url}: {last}")


def run_string_ppi(targets: Dict[str, str], species: int, required_score: int,
                   out_dir: Path, audit: Audit) -> None:
    """Query STRING network endpoint and compute simple topology metrics."""
    out_dir.mkdir(parents=True, exist_ok=True)
    identifiers = "%0d".join(targets.keys())
    url = "https://string-db.org/api/tsv/network"
    params = {
        "identifiers": identifiers,
        "species": species,
        "required_score": required_score,
        "network_type": "functional",
        "caller_identity": "dinp_crc_structural_pipeline",
    }
    try:
        r = http_get(url, params=params)
        raw_path = out_dir / "string_network.tsv"
        raw_path.write_text(r.text, encoding="utf-8")
        rows = list(csv.DictReader(r.text.splitlines(), delimiter="\t"))
        degree = {g: 0 for g in targets}
        edges = set()
        for row in rows:
            a = row.get("preferredName_A", "")
            b = row.get("preferredName_B", "")
            if a in degree and b in degree and a != b:
                edge = tuple(sorted((a, b)))
                if edge not in edges:
                    edges.add(edge)
                    degree[a] += 1
                    degree[b] += 1
        with (out_dir / "string_node_degree.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["gene", "degree_within_candidate_set"])
            for gene, deg in sorted(degree.items(), key=lambda x: (-x[1], x[0])):
                w.writerow([gene, deg])
        audit.add("STRING_PPI", "PASS", f"Retrieved {len(rows)} STRING rows; {len(edges)} unique candidate-set edges")
    except Exception as exc:
        audit.add("STRING_PPI", "UNAVAILABLE", str(exc))


def rcsb_search_by_uniprot(uniprot: str) -> List[str]:
    """Search PDB entries containing the UniProt accession.

    Uses RCSB full-text search intentionally as a broad inventory step. Hits are
    later annotated; final docking structures should be manually reviewed for
    human protein identity, construct, ligand state, completeness, and pocket.
    """
    # Search the polymer-entity reference-sequence accession field rather
    # than using an unscoped full-text query.  The latter is rejected by the
    # current RCSB Search API and can also return irrelevant text matches.
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 100}},
    }
    r = http_post("https://search.rcsb.org/rcsbsearch/v2/query", json_body=query)
    payload = r.json()
    return [x["identifier"] for x in payload.get("result_set", [])]


def rcsb_entry_metadata(pdb_id: str) -> Tuple[str, str, Optional[float], str]:
    r = http_get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    x = r.json()
    title = x.get("struct", {}).get("title", "")
    methods = x.get("exptl", []) or []
    method = ";".join(m.get("method", "") for m in methods if m.get("method"))
    res = None
    values = x.get("rcsb_entry_info", {}).get("resolution_combined") or []
    if values:
        try:
            res = min(float(v) for v in values if v is not None)
        except Exception:
            res = None
    release = x.get("rcsb_accession_info", {}).get("initial_release_date", "")
    return title, method, res, release


def structure_inventory(targets: Dict[str, str], out_dir: Path, audit: Audit) -> List[StructureHit]:
    out_dir.mkdir(parents=True, exist_ok=True)
    hits: List[StructureHit] = []
    for gene, uniprot in targets.items():
        try:
            pdbs = rcsb_search_by_uniprot(uniprot)
            for pdb in pdbs:
                try:
                    title, method, res, release = rcsb_entry_metadata(pdb)
                    hits.append(StructureHit(gene, uniprot, pdb, title, method, res, release))
                except Exception as exc:
                    audit.add("RCSB_ENTRY", "WARN", f"{gene}/{pdb}: {exc}")
            audit.add("RCSB_SEARCH", "PASS", f"{gene}: {len(pdbs)} candidate PDB entries")
        except Exception as exc:
            audit.add("RCSB_SEARCH", "UNAVAILABLE", f"{gene}: {exc}")

    csv_path = out_dir / "structure_inventory.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        fields = list(asdict(StructureHit("", "", "", "", "", None, "")).keys())
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for hit in sorted(hits, key=lambda h: (h.gene, h.resolution_angstrom or math.inf, h.pdb_id)):
            w.writerow(asdict(hit))
    return hits


def choose_structures(hits: List[StructureHit], cfg: dict) -> Dict[str, str]:
    """Choose preliminary structures for download.

    Manual config overrides automatic choice. Automatic choice is the best
    reported resolution among inventory hits and is NOT considered final
    biological validation.
    """
    by_gene: Dict[str, List[StructureHit]] = {}
    for h in hits:
        by_gene.setdefault(h.gene, []).append(h)
    chosen = {}
    for gene in DEFAULT_TARGETS:
        manual = cfg.get("manual_pdb", {}).get(gene)
        if manual:
            chosen[gene] = manual.upper()
            continue
        candidates = by_gene.get(gene, [])
        if candidates:
            candidates.sort(key=lambda h: (h.resolution_angstrom is None, h.resolution_angstrom or math.inf, h.pdb_id))
            chosen[gene] = candidates[0].pdb_id
    return chosen


def download_pdbs(chosen: Dict[str, str], out_dir: Path, audit: Audit) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for gene, pdb_id in chosen.items():
        try:
            r = http_get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
            path = out_dir / f"{gene}_{pdb_id}.pdb"
            path.write_bytes(r.content)
            audit.add("PDB_DOWNLOAD", "PASS", f"{gene}: {pdb_id}", sha256=sha256_file(path))
        except Exception as exc:
            # Recent cryo-EM entries can be distributed only as mmCIF. Keep
            # the selected entry and preserve the failure of the legacy PDB
            # route instead of silently dropping the structure.
            try:
                r = http_get(f"https://files.rcsb.org/download/{pdb_id}.cif")
                path = out_dir / f"{gene}_{pdb_id}.cif"
                path.write_bytes(r.content)
                audit.add(
                    "MMCIF_DOWNLOAD", "PASS",
                    f"{gene}: {pdb_id} (mmCIF fallback after PDB unavailable: {exc})",
                    sha256=sha256_file(path),
                )
            except Exception as cif_exc:
                audit.add("PDB_DOWNLOAD", "UNAVAILABLE", f"{gene}/{pdb_id}: {exc}; mmCIF fallback: {cif_exc}")


def fetch_dinp(out_dir: Path, ligand_name: str, audit: Audit) -> Optional[Path]:
    """Resolve DINP by exact name through PubChem and download 3D SDF."""
    out_dir.mkdir(parents=True, exist_ok=True)
    from urllib.parse import quote
    try:
        name = quote(ligand_name)
        prop_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/"
            "Title,IUPACName,CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
        )
        props = http_get(prop_url).json()
        save_json(out_dir / "dinp_pubchem_properties.json", props)
        cid = props.get("PropertyTable", {}).get("Properties", [{}])[0].get("CID")
        sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/record/SDF?record_type=3d"
        three_d_error = None
        try:
            r = http_get(sdf_url)
        except Exception:
            if cid is None:
                raise
            # PubChem may resolve an exact DINP record but expose its 3-D
            # conformer only through the CID endpoint, so retry explicitly.
            try:
                r = http_get(
                    f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{int(cid)}/record/SDF?record_type=3d"
                )
            except Exception as exc:
                three_d_error = exc
                r = None
        if r is None:
            # Keep a chemically explicit 2-D source artifact for downstream
            # 3-D generation, but never return it as a docking-ready ligand.
            r2d = http_get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{int(cid)}/record/SDF?record_type=2d"
            )
            sdf2d = out_dir / "DINP_2D.sdf"
            sdf2d.write_bytes(r2d.content)
            audit.add("PUBCHEM_DINP_2D", "PASS", f"Retrieved 2-D fallback for CID {int(cid)}", sha256=sha256_file(sdf2d))
            audit.add("PUBCHEM_DINP_3D", "UNAVAILABLE", f"No PubChem 3-D conformer for CID {int(cid)}: {three_d_error}")
            return None
        sdf = out_dir / "DINP_3D.sdf"
        sdf.write_bytes(r.content)
        audit.add("PUBCHEM_DINP", "PASS", f"Retrieved exact-name ligand '{ligand_name}'", sha256=sha256_file(sdf))
        return sdf
    except Exception as exc:
        audit.add("PUBCHEM_DINP", "UNAVAILABLE", str(exc))
        return None


def command_exists(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def run_cmd(cmd: List[str], cwd: Optional[Path], log_path: Path, audit: Audit, step: str) -> bool:
    try:
        proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None, capture_output=True, text=True, check=False)
        log_path.write_text(
            "$ " + " ".join(cmd) + "\n\nSTDOUT\n" + proc.stdout + "\n\nSTDERR\n" + proc.stderr,
            encoding="utf-8",
        )
        if proc.returncode == 0:
            audit.add(step, "PASS", " ".join(cmd))
            return True
        audit.add(step, "FAIL", f"exit={proc.returncode}; see {log_path}")
        return False
    except Exception as exc:
        audit.add(step, "FAIL", str(exc))
        return False


def prepare_and_dock(cfg: dict, chosen: Dict[str, str], structures_dir: Path,
                     ligand_sdf: Optional[Path], out_dir: Path, audit: Audit) -> None:
    """Run Vina only for targets with explicit docking boxes.

    Automatic pocket guessing is intentionally avoided. A docking box must be
    supplied after structure/pocket review. This prevents uninterpretable
    whole-protein blind docking from becoming the primary result.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    boxes = cfg.get("docking_boxes", {})
    if not boxes:
        audit.add("DOCKING", "WAITING_CONFIG", "No docking_boxes configured; structure inventory completed only")
        return
    if ligand_sdf is None:
        audit.add("DOCKING", "SKIP", "DINP SDF unavailable")
        return
    needed = ["vina", "mk_prepare_ligand.py", "mk_prepare_receptor.py"]
    missing = [x for x in needed if not command_exists(x)]
    if missing:
        audit.add("DOCKING", "UNAVAILABLE", "Missing executables: " + ", ".join(missing))
        return

    ligand_pdbqt = out_dir / "DINP.pdbqt"
    run_cmd(["mk_prepare_ligand.py", "-i", str(ligand_sdf), "-o", str(ligand_pdbqt)], None,
            out_dir / "prepare_ligand.log", audit, "PREPARE_LIGAND")
    if not ligand_pdbqt.exists():
        return

    for gene, box in boxes.items():
        if gene not in chosen:
            audit.add("DOCKING", "SKIP", f"{gene}: no selected structure")
            continue
        pdb_id = chosen[gene]
        receptor_pdb = structures_dir / f"{gene}_{pdb_id}.pdb"
        if not receptor_pdb.exists():
            audit.add("DOCKING", "SKIP", f"{gene}: receptor PDB missing")
            continue
        center = box.get("center")
        size = box.get("size")
        if not (isinstance(center, list) and len(center) == 3 and isinstance(size, list) and len(size) == 3):
            audit.add("DOCKING", "SKIP", f"{gene}: invalid center/size")
            continue

        gene_dir = out_dir / gene
        gene_dir.mkdir(exist_ok=True)
        receptor_pdbqt = gene_dir / "receptor.pdbqt"
        ok = run_cmd([
            "mk_prepare_receptor.py", "--read_pdb", str(receptor_pdb), "-o", str(receptor_pdbqt),
        ], None, gene_dir / "prepare_receptor.log", audit, f"PREPARE_RECEPTOR_{gene}")
        if not ok or not receptor_pdbqt.exists():
            continue
        out_pdbqt = gene_dir / "DINP_docked.pdbqt"
        log_txt = gene_dir / "vina_scores.txt"
        vina_cmd = [
            "vina", "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_pdbqt),
            "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
            "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
            "--exhaustiveness", str(cfg.get("vina_exhaustiveness", 32)),
            "--num_modes", str(cfg.get("vina_num_modes", 20)),
            "--out", str(out_pdbqt), "--log", str(log_txt),
        ]
        run_cmd(vina_cmd, None, gene_dir / "vina_run.log", audit, f"VINA_{gene}")


def write_md_templates(cfg: dict, chosen: Dict[str, str], out_dir: Path, audit: Audit) -> None:
    """Generate per-target MD protocol templates, not fake completed simulations."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for gene, pdb_id in chosen.items():
        d = out_dir / gene
        d.mkdir(exist_ok=True)
        protocol = f"""# {gene} / DINP MD protocol template\n\nSelected RCSB entry: {pdb_id}\nRequested production length: {cfg.get('md_ns', 100)} ns\nTemperature: {cfg.get('temperature_k', 310)} K\nPressure: {cfg.get('pressure_bar', 1.0)} bar\nProtein force field: {cfg.get('force_field', 'amber99sb-ildn')}\nWater model: {cfg.get('water_model', 'tip3p')}\n\nRequired before execution:\n1. Manually validate receptor construct and binding pocket.\n2. Convert the selected docking pose into a protein–DINP complex.\n3. Generate DINP parameters with ACPYPE/GAFF (or an explicitly documented alternative).\n4. Merge ligand topology with the GROMACS protein topology.\n5. Energy minimization -> NVT -> NPT -> production MD.\n6. Report RMSD, RMSF, radius of gyration, SASA, H-bonds, and binding free energy (MM-PBSA/GBSA).\n\nSuggested commands (adapt paths/topologies after ligand parameterization):\n\n    gmx pdb2gmx -f complex.pdb -o processed.gro -water {cfg.get('water_model', 'tip3p')} -ff {cfg.get('force_field', 'amber99sb-ildn')}\n    gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt dodecahedron\n    gmx solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top\n    gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr\n    gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral\n    gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr && gmx mdrun -deffnm em\n    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr && gmx mdrun -deffnm nvt\n    gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr && gmx mdrun -deffnm npt\n    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr && gmx mdrun -deffnm md\n"""
        (d / "MD_PROTOCOL.md").write_text(protocol, encoding="utf-8")
    audit.add("MD_TEMPLATE", "PASS", f"Generated MD protocol templates for {len(chosen)} selected structures")


def write_summary(out_dir: Path, chosen: Dict[str, str], audit: Audit) -> None:
    lines = [
        "# DINP–CRC structural pipeline summary",
        "",
        "## Candidate proteins",
        "",
        "Protein-coding macrophage-prioritized candidates: MMP9, TIMP1, STAT3, PTGER4, PTGES3, CXCR4. NEAT1 is retained biologically but excluded from protein PPI/docking because it is an lncRNA.",
        "",
        "## Preliminary structure choices",
        "",
    ]
    if chosen:
        for gene, pdb in chosen.items():
            lines.append(f"- {gene}: {pdb}")
    else:
        lines.append("- No structures selected; inspect audit/structure inventory.")
    lines += [
        "",
        "## Interpretation",
        "",
        "RCSB hits and automatic resolution-based choices are an inventory/prioritization step only. Final docking requires manual validation of species, construct, pocket, missing residues, cofactors, and ligand state. Docking/MD provide structural plausibility, not evidence that DINP reaches or binds the target in vivo.",
        "",
        "## Audit",
        "",
        f"- Audit events: {len(audit.rows)}",
    ]
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="analysis/dinp_crc_structural_pipeline/outputs")
    ap.add_argument("--config", default="analysis/dinp_crc_structural_pipeline/config.json")
    ap.add_argument("--skip-network", action="store_true")
    ap.add_argument("--skip-docking", action="store_true")
    args = ap.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    config_path = Path(args.config)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    cfg = load_or_create_config(config_path)
    audit = Audit()

    save_json(out_dir / "targets.json", DEFAULT_TARGETS)
    save_json(out_dir / "run_config_snapshot.json", cfg)

    if not args.skip_network:
        run_string_ppi(DEFAULT_TARGETS, int(cfg["species"]), int(cfg["string_required_score"]), out_dir / "ppi", audit)

    hits = structure_inventory(DEFAULT_TARGETS, out_dir / "structures", audit)
    chosen = choose_structures(hits, cfg)
    save_json(out_dir / "structures" / "selected_structures.json", chosen)
    download_pdbs(chosen, out_dir / "structures" / "pdb", audit)

    ligand_sdf = fetch_dinp(out_dir / "ligand", cfg.get("ligand_name", "diisononyl phthalate"), audit)

    if not args.skip_docking:
        prepare_and_dock(
            cfg,
            chosen,
            out_dir / "structures" / "pdb",
            ligand_sdf,
            out_dir / "docking",
            audit,
        )
    else:
        audit.add("DOCKING", "SKIPPED", "--skip-docking requested; no docking box or docking result was generated")

    write_md_templates(cfg, chosen, out_dir / "md", audit)
    write_summary(out_dir, chosen, audit)
    audit.write(out_dir / "audit.json")
    print(f"Done. Outputs: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
