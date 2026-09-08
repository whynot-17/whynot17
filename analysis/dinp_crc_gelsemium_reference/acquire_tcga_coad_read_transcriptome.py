"""Acquire a reproducible TCGA-COAD/READ transcriptome for WGCNA.

The Xena matrix is queried in probe chunks, cached outside Git, and reduced
only after full-transcriptome MAD values have been calculated. The fresh
41-gene set is never used to select the network genes; it is overlaid later by
the R WGCNA runner.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import xenaPython as xena
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Install xenaPython==1.0.14 in the runtime used for this script") from exc


ROOT = Path(__file__).resolve().parents[2]
FRESH_DIR = ROOT / "analysis" / "dinp_crc_gelsemium_reference"
OUTPUT = FRESH_DIR / "outputs" / "tcga_coad_read_wgcna"
GENE_FILE = FRESH_DIR / "outputs" / "query_41_genes.csv"
XENA_HUB = "https://toil.xenahubs.net"
EXPRESSION_DATASET = "TcgaTargetGtex_rsem_gene_tpm"
PHENOTYPE_DATASET = "TcgaTargetGTEX_phenotype.txt"
PROBEMAP_URL = "https://toil.xenahubs.net/download/probeMap/gencode.v23.annotation.gene.probemap"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-genes", type=int, default=8000)
    parser.add_argument("--chunk-size", type=int, default=300)
    parser.add_argument(
        "--runtime-dir",
        type=Path,
        default=Path("D:/whynot17/work/tcga_coad_read_wgcna_runtime"),
    )
    parser.add_argument("--refresh", action="store_true", help="ignore cached Xena chunks")
    return parser.parse_args()


def load_samples() -> pd.DataFrame:
    rows = []
    for cohort in ["coad", "read"]:
        path = FRESH_DIR / "outputs" / f"tcga_{cohort}_bulk_41_genes" / f"tcga_{cohort}_sample_manifest.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing prior bulk sample manifest: {path}")
        frame = pd.read_csv(path)
        frame = frame[frame["include"].astype(bool) & frame["group"].astype(str).str.endswith("_primary_tumor")].copy()
        frame["cohort"] = cohort.upper()
        rows.append(frame[["sample_id", "cohort", "group", "disease", "sample_type"]])
    metadata = pd.concat(rows, ignore_index=True)
    if metadata["sample_id"].duplicated().any():
        raise ValueError("Duplicate TCGA sample IDs across COAD and READ")
    if set(metadata["cohort"]) != {"COAD", "READ"}:
        raise ValueError("Both COAD and READ primary tumor cohorts are required")
    if metadata["sample_id"].nunique() != len(metadata):
        raise ValueError("TCGA primary tumor sample IDs are not unique")
    return metadata


def load_probemap(runtime_dir: Path) -> pd.DataFrame:
    path = runtime_dir / "gencode.v23.annotation.gene.probemap.tsv"
    if not path.exists():
        import urllib.request

        with urllib.request.urlopen(PROBEMAP_URL, timeout=60) as response:
            path.write_bytes(response.read())
    probemap = pd.read_csv(path, sep="\t", dtype=str)
    required = {"id", "gene"}
    if not required.issubset(probemap.columns):
        raise ValueError(f"Probemap lacks required columns {required}; columns={list(probemap.columns)}")
    probemap = probemap[["id", "gene"]].rename(columns={"id": "probe_id", "gene": "gene_symbol"})
    probemap["gene_symbol"] = probemap["gene_symbol"].fillna("").astype(str).str.strip()
    return probemap


def mad(values: np.ndarray) -> np.ndarray:
    medians = np.nanmedian(values, axis=1)
    return np.nanmedian(np.abs(values - medians[:, None]), axis=1)


def main() -> None:
    args = parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    args.runtime_dir.mkdir(parents=True, exist_ok=True)
    metadata = load_samples()
    samples = metadata["sample_id"].tolist()
    probemap = load_probemap(args.runtime_dir)
    probes = [probe for probe in xena.dataset_field(XENA_HUB, EXPRESSION_DATASET) if probe != "sampleID"]
    if len(probes) < 50000:
        raise ValueError(f"Unexpectedly small Xena probe universe: {len(probes)}")
    probemap = probemap[probemap["probe_id"].isin(probes)].copy()
    if probemap["probe_id"].duplicated().any():
        raise ValueError("Xena probemap contains duplicate probe IDs")
    mapping = dict(zip(probemap["probe_id"], probemap["gene_symbol"]))

    rows = []
    chunk_records = []
    total = (len(probes) + args.chunk_size - 1) // args.chunk_size
    for chunk_no, start in enumerate(range(0, len(probes), args.chunk_size), start=1):
        chunk_probes = probes[start : start + args.chunk_size]
        cache = args.runtime_dir / f"chunk_{chunk_no:04d}.npz"
        if cache.exists() and not args.refresh:
            with np.load(cache, allow_pickle=False) as loaded:
                cached_probes = loaded["probes"].astype(str).tolist()
                values = loaded["values"].astype(float)
            if cached_probes != chunk_probes or values.shape != (len(chunk_probes), len(samples)):
                cache.unlink()
                values = None
        else:
            values = None
        if values is None:
            raw = xena.dataset_fetch(XENA_HUB, EXPRESSION_DATASET, samples, chunk_probes)
            values = np.asarray(raw, dtype=float)
            if values.shape != (len(chunk_probes), len(samples)):
                raise ValueError(f"Chunk {chunk_no} returned shape {values.shape}, expected {(len(chunk_probes), len(samples))}")
            np.savez_compressed(cache, probes=np.asarray(chunk_probes, dtype="U"), values=values)
        chunk_mad = mad(values)
        for probe_id, value in zip(chunk_probes, chunk_mad):
            rows.append({"probe_id": probe_id, "gene_symbol": mapping.get(probe_id, ""), "mad": float(value)})
        chunk_records.append({"chunk": chunk_no, "start": start, "n_probes": len(chunk_probes), "cache": str(cache)})
        if chunk_no == 1 or chunk_no % 10 == 0 or chunk_no == total:
            print(f"Transcriptome chunks: {chunk_no}/{total}", flush=True)

    probe_stats = pd.DataFrame(rows)
    probe_stats = probe_stats[np.isfinite(probe_stats["mad"]) & probe_stats["gene_symbol"].ne("")].copy()
    probe_stats = probe_stats.sort_values(["mad", "gene_symbol"], ascending=[False, True])
    probe_stats = probe_stats.drop_duplicates("gene_symbol", keep="first")
    selected = probe_stats.head(args.top_genes).copy()
    if len(selected) < args.top_genes:
        raise ValueError(f"Only {len(selected)} unique mapped genes available, expected {args.top_genes}")
    selected["selection_rank"] = np.arange(1, len(selected) + 1)
    selected_probe_ids = set(selected["probe_id"])

    matrix = np.full((len(samples), len(selected)), np.nan, dtype=float)
    gene_order = selected["gene_symbol"].tolist()
    probe_to_col = dict(zip(selected["probe_id"], range(len(selected))))
    for chunk_no, start in enumerate(range(0, len(probes), args.chunk_size), start=1):
        cache = args.runtime_dir / f"chunk_{chunk_no:04d}.npz"
        with np.load(cache, allow_pickle=False) as loaded:
            cached_probes = loaded["probes"].astype(str).tolist()
            values = loaded["values"].astype(float)
        for row_no, probe_id in enumerate(cached_probes):
            if probe_id in selected_probe_ids:
                matrix[:, probe_to_col[probe_id]] = values[row_no, :]
    if np.isnan(matrix).any():
        raise ValueError("Selected transcriptome matrix contains missing values")

    expression = pd.DataFrame(matrix, index=samples, columns=gene_order)
    expression.index.name = "sample_id"
    matrix_path = args.runtime_dir / "tcga_coad_read_expression_top_variable.csv.gz"
    expression.to_csv(matrix_path, compression="gzip")
    metadata.to_csv(OUTPUT / "tcga_coad_read_wgcna_sample_manifest.csv", index=False)
    selected.to_csv(OUTPUT / "tcga_coad_read_wgcna_selected_genes.csv", index=False)
    pd.DataFrame(chunk_records).to_csv(args.runtime_dir / "tcga_coad_read_wgcna_chunk_manifest.csv", index=False)
    probe_stats.to_csv(args.runtime_dir / "tcga_coad_read_wgcna_probe_stats.csv.gz", index=False, compression="gzip")

    target_genes = pd.read_csv(GENE_FILE)["gene_symbol"].astype(str).str.strip().drop_duplicates().tolist()
    present_targets = sorted(set(target_genes) & set(gene_order))
    manifest = {
        "analysis": "TCGA-COAD + READ full-transcriptome input acquisition for WGCNA",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(ROOT)),
        "xena_hub": XENA_HUB,
        "expression_dataset": EXPRESSION_DATASET,
        "phenotype_dataset": PHENOTYPE_DATASET,
        "probemap_url": PROBEMAP_URL,
        "probemap_path": str((args.runtime_dir / "gencode.v23.annotation.gene.probemap.tsv").resolve()),
        "xena_probe_count": len(probes),
        "mapped_unique_gene_count": int(len(probe_stats)),
        "selected_top_variable_genes": int(len(selected)),
        "samples": int(len(samples)),
        "coad_primary_tumors": int((metadata["cohort"] == "COAD").sum()),
        "read_primary_tumors": int((metadata["cohort"] == "READ").sum()),
        "chunk_size": args.chunk_size,
        "chunk_count": total,
        "runtime_dir": str(args.runtime_dir.resolve()),
        "matrix_path": str(matrix_path.resolve()),
        "matrix_sha256": sha256_file(matrix_path),
        "gene_input": str(GENE_FILE.relative_to(ROOT)),
        "gene_input_sha256": sha256_file(GENE_FILE),
        "fresh_41_gene_count": len(target_genes),
        "fresh_41_genes_present_in_network_input": len(present_targets),
        "fresh_41_genes_absent_from_top_variable_network_input": sorted(set(target_genes) - set(present_targets)),
        "expression_scale": "Xena-delivered log2(TPM+0.001), as documented in dataset metadata",
        "selection_rule": "full Xena probe universe -> gene-symbol deduplication by MAD -> top variable genes; fresh 41-gene set not used for selection",
        "python": platform.python_version(),
    }
    (OUTPUT / "tcga_coad_read_wgcna_input_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"samples": len(samples), "selected_genes": len(selected), "fresh_targets_present": len(present_targets), "matrix": str(matrix_path), "manifest": str(OUTPUT / "tcga_coad_read_wgcna_input_manifest.json")}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
