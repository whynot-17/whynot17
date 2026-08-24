"""Validate the official source H5AD against completed Census query caches."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "work" / "mcop_phase2g" / "census_staged_cache"
OUTPUT = ROOT / "outputs" / "mcop_phase2g_source_h5ad_qc.json"
H5AD = Path("/mnt/d/cellxgene_census/2025-11-08/16023185-de21-4c0d-a9c8-73abdd52d142.h5ad")
EXPECTED_BYTES = 3_759_324_000
CORE = ["PPARA", "PPARD", "PPARG", "NR1I2", "NR1I3", "NR1H2", "NR1H3", "RELA", "STAT3"]


def main() -> None:
    adata = ad.read_h5ad(H5AD, backed="r")
    try:
        obs = adata.obs
        var = adata.raw.var
        donor_counts = obs.loc[obs["is_primary_data"].eq(True), "donor_id"].astype(str).value_counts()
        feature = var["feature_name"].astype(str)
        gene_cols = [int(np.flatnonzero(feature.eq(g).to_numpy())[0]) for g in CORE]
        c106_mask = obs["is_primary_data"].eq(True) & obs["donor_id"].astype(str).eq("C106")
        source_core = adata.raw.X[c106_mask.to_numpy(), :][:, gene_cols]
        source_nnz = int(source_core.nnz)
        source_sum = float(source_core.sum())
    finally:
        adata.file.close()

    cache_rows = []
    for path in sorted((CACHE / "stage1_C106").glob("raw_chunk_*.csv.gz")):
        cache_rows.append(pd.read_csv(path))
    cached = pd.concat(cache_rows, ignore_index=True)
    manifest_checks = []
    for path in sorted(CACHE.glob("stage3_*_manifest.json")):
        item = json.loads(path.read_text(encoding="utf-8"))
        if item.get("status") != "success":
            continue
        donor = str(item["donor"])
        observed = int(donor_counts.get(donor, 0))
        expected = int(item["n_obs"])
        manifest_checks.append({"donor": donor, "cache_n_obs": expected, "source_n_obs": observed, "exact": observed == expected})

    payload = {
        "source_h5ad": str(H5AD),
        "file_bytes": int(H5AD.stat().st_size),
        "expected_bytes": EXPECTED_BYTES,
        "file_size_exact": int(H5AD.stat().st_size) == EXPECTED_BYTES,
        "source_shape": [int(adata.n_obs), int(adata.n_vars)],
        "source_raw_shape": [int(adata.raw.n_obs), int(adata.raw.n_vars)],
        "all_primary": bool(adata.obs["is_primary_data"].eq(True).all()),
        "core_genes_unique": bool(all(int(feature.eq(g).sum()) == 1 for g in CORE)),
        "c106_source_n_obs": int(c106_mask.sum()),
        "c106_cache_n_obs": int(json.loads((CACHE / "stage1_C106_manifest.json").read_text(encoding="utf-8"))["n_obs"]),
        "c106_source_core_nnz": source_nnz,
        "c106_cache_core_nnz": int(len(cached)),
        "c106_source_core_sum": source_sum,
        "c106_cache_core_sum": float(pd.to_numeric(cached["soma_data"]).sum()),
        "c106_nnz_exact": source_nnz == len(cached),
        "c106_sum_exact": bool(np.isclose(source_sum, pd.to_numeric(cached["soma_data"]).sum(), rtol=0, atol=1e-6)),
        "successful_stage3_manifests": len(manifest_checks),
        "stage3_donor_counts_all_exact": bool(manifest_checks and all(x["exact"] for x in manifest_checks)),
        "stage3_donor_checks": manifest_checks,
    }
    payload["qc_pass"] = bool(
        payload["file_size_exact"]
        and payload["all_primary"]
        and payload["core_genes_unique"]
        and payload["c106_source_n_obs"] == payload["c106_cache_n_obs"]
        and payload["c106_nnz_exact"]
        and payload["c106_sum_exact"]
        and payload["stage3_donor_counts_all_exact"]
    )
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in payload.items() if k != "stage3_donor_checks"}, indent=2))
    if not payload["qc_pass"]:
        raise SystemExit("Source H5AD equivalence QC failed; formal Phase 2G analysis was not started.")


if __name__ == "__main__":
    main()
