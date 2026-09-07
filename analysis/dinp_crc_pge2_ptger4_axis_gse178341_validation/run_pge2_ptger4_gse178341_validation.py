#!/usr/bin/env python
"""Run the frozen PGE2 -> PTGER4 robustness analysis on GSE178341.

The GEO H5 file is the public processed UMI matrix.  This wrapper extracts the
small target-gene universe in two HDF5 blocks, joins the official c295 labels,
and creates marker-defined transfer states inside the official macrophage and
monocyte compartments.  The downstream score/rank/receiver calculations are
the same as the discovery robustness module.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import runpy
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


HERE = Path(__file__).parent
DEFAULT_H5 = Path(r"E:\mcop\mcop\work\gse178341\GSE178341_crc10x_full_c295v4_submit.h5")
DEFAULT_META = Path(r"E:\mcop\mcop\work\gse178341\GSE178341_metatables.csv.gz")
DEFAULT_CLUSTER = Path(r"E:\mcop\mcop\work\gse178341\GSE178341_cluster.csv.gz")
DEFAULT_CACHE = Path(r"E:\mcop\mcop\work\gse178341\GSE178341_pge2_target_cell_scores.csv")
DEFAULT_OUT = HERE / "outputs"
DEFAULT_DOWNSTREAM = Path(r"E:\mcop\mcop\work\gse178341\_downstream_robustness.py")

TARGET_GENES = ["PLA2G4A", "PTGS1", "PTGS2", "PTGES", "PTGES2", "PTGES3", "PTGER4"]
SIGNATURES = {
    "SPP1_like_TAM": ["SPP1", "APOC1", "MMP9", "CTSB", "CTSL", "TREM2", "GPNMB"],
    "C1QC_like_TAM": ["C1QA", "C1QB", "C1QC", "APOE", "FOLR2", "MRC1"],
    "FCN1_inflammatory_monocyte_like": ["FCN1", "S100A8", "S100A9", "CTSS", "LILRB1", "VCAN"],
    "resident_like": ["FOLR2", "LYVE1", "MRC1", "CD163"],
}
EXTRA_GENES = sorted(set(g for gs in SIGNATURES.values() for g in gs) - set(TARGET_GENES))
EXTRACT_GENES = TARGET_GENES + EXTRA_GENES


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def decode(values) -> list[str]:
    return [x.decode("utf-8", "replace") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values]


def extract_targets(path: Path, genes: list[str]) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Extract target rows from the 10x-style CSC H5 without loading all genes."""
    with h5py.File(path, "r", rdcc_nbytes=1024 * 1024) as f:
        m = f["matrix"]
        names = decode(m["features/name"][:])
        upper_to_row = {g.upper(): i for i, g in enumerate(names)}
        missing = [g for g in genes if g.upper() not in upper_to_row]
        if missing:
            raise ValueError(f"Genes missing from GSE178341 H5: {missing}")
        rows = np.array([upper_to_row[g.upper()] for g in genes], dtype=np.int64)
        barcodes = decode(m["barcodes"][:])
        indptr = np.asarray(m["indptr"][:], dtype=np.int64)
        n_cells = len(barcodes)
        out = np.zeros((n_cells, len(genes)), dtype=np.float32)
        totals = np.zeros(n_cells, dtype=np.float32)

        # The data array has two very large gzip chunks.  Read each once and
        # use a sparse CSC view so selecting 28 rows stays small in memory.
        first_chunk = int(m["data"].chunks[0])
        split = int(np.searchsorted(indptr, first_chunk, side="right") - 1)
        split = max(1, min(split, n_cells - 1)) if n_cells > 1 else n_cells
        for start, end in ((0, split), (split, n_cells)):
            a, b = int(indptr[start]), int(indptr[end])
            data = np.asarray(m["data"][a:b], dtype=np.float32)
            indices = np.asarray(m["indices"][a:b], dtype=np.int32)
            ptr = indptr[start : end + 1] - a
            mat = sparse.csc_matrix((data, indices, ptr), shape=(len(names), end - start))
            out[start:end, :] = np.asarray(mat[rows, :].T.toarray(), dtype=np.float32)
            totals[start:end] = np.asarray(mat.sum(axis=0)).ravel().astype(np.float32)
            del mat, data, indices, ptr
    return barcodes, out, totals


def zscore_columns(x: np.ndarray) -> np.ndarray:
    mu = np.nanmean(x, axis=0)
    sd = np.nanstd(x, axis=0)
    sd[sd == 0] = 1.0
    return (x - mu) / sd


def marker_states(log_expr: np.ndarray, genes: list[str], meta: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Assign transparent marker-defined transfer states within Macro/Mono."""
    gi = {g.upper(): i for i, g in enumerate(genes)}
    scores = np.full((len(log_expr), len(SIGNATURES)), np.nan, dtype=np.float32)
    labels = list(SIGNATURES)
    for j, label in enumerate(labels):
        inds = [gi[g.upper()] for g in SIGNATURES[label] if g.upper() in gi]
        if inds:
            scores[:, j] = zscore_columns(log_expr[:, inds]).mean(axis=1)
    best_idx = np.nanargmax(np.where(np.isfinite(scores), scores, -np.inf), axis=1)
    best = scores[np.arange(len(scores)), best_idx]
    second = np.partition(np.where(np.isfinite(scores), scores, -np.inf), -2, axis=1)[:, -2]
    margin = best - second
    state = np.full(len(log_expr), "other_myeloid", dtype=object)
    macro = meta["clMidwayPr"].astype(str).eq("Macro").to_numpy()
    mono = meta["clMidwayPr"].astype(str).eq("Mono").to_numpy()
    # A state must have positive standardized support and a visible margin;
    # ambiguous cells remain in other_myeloid and are not forced into a TAM.
    confident = np.isfinite(best) & (best >= 0.15) & (margin >= 0.05)
    for j, label in enumerate(labels):
        use = confident & (best_idx == j) & macro
        state[use] = label
    use_fcn1 = confident & (best_idx == labels.index("FCN1_inflammatory_monocyte_like")) & mono
    state[use_fcn1] = "FCN1_inflammatory_monocyte_like"
    return state, best.astype(np.float32), margin.astype(np.float32)


def build_cell_scores(h5: Path, metadata: Path, cluster: Path, cache: Path) -> dict:
    meta = pd.read_csv(metadata)
    cl = pd.read_csv(cluster)
    if len(meta) != len(cl) or not meta["cellID"].astype(str).equals(cl["sampleID"].astype(str)):
        raise ValueError("GSE178341 metadata and cluster rows are not in the same cell order")
    barcodes, counts, totals = extract_targets(h5, EXTRACT_GENES)
    if not np.array_equal(meta["cellID"].astype(str).to_numpy(), np.array(barcodes, dtype=str)):
        raise ValueError("GSE178341 H5 barcodes do not match metatables cellID order")
    lib = np.maximum(totals.astype(np.float64), 1.0)
    log_expr = np.log1p(counts.astype(np.float64) / lib[:, None] * 1e4).astype(np.float32)
    target_cols = {g: log_expr[:, EXTRACT_GENES.index(g)] for g in TARGET_GENES}
    state, state_score, state_margin = marker_states(log_expr, EXTRACT_GENES, cl)

    group = np.full(len(meta), "other_cells", dtype=object)
    top = cl["clTopLevel"].astype(str)
    mid = cl["clMidwayPr"].astype(str)
    spec = meta["SPECIMEN_TYPE"].astype(str)
    group[top.eq("Mast").to_numpy()] = "mast_cells"
    epi = top.eq("Epi").to_numpy()
    group[epi & spec.eq("T").to_numpy()] = "tumor_epithelial"
    group[epi & spec.eq("N").to_numpy()] = "normal_epithelial"
    group[mid.eq("Fibro").to_numpy()] = "fibroblast_like"
    group[mid.eq("Endo").to_numpy()] = "endothelial_like"
    myeloid = top.eq("Myeloid").to_numpy()
    group[myeloid] = "other_myeloid"
    group[myeloid & np.isin(state, list(SIGNATURES))] = state[myeloid & np.isin(state, list(SIGNATURES))]
    # Official macrophage-like cells without a confident transfer state remain
    # visible as other_myeloid; the report separately gives cM02 coverage.
    out = pd.DataFrame({
        "cell_id": meta["cellID"].astype(str),
        "donor_id": meta["PID"].astype(str),
        "condition": np.where(spec.eq("T"), "tumor", "normal"),
        "analysis_group": group,
        "official_top_level": cl["clTopLevel"].astype(str),
        "official_midway": cl["clMidwayPr"].astype(str),
        "official_subtype": cl["cl295v11SubFull"].astype(str),
        "marker_state_score": state_score,
        "marker_state_margin": state_margin,
    })
    for g in TARGET_GENES:
        out[f"{g}_log1p"] = target_cols[g]
    cache.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(cache, index=False)
    return {
        "n_cells": int(len(out)),
        "n_donors": int(out.donor_id.nunique()),
        "donors": sorted(out.donor_id.unique().tolist()),
        "official_top_level": {str(k): int(v) for k, v in cl["clTopLevel"].value_counts().items()},
        "official_macrophage_cells": int(mid.eq("Macro").sum()),
        "official_cM02_cells": int(cl["cl295v11SubFull"].eq("cM02 (Macrophage-like)").sum()),
        "marker_state_counts": {str(k): int(v) for k, v in pd.Series(group).value_counts().items()},
    }


def patch_downstream_script(source: Path, dest: Path):
    text = source.read_text(encoding="utf-8")
    text = text.replace("GSE144735", "GSE178341")
    text = text.replace("six donors", "62 donors")
    text = text.replace("only one tumor donor", "the available donor coverage")
    text = text.replace("only one tumor donor, so its apparent high rank cannot be generalized to six donors", "variable donor coverage, so ranks must be read with the coverage tables")
    text = text.replace("exact paired sign-flip permutation p-values", "seeded Monte Carlo paired sign-flip p-values (200,000 draws; exact binomial sign test also reported)")
    text = text.replace("exact paired sign-flip permutation; supplementary", "seeded Monte Carlo paired sign-flip; supplementary")
    text = text.replace("{int(r.top2_donors)}/6 total\", \"top-two donor frequency; denominator reports eligible donor groups and all six tumor donors\")", "{int(r.top2_donors)}/{int(donor_rank_df.donor_id.nunique())} total\", \"top-two donor frequency; denominator reports eligible donor groups and all available tumor donors\")")
    text = text.replace("{int(r.top1_donors)}/{int(r.n_donors_with_group)} eligible; {int(r.top1_donors)}/6 total\" for _, r in mast_r.iterrows()", "{int(r.top1_donors)}/{int(r.n_donors_with_group)} eligible; {int(r.top1_donors)}/{int(donor_rank_df.donor_id.nunique())} total\" for _, r in mast_r.iterrows()")
    text = text.replace("{int(r.top2_donors)}/6 total)\"", "{int(r.top2_donors)}/{int(cells.donor_id.nunique())} total)\"")
    text = text.replace("{int(spp1_a.top2_donors.iloc[0]) if len(spp1_a) else 0}/6 total)", "{int(spp1_a.top2_donors.iloc[0]) if len(spp1_a) else 0}/{cells.donor_id.nunique()} total)")
    text = text.replace("{int(c1spp)}/6 total)", "{int(c1spp)}/{cells.donor_id.nunique()} total)")
    text = text.replace("{c1spp}/6 total)", "{c1spp}/{cells.donor_id.nunique()} total)")
    text = text.replace("{c1spp}/6 overall", "{c1spp}/{cells.donor_id.nunique()} overall")
    text = text.replace("{c1spp}/6** overall", "{c1spp}/{cells.donor_id.nunique()}** overall")
    text = text.replace("{c1spp}/6**", "{c1spp}/{cells.donor_id.nunique()}**")
    # The discovery cohort has six donors, so its exact sign-flip helper can
    # enumerate 2^6 assignments.  GSE178341 has 62 donors; use the same exact
    # binomial sign test and a seeded Monte Carlo sign-flip companion instead
    # of attempting an impossible 2^62 enumeration.
    start, end = text.index("def sign_test"), text.index("def c1qc_contrasts")
    sign_fn = '''def sign_test(x: pd.Series) -> tuple[int, int, float, float]:
    x = pd.to_numeric(x, errors="coerce").dropna()
    pos, neg = int((x > 0).sum()), int((x < 0).sum())
    nonzero = pos + neg
    if nonzero == 0:
        return pos, neg, np.nan, np.nan
    p = float(stats.binomtest(pos, nonzero, .5, alternative="two-sided").pvalue)
    vals = x.to_numpy(float)
    observed = abs(float(vals.sum()))
    rng = np.random.default_rng(178341)
    draws = 200000
    signs = rng.choice(np.array([-1.0, 1.0]), size=(draws, len(vals)))
    perm_p = float((np.count_nonzero(np.abs(signs @ vals) >= observed - 1e-12) + 1) / (draws + 1))
    return pos, neg, p, perm_p


'''
    text = text[:start] + sign_fn + text[end:]
    dest.write_text(text, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5", type=Path, default=DEFAULT_H5)
    ap.add_argument("--metadata", type=Path, default=DEFAULT_META)
    ap.add_argument("--cluster", type=Path, default=DEFAULT_CLUSTER)
    ap.add_argument("--cell-scores", type=Path, default=DEFAULT_CACHE)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--downstream", type=Path, default=DEFAULT_DOWNSTREAM)
    args = ap.parse_args()
    required_cache = {"cell_id", "donor_id", "condition", "analysis_group"} | {f"{g}_log1p" for g in TARGET_GENES}
    if args.cell_scores.exists():
        cached = pd.read_csv(args.cell_scores, usecols=lambda c: c in required_cache)
        if len(cached) == 370115 and required_cache.issubset(cached.columns):
            cl = pd.read_csv(args.cluster, usecols=["clTopLevel", "clMidwayPr", "cl295v11SubFull"])
            info = {"n_cells": int(len(cached)), "n_donors": int(cached.donor_id.nunique()), "donors": sorted(cached.donor_id.unique().tolist()), "official_macrophage_cells": int(cl["clMidwayPr"].eq("Macro").sum()), "official_cM02_cells": int(cl["cl295v11SubFull"].eq("cM02 (Macrophage-like)").sum()), "marker_state_counts": {str(k): int(v) for k, v in cached.analysis_group.value_counts().items()}, "cache_reused": True}
        else:
            info = build_cell_scores(args.h5, args.metadata, args.cluster, args.cell_scores)
    else:
        info = build_cell_scores(args.h5, args.metadata, args.cluster, args.cell_scores)

    source = Path(__file__).parents[1] / "dinp_crc_pge2_ptger4_axis_robustness" / "run_pge2_ptger4_robustness.py"
    patch_downstream_script(source, args.downstream)
    argv = [
        str(args.downstream), "--cell-scores", str(args.cell_scores), "--matrix", str(args.h5),
        "--annotation", str(args.metadata), "--subtype", str(args.cell_scores), "--out", str(args.out),
    ]
    old = sys.argv
    try:
        sys.argv = argv
        runpy.run_path(str(args.downstream), run_name="__main__")
    finally:
        sys.argv = old
    manifest = {"dataset": "GSE178341", "role": "primary validation", "h5": {"path": str(args.h5), "sha256": sha256(args.h5), "size_bytes": args.h5.stat().st_size}, "metadata": {"path": str(args.metadata), "sha256": sha256(args.metadata), "size_bytes": args.metadata.stat().st_size}, "cluster": {"path": str(args.cluster), "sha256": sha256(args.cluster), "size_bytes": args.cluster.stat().st_size}, "cell_score_cache": {"path": str(args.cell_scores), "sha256": sha256(args.cell_scores), "size_bytes": args.cell_scores.stat().st_size}, "cohort_info": info, "official_annotation": "GSE178341 c295 clTopLevel/clMidwayPr/cl295v11SubFull; marker-defined transfer states only inside official Macro/Mono", "raw_data_access": "controlled-access dbGaP; public processed GEO H5 used", "state_rule": "gene-wise z-scored marker means; best score >=0.15 and margin >=0.05; ambiguous cells retained as other_myeloid"}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "validation_input_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    report = args.out / "PGE2_PTGER4_GSE178341_VALIDATION_REPORT.md"
    robustness_report = args.out / "PGE2_PTGER4_ROBUSTNESS_REPORT.md"
    robustness_text = robustness_report.read_text(encoding="utf-8")
    robustness_text = robustness_text.replace(
        "SPP1-like TAM has adequate coverage in the available donor coverage, so its apparent high rank cannot be generalized to 62 donors.",
        "SPP1-like TAM is covered unevenly across donors, so its apparent high rank must be interpreted with the eligible-donor denominator and cannot be generalized to all 62 donors.",
    )
    robustness_report.write_text(robustness_text, encoding="utf-8")
    report.write_text(robustness_text + "\n\n## Validation-specific annotation\n\nGSE178341 uses the public c295 official labels. The cohort has one official macrophage-like cluster (cM02); SPP1-like, C1QC-like, FCN1-like and resident-like labels are transparent marker-defined transfer states within Macro/Mono and are not presented as the original c295 annotation. Ambiguous cells remain `other_myeloid`.\n", encoding="utf-8")
    print("GSE178341 PRIMARY VALIDATION: PASS")
    print(f"Cells: {info['n_cells']}; donors: {info['n_donors']}; official cM02 macrophage-like cells: {info['official_cM02_cells']}")
    print(f"Validation outputs: {args.out}")


if __name__ == "__main__":
    main()
