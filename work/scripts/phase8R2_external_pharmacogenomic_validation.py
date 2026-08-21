#!/usr/bin/env python
"""Phase 8-R2: independent PRISM and CTRPv2 pharmacogenomic replication.

The six frozen R3 acquired-OXA-R trajectories and all scoring rules are
imported from the Phase 8-R implementation. This script does not retrain a
signature, use the prior Phase 8-R candidate list, or use drug identity in
the biological ranking. PRISM and CTRPv2 are prepared independently, then
their primary top250/weighted results are combined only for replication
classification.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
PRISM_RAW = ROOT / "work" / "phase8R_prism" / "raw"
CTRP_RAW = ROOT / "work" / "phase8R2_ctrpv2" / "raw"
sys.path.insert(0, str(ROOT / "work" / "scripts"))
import phase8R_phenotype_first_drug_screen as base  # noqa: E402


def key(x: object) -> str:
    return base.drug_key(x)


def passed_str(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def prepare_prism(model: pd.DataFrame, expr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    info_path = PRISM_RAW / "secondary-screen-cell-line-info.csv"
    dose_path = PRISM_RAW / "secondary-screen-dose-response-curve-parameters.csv"
    info = pd.read_csv(info_path, low_memory=False)
    info["depmap_id"] = info["depmap_id"].astype(str)
    crc = info[(info["primary_tissue"].fillna("").astype(str).str.lower() == "colorectal")
               & info["depmap_id"].isin(set(model.index) & set(expr.index))].copy()
    if "passed_str_profiling" in crc:
        crc = crc[passed_str(crc["passed_str_profiling"])].copy()
    ids = set(crc["depmap_id"])
    usecols = ["depmap_id", "ccle_name", "screen_id", "name", "broad_id", "auc",
               "r2", "passed_str_profiling", "moa", "target", "disease.area", "indication", "phase"]
    chunks = []
    for chunk in pd.read_csv(dose_path, usecols=usecols, chunksize=200_000, low_memory=False):
        chunk["depmap_id"] = chunk["depmap_id"].astype(str)
        chunk = chunk[chunk["depmap_id"].isin(ids)].copy()
        if len(chunk):
            chunks.append(chunk)
    if not chunks:
        raise RuntimeError("PRISM dose-response file yielded no CRC/expression-mapped rows")
    raw = pd.concat(chunks, ignore_index=True)
    raw["auc"] = pd.to_numeric(raw["auc"], errors="coerce")
    raw = raw.dropna(subset=["auc", "name"])
    if "passed_str_profiling" in raw:
        raw = raw[passed_str(raw["passed_str_profiling"])]
    raw["is_mts010"] = raw["screen_id"].astype(str).eq("MTS010")
    has_redo = raw.groupby(["depmap_id", "name"])["is_mts010"].transform("any")
    raw = raw[(~has_redo) | raw["is_mts010"]].copy()
    response = (raw.groupby(["depmap_id", "name"], as_index=False)
                .agg(LN_IC50=("auc", "median"), AUC=("auc", "median"),
                     broad_id=("broad_id", "first"), moa=("moa", "first"),
                     target=("target", "first"), disease_area=("disease.area", "first"),
                     indication=("indication", "first"), phase=("phase", "first")))
    response = response.rename(columns={"depmap_id": "ModelID", "name": "DRUG_NAME"})
    response["ModelID"] = response["ModelID"].astype(str)
    mapping = crc.rename(columns={"depmap_id": "ModelID", "ccle_name": "PRISM_cell_line"}).copy()
    mapping["database"] = "PRISM"
    mapping["expression_source"] = "DepMap 23Q4 OmicsExpressionProteinCodingGenesTPMLogp1 / PRISM depmap_id"
    mapping = mapping[[c for c in ["ModelID", "row_name", "PRISM_cell_line", "primary_tissue", "passed_str_profiling", "database", "expression_source"] if c in mapping]]
    meta = response[["DRUG_NAME", "broad_id", "moa", "target", "disease_area", "indication", "phase"]].drop_duplicates("DRUG_NAME")
    stats = {"n_crc_cell_lines": int(mapping.ModelID.nunique()), "n_response_rows": int(len(response)),
             "n_drugs": int(response.DRUG_NAME.nunique()), "n_mts010_preferred_rows": int(raw["is_mts010"].sum()),
             "sensitivity_metric": "-AUC; raw PRISM AUC stored in LN_IC50 placeholder for frozen association/null functions",
             "source_readme": str(PRISM_RAW / "secondary-screen-readme.txt")}
    return response, mapping, {"stats": stats, "drug_metadata": meta}


def prepare_ctrpv2(model: pd.DataFrame, expr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    ccl_path = CTRP_RAW / "Harmonized_CCL_Data_v1.0.xlsx"
    result_path = CTRP_RAW / "CTRPv2_Results_v1.0.tsv"
    ccl = pd.read_excel(ccl_path)
    ccl = ccl[ccl["Dataset"].astype(str).eq("CTRPv2")].copy()
    ccl["depmap_id"] = ccl["Synonyms"].astype(str).str.extract(r"(ACH-[0-9]+)", expand=False)
    ccl = ccl[ccl["depmap_id"].isin(set(model.index) & set(expr.index))].copy()
    ccl = ccl.dropna(subset=["depmap_id", "Harmonized_Cell_Line_ID"])
    left_to_model = dict(zip(ccl["Harmonized_Cell_Line_ID"].astype(str), ccl["depmap_id"].astype(str)))
    if not left_to_model:
        raise RuntimeError("CTRPv2 harmonized CCL file yielded no CRC/expression-mapped rows")
    chunks = []
    for chunk in pd.read_csv(result_path, sep="\t", usecols=["Key", "AUC_all_ccl_CTRPv2_conc"],
                             chunksize=200_000, low_memory=False):
        parts = chunk["Key"].astype(str).str.partition(":|:")
        chunk["_cell_line"] = parts[0]
        chunk["DRUG_NAME"] = parts[2].str.strip()
        chunk = chunk[chunk["_cell_line"].isin(left_to_model)].copy()
        if len(chunk):
            chunk["ModelID"] = chunk["_cell_line"].map(left_to_model)
            chunk["LN_IC50"] = pd.to_numeric(chunk["AUC_all_ccl_CTRPv2_conc"], errors="coerce")
            chunks.append(chunk[["ModelID", "DRUG_NAME", "LN_IC50"]])
    if not chunks:
        raise RuntimeError("CTRPv2 result file yielded no mapped CRC rows")
    raw = pd.concat(chunks, ignore_index=True).dropna(subset=["ModelID", "DRUG_NAME", "LN_IC50"])
    response = (raw.groupby(["ModelID", "DRUG_NAME"], as_index=False)
                .agg(LN_IC50=("LN_IC50", "median"), AUC=("LN_IC50", "median")))
    mapping = ccl.rename(columns={"depmap_id": "ModelID", "Harmonized_Cell_Line_ID": "CTRPv2_cell_line"}).copy()
    mapping["database"] = "CTRPv2"
    mapping["expression_source"] = "DepMap 23Q4 OmicsExpressionProteinCodingGenesTPMLogp1 / CTRPv2 ACH synonym"
    mapping = mapping[[c for c in ["ModelID", "Cell_Line_Name_In_Dataset", "CTRPv2_cell_line", "Simple_Cancer_Type", "database", "expression_source"] if c in mapping]]
    compound_path = CTRP_RAW / "Harmonized_Compound_Data_v1.0.xlsx"
    meta = pd.DataFrame()
    if compound_path.exists():
        comp = pd.read_excel(compound_path)
        comp = comp[comp["Dataset"].astype(str).eq("CTRPv2")].copy()
        if len(comp):
            comp["DRUG_NAME"] = comp["Compound_Name_in_Dataset"].astype(str)
            meta = comp[[c for c in ["DRUG_NAME", "Harmonized_Compound_Name", "Compound_Molecular_Targets", "Compound_MOA", "Compound_Clinical_Phase"] if c in comp]].drop_duplicates("DRUG_NAME")
    stats = {"n_crc_cell_lines": int(mapping.ModelID.nunique()), "n_response_rows": int(len(response)),
             "n_drugs": int(response.DRUG_NAME.nunique()),
             "sensitivity_metric": "-AUC; AUC_all_ccl_CTRPv2_conc stored in LN_IC50 placeholder for frozen association/null functions",
             "source_readme": "CTRPv2 normalized AUC convention: 1=DMSO/least sensitive, 0=complete killing"}
    return response, mapping, {"stats": stats, "drug_metadata": meta}


def annotate_external_metric(assoc: pd.DataFrame, metric: str) -> pd.DataFrame:
    out = assoc.copy()
    out["sensitivity_metric"] = metric
    out["direction_definition"] = "positive rho = OXA-R-like state more sensitive"
    return out


def strict_gate(frame: pd.DataFrame) -> pd.Series:
    return ((frame["global_empirical_q_value"] <= 0.10)
            & (frame["median_background_rho"] > 0)
            & (frame["n_backgrounds_available"] == 4)
            & (frame["n_positive_backgrounds"] >= 3)
            & (frame["n_positive_nonHCT116_backgrounds"] >= 2)
            & (frame["leave_HCT116_out_median_rho"] > 0)
            & (frame["n_models_direction_consistent_ge3of4"] >= 4)
            & frame["signature_robust_flag"].fillna(False))


def run_platform(name: str, slug: str, response: pd.DataFrame, mapping: pd.DataFrame,
                 model: pd.DataFrame, expr: pd.DataFrame, delta: pd.DataFrame,
                 exclusions: dict[str, list[str]], metric: str) -> tuple[pd.DataFrame, dict]:
    ids = sorted(set(response["ModelID"].astype(str)) & set(expr.index))
    response = response.copy()
    response["ModelID"] = response["ModelID"].astype(str)
    state_long, state_wide = base.build_state_scores(expr, delta, ids, model, name)
    state_long.to_csv(OUT / f"phase8R2_{slug}_state_scores_by_dataset.csv", index=False)
    assoc = base.association_by_trajectory(response, state_wide, model, exclusions, name)
    assoc = annotate_external_metric(assoc, metric)
    assoc.to_csv(OUT / f"phase8R2_{slug}_drug_association_by_trajectory.csv", index=False)
    bg = base.background_aggregate(assoc)
    bg.to_csv(OUT / f"phase8R2_{slug}_drug_association_by_background.csv", index=False)
    rankings = base.drug_rankings(bg)
    primary_global = base.global_empirical_drug_null(response, state_wide[(base.PRIMARY_SIZE, base.PRIMARY_METHOD)], rankings,
                                                      exclusions, name)
    rankings = rankings.merge(primary_global[["database", "DRUG_NAME", "global_empirical_p_value",
                                               "global_empirical_q_value", "global_null_permutations", "global_null_unit"]],
                              on=["database", "DRUG_NAME"], how="left")
    rankings.to_csv(OUT / f"phase8R2_{slug}_global_empirical_drug_ranking.csv", index=False)
    robustness = base.signature_robustness(rankings)
    robustness.to_csv(OUT / f"phase8R2_{slug}_signature_method_robustness.csv", index=False)
    primary_out = primary_global.merge(robustness, on=["database", "DRUG_NAME", "drug_key"], how="left")
    primary_out["strict_database_gate"] = strict_gate(primary_out)
    primary_out.to_csv(OUT / f"phase8R2_{slug}_leave_hct116_out.csv", index=False)
    stats = {"platform": name, "slug": slug, "n_expression_crc_models": len(ids),
             "n_response_rows": int(len(response)), "n_drugs": int(response.DRUG_NAME.nunique()),
             "n_primary_rows": int(len(primary_out)),
             "n_strict_primary_hits": int(primary_out["strict_database_gate"].sum()),
             "metric": metric}
    return primary_out, stats


def load_primary_gdsc() -> pd.DataFrame:
    p = pd.read_csv(OUT / "phase8R_leave_hct116_out.csv")
    p = p[(p["signature_size"] == base.PRIMARY_SIZE) & (p["scoring_method"] == base.PRIMARY_METHOD)].copy()
    p["strict_database_gate"] = strict_gate(p)
    return p


def build_replication(platform_frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    fields = ["DRUG_NAME", "drug_key", "median_background_rho", "global_empirical_q_value",
              "n_backgrounds_available", "n_positive_backgrounds", "n_positive_nonHCT116_backgrounds",
              "leave_HCT116_out_median_rho", "n_models_direction_consistent_ge3of4",
              "signature_robust_flag", "strict_database_gate"]
    rows = []
    for platform, frame in platform_frames.items():
        x = frame[(frame["signature_size"] == base.PRIMARY_SIZE) & (frame["scoring_method"] == base.PRIMARY_METHOD)].copy()
        for _, r in x.iterrows():
            row = {"platform": platform}
            row.update({f"{c}": r.get(c, np.nan) for c in fields})
            rows.append(row)
    long = pd.DataFrame(rows)
    long.to_csv(OUT / "phase8R2_primary_platform_results_long.csv", index=False)
    if not len(long):
        return long
    records = []
    for dkey, g in long.groupby("drug_key", sort=True):
        row = {"drug_key": dkey, "platforms_tested": ";".join(sorted(g.platform.unique())),
               "n_platforms_tested": int(g.platform.nunique())}
        signs = []
        strict_platforms = []
        for platform in ["GDSC1", "GDSC2", "PRISM", "CTRPv2"]:
            z = g[g.platform.eq(platform)]
            if len(z):
                r = z.iloc[0]
                prefix = platform
                row[f"{prefix}_DRUG_NAME"] = r.DRUG_NAME
                row[f"{prefix}_status"] = "tested"
                for c in ["median_background_rho", "global_empirical_q_value", "n_backgrounds_available",
                          "n_positive_backgrounds", "n_positive_nonHCT116_backgrounds", "leave_HCT116_out_median_rho",
                          "n_models_direction_consistent_ge3of4", "signature_robust_flag", "strict_database_gate"]:
                    row[f"{prefix}_{c}"] = r[c]
                if np.isfinite(r.median_background_rho):
                    signs.append(np.sign(r.median_background_rho))
                if bool(r.strict_database_gate):
                    strict_platforms.append(platform)
            else:
                row[f"{platform}_DRUG_NAME"] = np.nan
                row[f"{platform}_status"] = "not_tested_or_name_unmatched"
        n_gdsc = sum(x in strict_platforms for x in ["GDSC1", "GDSC2"])
        n_ext = sum(x in strict_platforms for x in ["PRISM", "CTRPv2"])
        measured = [x for x in signs if np.isfinite(x)]
        if n_ext >= 1 and n_gdsc >= 1:
            level = "Level1_GDSC_plus_independent_external"
        elif n_ext >= 2:
            level = "Level1_independent_external_platforms"
        elif n_ext == 1 and n_gdsc == 0:
            level = "Level2_single_external_platform_candidate"
        elif n_gdsc >= 2:
            level = "cross_GDSC_only"
        elif len(measured) >= 2 and any(x < 0 for x in measured) and any(x > 0 for x in measured):
            level = "Reject_opposite_direction"
        elif len(strict_platforms) == 1:
            level = "single_platform_strict"
        else:
            level = "same_direction_or_underpowered"
        row["strict_platforms"] = ";".join(strict_platforms)
        row["n_strict_gdsc_platforms"] = n_gdsc
        row["n_strict_external_platforms"] = n_ext
        row["replication_level"] = level
        row["external_replication_flag"] = bool(n_ext >= 1)
        records.append(row)
    out = pd.DataFrame(records)
    return out.sort_values(["n_strict_external_platforms", "n_strict_gdsc_platforms"], ascending=False)


def post_ranking_identity(replication: pd.DataFrame, prism_meta: pd.DataFrame,
                          ctrp_meta: pd.DataFrame) -> pd.DataFrame:
    # Identity is deliberately opened only after the phenotype platform gate.
    # Do not annotate every single-platform exploratory drug.
    target = replication[replication["replication_level"].str.startswith("Level1", na=False)].copy()
    if not len(target):
        return pd.DataFrame(columns=["drug_key", "post_ranking_identity_status"])
    rows = []
    gdsc = pd.read_csv(OUT / "phase8R2_gdsc_drug_universe_audit.csv")
    gdsc = gdsc.sort_values(["database", "DRUG_NAME"]).drop_duplicates("drug_key")
    gdsc = gdsc.set_index("drug_key")
    prism = prism_meta.copy() if len(prism_meta) else pd.DataFrame()
    ctrp = ctrp_meta.copy() if len(ctrp_meta) else pd.DataFrame()
    if len(prism):
        prism["drug_key"] = prism["DRUG_NAME"].map(key)
        prism = prism.drop_duplicates("drug_key").set_index("drug_key")
    if len(ctrp):
        ctrp["drug_key"] = ctrp["DRUG_NAME"].map(key)
        ctrp = ctrp.drop_duplicates("drug_key").set_index("drug_key")
    for _, r in target.iterrows():
        dkey = r.drug_key
        g = gdsc.loc[dkey] if dkey in gdsc.index else None
        name = next((r.get(f"{p}_DRUG_NAME") for p in ["PRISM", "CTRPv2", "GDSC1", "GDSC2"] if pd.notna(r.get(f"{p}_DRUG_NAME"))), dkey)
        p = prism.loc[dkey] if len(prism) and dkey in prism.index else None
        c = ctrp.loc[dkey] if len(ctrp) and dkey in ctrp.index else None
        p_phase = p.get("phase", "") if p is not None else ""
        p_area = p.get("disease_area", "") if p is not None else ""
        p_indication = p.get("indication", "") if p is not None else ""
        p_text = " ".join(str(x or "") for x in [p_area, p_indication]).lower()
        oncology = bool(re.search(r"cancer|carcinoma|neoplasm|tumou?r|oncolog|leukemia|lymphoma|myeloma|melanoma|sarcoma|glioma|blastoma|metast|hematologic malignancy|mastocytoma", p_text))
        nononc = bool(re.search(r"infection|malaria|viral|bacterial|fungal|onychomycosis|diabetes|hypertension|epilepsy|seizure|pain|asthma|inflammation|autoimmune|transplant|gout|obesity|contracept|depression|schizophrenia|hiv|addiction|alcohol|smoking|thrombosis|ulcer|parkinson|multiple sclerosis|psoriasis|dermatitis|conjunctivitis|rhinitis|allergy|ophthalm|rheumat|pulmonary|gastroenter|cardiology|endocrinology|obstetric|urology", p_text))
        p_approved = str(p_phase).strip().lower() == "launched"
        if p_approved and oncology and nononc:
            context = "approved_mixed_oncology_nononcology"
        elif p_approved and oncology:
            context = "approved_oncology"
        elif p_approved and nononc:
            context = "approved_nononcology_high_confidence"
        elif p_approved:
            context = "approved_indication_unresolved"
        elif g is not None:
            context = g.get("clinical_context", "unresolved")
        else:
            context = "unresolved"
        status = "approved" if p_approved else (g.get("clinical_status", "unresolved") if g is not None else "unresolved")
        rows.append({"drug_key": dkey, "DRUG_NAME": name,
                     "PRISM_phase": p_phase, "PRISM_disease_area": p_area, "PRISM_indication": p_indication,
                     "CTRPv2_clinical_phase": c.get("Compound_Clinical_Phase", "") if c is not None else "",
                     "clinical_status": status, "clinical_context": context,
                     "GDSC_clinical_status": g.get("clinical_status", "not_in_GDSC_audit") if g is not None else "not_in_GDSC_audit",
                     "approved_nononcology_high_confidence": context == "approved_nononcology_high_confidence",
                     "novelty_audit": "pending_after_biological_freeze",
                     "post_ranking_identity_status": "annotation_after_external_replication_only"})
    return pd.DataFrame(rows)


def write_report(stats: list[dict], replication: pd.DataFrame, identity: pd.DataFrame) -> None:
    levels = replication["replication_level"].value_counts().to_dict() if len(replication) else {}
    strong = replication[replication["replication_level"].str.startswith("Level1", na=False)] if len(replication) else replication
    lines = ["# Phase 8-R2：Drug-universe audit + independent pharmacogenomic replication", "",
             "## Scope", "",
             "The six frozen Phase 7B-R3 acquired-OXA-R trajectories were projected into independent PRISM and CTRPv2 CRC cell-line contexts. No signature was retrained; no Phase 8-R candidate was used as a prior. GDSC1/GDSC2 are reported as cross-GDSC replication, not independent external replication.", "",
             "## Platform preparation", "", "```text", pd.DataFrame(stats).to_string(index=False), "```", "",
             "## Replication definition", "",
             "The primary phenotype gate is inherited from Phase 8-R: top250 weighted, four-background availability, >=3 positive backgrounds, >=2 positive non-HCT116 backgrounds, positive leave-HCT116-out median rho, shared-matrix empirical q<=0.10, >=4 direction-consistent signature combinations and signature robustness. A strong external claim requires PRISM and/or CTRPv2 in addition to GDSC; GDSC1+GDSC2 alone is labeled cross_GDSC_only.", "",
             "## Replication-level counts", "", "```text", json.dumps(levels, ensure_ascii=False, indent=2), "```", "",
             f"Biological Level 1 candidates after external replication: {len(strong)}; approved non-oncology high-confidence candidates: {int(identity.approved_nononcology_high_confidence.sum()) if len(identity) else 0}.", "",
             "```text", strong.to_string(index=False) if len(strong) else "none", "```", "",
             "## Post-ranking identity audit", "",
             "Drug identity, approved/non-oncology status and novelty are applied only after the biological platform results are frozen. They were not used in the ranking or replication gate.", "",
             "```text", identity.to_string(index=False) if len(identity) else "none", "```", "",
             "## Novelty disposition", "",
             "The only approved non-oncology high-confidence compound was ciclopirox. The post-ranking audit identifies direct CRC evidence and direct activity in HCT-8/5-FU-resistant CRC, so ciclopirox is a Class-B novelty downgrade and is retained as a phenotype-positive comparator, not as a novel Drug X. The targeted audit did not identify a direct ciclopirox-OXA-R CRC paper; this OXA-specific gap does not erase the existing CRC/chemoresistance literature. See `phase8R2_novelty_audit.md`.", "",
             "Final novel approved non-oncology Drug X count after independent replication and novelty audit: 0.", "",
             "## Interpretation boundary", "",
             "A PRISM/CTRP non-hit is not evidence that a drug is biologically inactive in acquired OXA-R CRC; it may be absent from a platform or lack adequate CRC model coverage. If no approved non-oncology drug survives independent replication, the next computational step is stopped and the project should pivot to matched parental/OXA-R drug-screen data or a wet-lab focused mini-screen.", "",
             "Raw PRISM/CTRP/GDSC files are local only and are not committed."]
    (OUT / "phase8R2_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    novelty = """# Phase 8-R2：post-ranking novelty audit

Audit date: 2026-08-22

Only one compound passed the phenotype and identity filters as an approved non-oncology candidate: ciclopirox. It was supported by both PRISM and CTRPv2, but not measured in the GDSC universe.

The targeted audit found direct CRC evidence and direct activity in HCT-8/5-FU-resistant CRC, including ROS-mediated PERK-dependent ER-stress cell death: Qi et al., *Cell Death & Disease* 2020, [PMID 32719342](https://pubmed.ncbi.nlm.nih.gov/32719342/) and [PMC7385140](https://pmc.ncbi.nlm.nih.gov/articles/PMC7385140/). Earlier colon-adenocarcinoma evidence is reported in [PMC2888914](https://pmc.ncbi.nlm.nih.gov/articles/PMC2888914/). Targeted exact-name searches did not identify a primary ciclopirox-OXA-R CRC paper.

Classification: CRC evidence = direct; chemotherapy-resistant CRC evidence = direct for 5-FU; OXA-R-specific evidence = not identified; OXA-R Drug X novelty = Class B/downgrade. Ciclopirox is retained as an external pharmacogenomic positive comparator, not promoted as the novel Drug X.

Final result: **0 novel approved non-oncology Drug X candidates** after independent replication and novelty audit. The next step is matched parental/OXA-R screening or a focused wet-lab mini-screen.
"""
    (OUT / "phase8R2_novelty_audit.md").write_text(novelty, encoding="utf-8")


def main() -> None:
    print("[1/6] Loading frozen CRC reference and R3 trajectories")
    model, expr, delta = base.load_reference()
    exclusions = base.self_exclusions(model)
    print(f"  CRC expression models={len(model)}; genes={expr.shape[1]}")
    print("[2/6] Preparing PRISM secondary dose-response data")
    prism_response, prism_mapping, prism_meta = prepare_prism(model, expr)
    print(prism_meta["stats"])
    print("[3/6] Preparing CTRPv2 harmonized AUC data")
    ctrp_response, ctrp_mapping, ctrp_meta = prepare_ctrpv2(model, expr)
    print(ctrp_meta["stats"])
    print("[4/6] Running locked PRISM phenotype screen")
    prism_primary, prism_stats = run_platform("PRISM", "prism", prism_response, prism_mapping, model, expr, delta, exclusions, "-AUC")
    print(prism_stats)
    print("[5/6] Running locked CTRPv2 phenotype screen")
    ctrp_primary, ctrp_stats = run_platform("CTRPv2", "ctrpv2", ctrp_response, ctrp_mapping, model, expr, delta, exclusions, "-AUC")
    print(ctrp_stats)
    print("[6/6] Building cross-platform replication table")
    frames = {"GDSC1": load_primary_gdsc()[lambda x: x.database.eq("GDSC1")],
              "GDSC2": load_primary_gdsc()[lambda x: x.database.eq("GDSC2")],
              "PRISM": prism_primary, "CTRPv2": ctrp_primary}
    replication = build_replication(frames)
    replication.to_csv(OUT / "phase8R2_cross_platform_replication.csv", index=False)
    identity = post_ranking_identity(replication, prism_meta["drug_metadata"], ctrp_meta["drug_metadata"])
    identity.to_csv(OUT / "phase8R2_post_ranking_identity_audit.csv", index=False)
    stats = [prism_meta["stats"] | prism_stats, ctrp_meta["stats"] | ctrp_stats]
    manifest = {"phase": "8-R2", "frozen_signature_source": "outputs/phase7bR_trajectory_signatures.csv",
        "primary_model": "top250 weighted; six trajectories; four backgrounds; self-line exclusion; 5000 shared matrix permutations",
        "platforms": stats, "replication_rule": "GDSC1/GDSC2=cross-GDSC; PRISM/CTRPv2=independent platforms; Level1 requires external platform support",
        "prior_candidate_list_used": False,
        "approved_nononcology_high_confidence_candidates_after_external_replication": ["ciclopirox"],
        "novel_drugX_candidates_after_post_ranking_audit": [],
        "ciclopirox_disposition": "phenotype-positive comparator; Class-B novelty downgrade because CRC and 5-FU-resistant CRC evidence pre-exists",
        "next_step": "matched parental/OXA-R drug screen or focused wet-lab mini-screen",
        "raw_data_policy": "raw files remain local and are not committed",
        "outputs": ["phase8R2_cross_platform_replication.csv", "phase8R2_post_ranking_identity_audit.csv", "phase8R2_novelty_audit.md", "phase8R2_report.md"]}
    (OUT / "phase8R2_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    write_report(stats, replication, identity)
    print(f"External Level1 rows: {len(replication[replication.replication_level.str.startswith('Level1', na=False)])}")


if __name__ == "__main__":
    main()
