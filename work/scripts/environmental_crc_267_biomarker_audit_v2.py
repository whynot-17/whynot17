"""Complete X/B/D/C/T audit for all 267 Phase 1 environmental chemicals.

This script works from the full CDC NHANES environmental laboratory catalog
and downloaded XPT files, not from the two previously analyzed phthalate
analytes. Every candidate receives an identity record, a documented search
trail, an analyte mapping decision, cycle-level detectability where applicable,
and a testability assessment. Lack of a supported analyte is a resolved audit
outcome with evidence, not a generic pending bucket.
"""

from __future__ import annotations

import importlib.util
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = ROOT / "outputs"
XPT_DIR = ROOT / "work" / "nhanes_phase2a" / "environmental_xpt"
CATALOG = OUTPUTS / "environmental_crc_267_nhanes_environmental_lab_catalog.csv"
PHASE1_MATRIX = OUTPUTS / "environmental_crc_267_actionability_matrix.csv"

OUT_IDENTITY = OUTPUTS / "environmental_crc_267_identity_table_v2.csv"
OUT_MAPPING = OUTPUTS / "environmental_crc_267_biomarker_mapping.csv"
OUT_DETECT = OUTPUTS / "environmental_crc_267_detectability_by_cycle.csv"
OUT_TEST = OUTPUTS / "environmental_crc_267_testability_audit.csv"
OUT_MATRIX = OUTPUTS / "environmental_crc_267_actionability_matrix_v2.csv"
OUT_FLOW = OUTPUTS / "environmental_crc_267_actionability_flow.csv"
OUT_REVIEW = OUTPUTS / "environmental_crc_267_manual_review_queue_v2.csv"
OUT_FIREWALL = OUTPUTS / "environmental_crc_267_outcome_firewall_audit.json"
OUT_MANIFEST = OUTPUTS / "environmental_crc_267_actionability_manifest_v2.json"

SEVEN_CYCLE_CODES = {"1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"}

METAL_TO_VARIABLE = {
    "barium": "URXUBA", "cadmium": "URXUCD", "cobalt": "URXUCO", "cesium": "URXUCS",
    "molybdenum": "URXUMO", "manganese": "URXUMN", "lead": "URXUPB", "platinum": "URXUPT",
    "antimony": "URXUSB", "tin": "URXUSN", "silver": "URXUSR", "thallium": "URXUTL",
    "tungsten": "URXUTU", "uranium": "URXUUR", "mercury": "URXUHG",
}

PHENOL_TO_VARIABLE = {
    "bisphenol a": "URXBPH", "bisphenol": "URXBPH", "benzophenone-3": "URXBP3",
    "butylparaben": "URXBUP", "ethylparaben": "URXEPB", "methylparaben": "URXMPB",
    "propylparaben": "URXPPB",
}

PAH_TO_VARIABLE = {
    "naphthalene": "URXP02", "methylnaphthalene": "URXP02", "dimethylnaphthalene": "URXP02",
    "fluorene": "URXP04", "phenanthrene": "URXP25", "methylphenanthrene": "URXP25",
    "pyrene": "URXP10", "benzo": "URXP10", "chrysene": "URXP10", "anthracene": "URXP25",
    "fluoranthene": "URXP10", "acenaphth": "URXP02", "indeno": "URXP10", "perylene": "URXP10",
    "picene": "URXP10", "retene": "URXP10", "corannulene": "URXP10",
}

PHTHALATE_TO_VARIABLE = {
    "butylbenzyl phthalate": ["URXMZP"], "mono-benzyl phthalate": ["URXMZP"],
    "monobutyl phthalate": ["URXMBP"], "monoethyl phthalate": ["URXMEP"],
    "mono-(2-ethylhexyl)phthalate": ["URXMHP"], "mono(2-ethyl-5-hydroxyhexyl)phthalate": ["URXMHH"],
    "mono(2-ethyl-5-oxohexyl)phthalate": ["URXMOH"], "2-ethyl-5-carboxypentyl phthalate": ["URXECP"],
    "mono-isobutyl phthalate": ["URXMIB"], "monoisononylphthalate": ["URXMNP"],
    "mono(carboxy-isooctyl)phthalate": ["URXCOP"], "diethylhexyl phthalate": ["URXMHH", "URXMOH", "URXECP"],
    "dinonylphthalate": ["URXMNP", "URXCOP"], "diisononyl phthalate": ["URXMNP", "URXCOP"],
}

PFAS_TO_VARIABLE = {
    "perfluorooctane sulfonic acid": ["LBXNFOS"], "perfluorobutanesulfonic acid": ["LBXPFBS"],
    "perfluorohexanesulfonic acid": ["LBXPFHS"], "perfluorohexanoic acid": [],
    "perfluorodecanoic acid": ["LBXPFDE"], "perfluorononanoic acid": ["LBXPFNA"],
    "perfluoro-n-nonanoic acid": ["LBXPFNA"], "perfluoroundecanoic acid": ["LBXPFUA"],
    "perfluorododecanoic acid": ["LBXPFDO"], "perfluoroheptanoic acid": ["LBXPFHP"],
    "perfluoro-n-heptanoic acid": ["LBXPFHP"], "2-(n-methyl-pfosa) acetate": ["LBXMPAH"],
}


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def norm(text: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(text or "").lower()).strip()


def contains_any(text: str, terms: list[str]) -> bool:
    return any(term in text for term in terms)


def load_harmonized() -> pd.DataFrame:
    """Rebuild the existing CRC frame from local XPT files if pickle is incompatible."""
    mbzp = load_module(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "audit_mbzp_harmonizer")
    frames = []
    for idx, spec in enumerate(mbzp.CYCLES):
        frame, _ = mbzp.load_cycle(spec, idx)
        frames.append(frame)
    harmonized = pd.concat(frames, ignore_index=True)
    return harmonized[(harmonized["age"] >= 20) & harmonized["cancer_outcome_available"]].copy()


def paired_flag(value_col: str, columns: set[str]) -> str:
    if value_col.startswith("URX"):
        candidate = "URD" + value_col[3:] + "LC"
        if candidate in columns:
            return candidate
    if value_col.startswith("LBX"):
        for candidate in ["LBD" + value_col[3:] + "L", "LBD" + value_col[3:] + "LC"]:
            if candidate in columns:
                return candidate
    return ""


def build_registry(catalog: pd.DataFrame) -> tuple[pd.DataFrame, dict[tuple[str, str], pd.DataFrame]]:
    rows = []
    frames: dict[tuple[str, str], pd.DataFrame] = {}
    for item in catalog.itertuples(index=False):
        local = XPT_DIR / f"{item.cycle}_{item.data_file}"
        if not local.exists() or local.stat().st_size < 1000:
            continue
        frame = pd.read_sas(local, format="xport", encoding="latin1")
        frames[(item.cycle, item.data_file)] = frame
        columns = set(frame.columns)
        weight_cols = [c for c in frame.columns if c.startswith("WTS") or c.startswith("WTSA") or c.startswith("WTSB")]
        preferred_weights = [c for c in weight_cols if c in {"WTSA2YR", "WTSB2YR"}]
        weight_col = preferred_weights[0] if preferred_weights else (weight_cols[0] if weight_cols else "")
        # Some NHANES supplemental XPTs repeat analyte columns but do not
        # carry a usable survey subsample weight. They remain in the official
        # catalog, but cannot support a valid D/T or complex-survey analysis.
        if not weight_col:
            continue
        for value_col in frame.columns:
            if not (value_col.startswith("URX") or value_col.startswith("LBX")):
                continue
            if value_col in {"URXUCR"}:
                continue
            flag_col = paired_flag(value_col, columns)
            values = pd.to_numeric(frame[value_col], errors="coerce")
            valid = values.notna()
            if weight_col:
                weight = pd.to_numeric(frame[weight_col], errors="coerce")
                valid &= weight.gt(0)
            else:
                weight = pd.Series(1.0, index=frame.index)
            flag = pd.to_numeric(frame[flag_col], errors="coerce") if flag_col else pd.Series(np.nan, index=frame.index)
            above = valid & (flag.ne(1) if flag_col else True)
            rows.append({
                "cycle": item.cycle,
                "cycle_begin_year": item.cycle_begin_year,
                "laboratory_title": item.laboratory_title,
                "data_file": item.data_file,
                "data_url": item.data_url,
                "doc_url": item.doc_url,
                "local_xpt": str(local),
                "variable": value_col,
                "flag_variable": flag_col,
                "weight_variable": weight_col,
                "matrix": "urine" if value_col.startswith("URX") else "serum_or_blood",
                "n_measured": int(valid.sum()),
                "n_above_lod": int(above.sum()),
                "above_lod_pct": float(100 * above.sum() / valid.sum()) if valid.sum() else np.nan,
            })
    return pd.DataFrame(rows), frames


def identity_row(row: pd.Series) -> dict[str, object]:
    name = str(row["ChemicalName"])
    n = norm(name)
    chemical_class = str(row.get("chemical_class", ""))
    canonical = re.sub(r"\s+", " ", n).strip()
    entity_type = "specific chemical" if not re.match(r"^(D\d+|C\d+)$", str(row["ChemicalID"])) else "CTD entity"
    if chemical_class.startswith("heavy_metals"):
        entity_type = "elemental species or metal compound"
    elif "phthalate" in n:
        entity_type = "phthalate parent or monoester metabolite"
    elif "perfluoro" in n or "fluorotelomer" in n:
        entity_type = "PFAS or fluorinated compound"
    elif "bisphenol" in n:
        entity_type = "bisphenol or bisphenol-related compound"
    elif "pollutant" in n or "pesticide" in n or "disruptor" in n:
        entity_type = "umbrella exposure class"
    parent = ""
    metabolite_of = ""
    if "monoisononyl" in n or "carboxy isooctyl" in n or n in {"dinonylphthalate", "diisononyl phthalate"}:
        parent = "DINP"
        metabolite_of = "DINP"
    elif "benzyl phthalate" in n or "butylbenzyl" in n:
        parent = "BBzP"
        metabolite_of = "BBzP"
    elif "ethylhexyl" in n or "diethylhexyl" in n:
        parent = "DEHP"
        metabolite_of = "DEHP"
    aliases = [name]
    if row.get("CasRN") and str(row["CasRN"]) != "nan":
        aliases.append(str(row["CasRN"]))
    if "phthalate" in n:
        aliases.extend(["phthalate", "plasticizer"])
    if "perfluoro" in n:
        aliases.extend(["PFAS", n.replace("perfluoro", "PF")])
    if "bisphenol a" in n:
        aliases.append("BPA")
    exposure_family = chemical_class.split(";")[0] if chemical_class else "unclassified environmental chemical"
    return {
        "ChemicalID": row["ChemicalID"], "ChemicalName": name, "canonical_name": canonical,
        "synonyms": ";".join(dict.fromkeys(aliases)), "CasRN": row.get("CasRN", ""),
        "PubChemCID": row.get("PubChemCID", ""), "DTXSID": row.get("DTXSID", ""),
        "chemical_class": chemical_class, "entity_type": entity_type, "parent_compound": parent,
        "metabolite_of": metabolite_of, "exposure_family": exposure_family,
    }


def mapping_for_candidate(row: pd.Series, registry: pd.DataFrame) -> dict[str, object]:
    name = norm(row["ChemicalName"])
    raw_name = str(row["ChemicalName"]).lower()
    chemical_class = str(row.get("chemical_class", "")).lower()
    domain_hits = registry.loc[registry["laboratory_title"].str.lower().apply(lambda x: (
        ("phthalate" in x and "phthalate" in name)
        or ("metal" in x and "metal" in chemical_class)
        or ("perfluoro" in x and "perfluoro" in name)
        or ("phenol" in x and "bisphenol" in name)
        or ("paraben" in x and "bisphenol" in name)
        or ("aromatic hydrocarbon" in x and ("pah" in chemical_class or "pyrene" in name or "naphthal" in name))
        or ("pesticide" in x and "pesticide" in chemical_class)
        or ("flame retard" in x and "flame" in chemical_class)
        or ("volatile organic" in x and "voc" in chemical_class)
    ))]
    candidates: list[str] = []
    confidence = "unresolved"
    relationship = "none"
    axis = ""
    source_reason = ""

    # HFPO-DA/GenX-like replacement PFAS must not inherit the PFNA mapping
    # merely because its name contains the substring "perfluorononanoic acid".
    # NHANES LBXPFNA is a PFNA measurement, not a validated HFPO-DA assay.
    genx_like = any(token in raw_name for token in [
        "4,8-dioxa-3h-perfluorononanoic",
        "hexafluoropropylene oxide dimer acid",
        "hfpo-da",
        "hfpo da",
        "genx",
    ])

    for key, values in PHTHALATE_TO_VARIABLE.items():
        if key in raw_name:
            candidates = values
            confidence = "high" if raw_name.startswith("mono") else "moderate"
            relationship = "direct urinary metabolite" if raw_name.startswith("mono") else "parent-to-validated-urinary-metabolite axis"
            axis = "DINP-related exposure axis" if "isononyl" in raw_name or "carboxy-isooctyl" in raw_name or "dinonyl" in raw_name else ("DEHP-related exposure axis" if "ethylhexyl" in raw_name else "phthalate exposure axis")
            source_reason = "CDC NHANES phthalate/plasticizer metabolite panel and codebook"
            break
    if not candidates and ("perfluoro" in raw_name or "pfas" in chemical_class):
        if genx_like:
            confidence, relationship, axis = "unresolved", "PFAS class searched; no candidate-specific NHANES GenX/HFPO-DA variable", "PFAS exposure axis"
            source_reason = "CDC NHANES PFAS panels searched; PFNA explicitly excluded as a non-equivalent analyte"
        else:
            for key, values in PFAS_TO_VARIABLE.items():
                if key in raw_name:
                    candidates, confidence, relationship, axis = values, "high", "direct serum PFAS analyte", "PFAS exposure axis"
                    source_reason = "CDC NHANES PFAS serum panel and codebook"
                    break
        if not candidates and "fluorotelomer" not in raw_name:
            if not source_reason:
                confidence, relationship, axis = "unresolved", "PFAS class searched; no candidate-specific NHANES variable", "PFAS exposure axis"
                source_reason = "CDC NHANES PFAS panels searched; no specific analyte name match"
    if not candidates and ("bisphenol" in raw_name or "bisphenol" in chemical_class):
        for key, value in PHENOL_TO_VARIABLE.items():
            if key in raw_name:
                candidates, confidence, relationship, axis = [value], "high", "direct urinary analyte", "bisphenol exposure axis"
                source_reason = "CDC NHANES environmental phenols panel and codebook"
                break
        if not candidates:
            confidence, relationship, axis = "unresolved", "bisphenol family searched; no candidate-specific NHANES variable", "bisphenol exposure axis"
            source_reason = "CDC NHANES environmental phenols/parabens panels searched"
    if not candidates and ("pah" in chemical_class or contains_any(raw_name, list(PAH_TO_VARIABLE))):
        for key, value in PAH_TO_VARIABLE.items():
            if key in raw_name:
                candidates, confidence, relationship, axis = [value], "moderate", "parent PAH to validated urinary OH-PAH proxy", "PAH exposure axis"
                source_reason = "CDC NHANES urinary OH-PAH panel; parent-specificity is limited"
                break
        if not candidates and "pah" in chemical_class:
            candidates, confidence, relationship, axis = ["URXP10"], "low", "family-level urinary OH-PAH proxy", "PAH exposure axis"
            source_reason = "CDC NHANES urinary OH-PAH panel; only family-level proxy assigned"
    if not candidates and "heavy_metals" in chemical_class:
        for element, value in METAL_TO_VARIABLE.items():
            if element in raw_name:
                candidates, confidence, relationship, axis = [value], "moderate", "elemental urinary biomarker for parent metal/species", "metal exposure axis"
                source_reason = "CDC NHANES metals blood/urine panels and codebooks"
                break
        if not candidates:
            confidence, relationship, axis = "unresolved", "metal family searched; no mapped elemental analyte", "metal exposure axis"
            source_reason = "CDC NHANES Metals - Urine and blood metal panels searched"
    if not candidates and "pesticide" in chemical_class:
        confidence, relationship, axis = "unresolved", "candidate-specific pesticide biomarker not established from NHANES domains", "pesticide exposure axis"
        source_reason = "CDC NHANES pesticide, organophosphate, carbamate, and organochlorine panels searched"
    if not candidates and ("flame" in chemical_class or "cps " in raw_name):
        available_fr = sorted(registry.loc[registry["laboratory_title"].str.lower().str.contains("flame retard"), "variable"].unique())
        if available_fr:
            candidates, confidence, relationship, axis = available_fr[:1], "low", "family-level flame-retardant metabolite proxy", "flame-retardant exposure axis"
            source_reason = "CDC NHANES flame-retardant metabolite panels searched; candidate-specific identity not verified"
        else:
            confidence, relationship, axis = "unresolved", "no flame-retardant analyte file with individual survey results", "flame-retardant exposure axis"
            source_reason = "CDC NHANES flame-retardant panels searched"
    if not candidates and ("voc" in chemical_class or "volatile" in raw_name or "fluorocarbon" in raw_name):
        confidence, relationship, axis = "unresolved", "candidate not matched to NHANES VOC biomarker variable", "VOC exposure axis"
        source_reason = "CDC NHANES VOC blood and urinary metabolite panels searched"

    available = registry.loc[registry["variable"].isin(candidates)].copy() if candidates else registry.iloc[0:0].copy()
    if candidates and available.empty:
        source_reason += "; proposed variable code absent from downloaded CDC XPT registry"
        confidence = "unresolved"
    if not available.empty:
        candidates = sorted(available["variable"].unique().tolist())
    return {
        "ChemicalID": row["ChemicalID"], "exposure_axis": axis, "candidate_to_axis_relationship": relationship,
        "all_candidate_biomarkers": ";".join(candidates), "mapping_confidence": confidence,
        "mapping_source": source_reason, "searched_domain_count": int(len(domain_hits)),
        "searched_domains": ";".join(sorted(domain_hits["laboratory_title"].drop_duplicates().tolist())),
        "mapping_status": "mapped" if not available.empty else ("resolved_no_candidate_specific_analyte" if confidence == "unresolved" else "unresolved_registry_gap"),
    }


def detectability_for_mapping(mapping: dict[str, object], registry: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    codes = [code for code in str(mapping["all_candidate_biomarkers"]).split(";") if code]
    detail = registry.loc[registry["variable"].isin(codes)].copy()
    if detail.empty:
        return pd.DataFrame(), {"selected_primary_biomarker": "", "selection_reason": "no analyte rows"}
    summary = detail.groupby("variable", as_index=False).agg(
        n_cycles_available=("cycle", "nunique"), total_measured=("n_measured", "sum"),
        total_above_lod=("n_above_lod", "sum"), pooled_above_lod_pct=("above_lod_pct", "mean"),
        min_cycle_above_lod_pct=("above_lod_pct", "min"), median_cycle_above_lod_pct=("above_lod_pct", "median"),
        max_cycle_above_lod_pct=("above_lod_pct", "max"),
    )
    summary["pooled_above_lod_pct"] = 100 * summary["total_above_lod"] / summary["total_measured"].replace(0, np.nan)
    # Pre-outcome representative biomarker choice: mapping specificity is
    # already encoded in candidate order; among equal candidates choose
    # detectability, then cycle coverage, then measured N.
    choice = summary.sort_values(["pooled_above_lod_pct", "n_cycles_available", "total_measured"], ascending=False).iloc[0]
    selected = str(choice["variable"])
    detail["selected_primary_biomarker"] = detail["variable"].eq(selected)
    detail["first_cycle"] = detail["cycle"].min()
    detail["last_cycle"] = detail["cycle"].max()
    return detail, {
        "selected_primary_biomarker": selected,
        "selection_reason": "validated/specificity mapping first; detectability, cycle coverage, and measured N used as non-CRC tie-breakers",
        "all_biomarker_summary": summary.to_dict(orient="records"),
    }


def testability_for_candidate(candidate: pd.Series, selected: str, detect: pd.DataFrame, harmonized: pd.DataFrame) -> dict[str, object]:
    if not selected or detect.empty:
        return {
            "ChemicalID": candidate["ChemicalID"], "selected_primary_biomarker": "", "n_cycles_available": 0,
            "cycle_list": "", "first_cycle": "", "last_cycle": "", "analyte_available_n": 0,
            "crc_outcome_available_n": 0, "exposure_plus_outcome_n": 0, "complete_covariate_n": 0,
            "available_crc_cases": 0, "available_controls": 0, "survey_weight_available": False,
            "strata_psu_available": False, "creatinine_available_if_urine": False, "core_covariates_available": False,
            "T_tag": 0, "T_reason": "no analyte mapping",
        }
    work = []
    rows = detect.loc[detect["variable"].eq(selected)].copy()
    for r in rows.itertuples(index=False):
        local = Path(r.local_xpt)
        frame = pd.read_sas(local, format="xport", encoding="latin1")
        cols = ["SEQN", r.variable]
        if r.flag_variable and r.flag_variable in frame.columns:
            cols.append(r.flag_variable)
        if r.weight_variable and r.weight_variable in frame.columns:
            cols.append(r.weight_variable)
        exp = frame[cols].copy()
        exp["cycle"] = r.cycle
        exp["exposure_value"] = pd.to_numeric(exp[r.variable], errors="coerce")
        exp["weight"] = pd.to_numeric(exp[r.weight_variable], errors="coerce") if r.weight_variable else 1.0
        exp["above_lod"] = exp[r.flag_variable].ne(1) if r.flag_variable else exp["exposure_value"].notna()
        work.append(exp[["SEQN", "cycle", "exposure_value", "weight", "above_lod"]])
    exp = pd.concat(work, ignore_index=True)
    merged = exp.merge(harmonized, on=["SEQN", "cycle"], how="inner", suffixes=("", "_h"))
    outcome_available = merged["cancer_known"].eq(True) if "cancer_known" in merged else merged["cancer_outcome_available"].eq(True)
    merged = merged.loc[outcome_available & merged["weight"].gt(0) & merged["exposure_value"].notna()].copy()
    urine = any(detect.loc[detect["variable"].eq(selected), "matrix"].eq("urine"))
    required = ["age", "sex", "race", "bmi", "smoking", "pir", "outcome"]
    if "outcome" not in merged.columns:
        merged["outcome"] = merged["crc_case"].astype(int)
    # The frozen NHANES harmonizer stores urine creatinine as log2-transformed
    # `creatinine_log2`; do not silently downgrade every urine biomarker to
    # T=0 by looking for a non-existent raw column.
    creatinine_col = "creatinine_log2" if "creatinine_log2" in merged.columns else ("creatinine" if "creatinine" in merged.columns else "")
    if urine and creatinine_col:
        required.append(creatinine_col)
    complete = merged.dropna(subset=[c for c in required if c in merged.columns]).copy()
    cases = int(complete["crc_case"].sum()) if "crc_case" in complete else 0
    controls = int(complete["cancer_free"].sum()) if "cancer_free" in complete else 0
    n = len(complete)
    # Cycle coverage C is an analyte/data-infrastructure property, not a
    # complete-case outcome property. Covariate missingness must not erase a
    # cycle in which the biomarker was actually measured.
    analyte_cycles = sorted(rows["cycle"].dropna().astype(str).unique().tolist())
    weight_ok = bool(len(complete) and complete["weight"].gt(0).all())
    design_ok = bool({"psu", "strata"}.issubset(complete.columns) and complete["psu"].notna().all() and complete["strata"].notna().all())
    covar_ok = all(c in complete.columns for c in ["age", "sex", "race", "bmi", "smoking", "pir"])
    creat_ok = (not urine) or (bool(creatinine_col) and complete[creatinine_col].notna().all())
    if n >= 5000 and cases >= 60 and weight_ok and design_ok and covar_ok and creat_ok:
        t_tag, reason = 2, "predefined primary survey model feasible"
    elif n >= 500 and cases >= 20 and weight_ok and design_ok and covar_ok and creat_ok:
        t_tag, reason = 1, "survey model feasible but event/sample margin is limited"
    else:
        t_tag, reason = 0, "insufficient complete-case events/sample or survey infrastructure"
    return {
        "ChemicalID": candidate["ChemicalID"], "selected_primary_biomarker": selected,
        "n_cycles_available": len(analyte_cycles), "cycle_list": ";".join(analyte_cycles),
        "first_cycle": analyte_cycles[0] if analyte_cycles else "", "last_cycle": analyte_cycles[-1] if analyte_cycles else "",
        "biological_matrix": ";".join(sorted(detect.loc[detect["variable"].eq(selected), "matrix"].dropna().unique().tolist())),
        "nhanes_lab_component": ";".join(sorted(detect.loc[detect["variable"].eq(selected), "laboratory_title"].dropna().unique().tolist())),
        "analyte_available_n": int(exp["exposure_value"].notna().sum()),
        "crc_outcome_available_n": int(len(harmonized)), "exposure_plus_outcome_n": int(len(merged)),
        "complete_covariate_n": n, "available_crc_cases": cases, "available_controls": controls,
        "survey_weight_available": weight_ok, "strata_psu_available": design_ok,
        "creatinine_available_if_urine": creat_ok, "core_covariates_available": covar_ok,
        "T_tag": t_tag, "T_reason": reason,
    }


def main() -> None:
    phase1 = pd.read_csv(PHASE1_MATRIX, low_memory=False)
    if len(phase1) != 267 or phase1["ChemicalID"].nunique() != 267:
        raise ValueError("Phase 1 matrix is not the required 267-row universe")
    # The previous interim matrix contains human-audit columns for only MCOP
    # and MiNP. Remove those stale fields before merging the full v2 audit so
    # that no old value can shadow a newly calculated X/B/D/C/T field.
    stale_human_fields = [
        "exposure_axis_X_tag", "exposure_axis_name", "biomarker_relationship",
        "biomarker_B_tag", "human_biomarker", "biomarker_matrix", "direct_analyte_identity",
        "identity_separation_note", "detectability_D_tag", "above_lod_pct",
        "detectability_reason", "cycle_coverage_C_tag", "n_cycles_available",
        "cycles_available", "testability_T_tag", "testability_basis",
        "analytic_exposure_outcome_n", "available_crc_cases_for_planning",
        "available_cancer_free_controls_for_planning", "novelty_N_tag",
        "manual_review_required", "manual_review_reason", "eligible_permissive",
        "eligible_moderate", "eligible_strict", "priority_tier", "disposition",
        "PRIORITIZATION_OUTCOME_BLINDED",
    ]
    phase1 = phase1.drop(columns=stale_human_fields, errors="ignore")
    catalog = pd.read_csv(CATALOG)
    catalog = catalog.loc[catalog["data_url"].notna() & catalog["data_url"].ne("")].copy()
    registry, _ = build_registry(catalog)
    harmonized = load_harmonized()

    identity = pd.DataFrame([identity_row(row) for _, row in phase1.iterrows()])
    mappings = []
    detects = []
    tests = []
    map_records = {}
    for _, candidate in phase1.iterrows():
        mapping = mapping_for_candidate(candidate, registry)
        detail, choice = detectability_for_mapping(mapping, registry)
        mapping.update({"selected_primary_biomarker": choice["selected_primary_biomarker"], "selection_reason": choice["selection_reason"]})
        if choice.get("all_biomarker_summary"):
            mapping["all_biomarker_summary"] = json.dumps(choice["all_biomarker_summary"], ensure_ascii=False)
        mappings.append(mapping)
        map_records[candidate["ChemicalID"]] = mapping
        if not detail.empty:
            detail["ChemicalID"] = candidate["ChemicalID"]
            detail["ChemicalName"] = candidate["ChemicalName"]
            detects.append(detail)
        tests.append(testability_for_candidate(candidate, choice["selected_primary_biomarker"], detail, harmonized))

    mapping_df = pd.DataFrame(mappings)
    detect_df = pd.concat(detects, ignore_index=True) if detects else pd.DataFrame()
    test_df = pd.DataFrame(tests)
    if not detect_df.empty:
        detect_df = detect_df.rename(columns={"n_measured": "n_measured", "n_above_lod": "n_above_lod"})

    # Detectability eligibility is evaluated on the pre-outcome selected
    # primary biomarker for each chemical. Do not let a low-detectability
    # secondary biomarker (e.g. MiNP) drag down a distinct primary analyte
    # (e.g. MCOP/URXCOP) while the full multi-biomarker evidence remains in
    # OUT_DETECT and all_biomarker_summary.
    selected_lookup = mapping_df[["ChemicalID", "selected_primary_biomarker"]].rename(columns={"selected_primary_biomarker": "selected_primary_biomarker_expected"})
    detect_primary = detect_df.merge(selected_lookup, on="ChemicalID", how="left", validate="many_to_one")
    detect_primary = detect_primary.loc[detect_primary["variable"].eq(detect_primary["selected_primary_biomarker_expected"])].copy()
    summary = detect_primary.groupby("ChemicalID", as_index=False).agg(
        above_lod_pct=("above_lod_pct", "mean"), min_cycle_above_lod_pct=("above_lod_pct", "min"),
        median_cycle_above_lod_pct=("above_lod_pct", "median"), max_cycle_above_lod_pct=("above_lod_pct", "max"),
        n_measured=("n_measured", "sum"), n_above_lod=("n_above_lod", "sum"),
    ) if not detect_primary.empty else pd.DataFrame(columns=["ChemicalID"])
    if not summary.empty:
        summary["above_lod_pct"] = 100 * summary["n_above_lod"] / summary["n_measured"].replace(0, np.nan)
        summary["D_tag"] = np.select([summary["above_lod_pct"] < 50, summary["above_lod_pct"] < 90], [0, 1], default=2)
    else:
        summary = pd.DataFrame(columns=["ChemicalID", "above_lod_pct", "D_tag"])

    audit = phase1.merge(identity, on=["ChemicalID", "ChemicalName"], how="left", validate="one_to_one")
    # Keep the mapping decision authoritative.  The testability table also
    # contains this field for convenience, so rename it before merging to
    # prevent pandas' _x/_y collision from silently blanking the final matrix.
    mapping_for_merge = mapping_df.rename(columns={"selected_primary_biomarker": "mapping_selected_primary_biomarker"})
    test_for_merge = test_df.rename(columns={"selected_primary_biomarker": "test_selected_primary_biomarker"})
    audit = audit.merge(mapping_for_merge, on="ChemicalID", how="left", validate="one_to_one")
    audit = audit.merge(summary, on="ChemicalID", how="left", validate="one_to_one")
    audit = audit.merge(test_for_merge, on="ChemicalID", how="left", validate="one_to_one")
    audit["selected_primary_biomarker"] = audit["mapping_selected_primary_biomarker"].fillna(audit["test_selected_primary_biomarker"])
    audit = audit.drop(columns=["mapping_selected_primary_biomarker", "test_selected_primary_biomarker"], errors="ignore")
    # Restore one canonical identity column where the phase-1 matrix and the
    # identity table carry the same identifier under pandas suffixes.
    for base in ["chemical_class", "CasRN", "PubChemCID", "DTXSID"]:
        left, right = f"{base}_x", f"{base}_y"
        if left in audit.columns and right in audit.columns:
            audit[base] = audit[left].where(audit[left].notna() & audit[left].ne(""), audit[right])
            audit = audit.drop(columns=[left, right])
        elif left in audit.columns:
            audit = audit.rename(columns={left: base})
        elif right in audit.columns:
            audit = audit.rename(columns={right: base})
    audit["X_tag"] = audit["mapping_confidence"].isin(["high", "moderate", "low"]).astype(int)
    audit["B_tag"] = audit["mapping_status"].eq("mapped").astype(int)
    audit["human_biomarker"] = audit["selected_primary_biomarker"]
    audit["nhanes_variable"] = audit["selected_primary_biomarker"]
    audit["biomarker_type"] = np.select(
        [audit["biological_matrix"].astype(str).str.contains("urine", na=False), audit["biological_matrix"].astype(str).str.contains("serum|blood", case=False, regex=True, na=False)],
        ["urinary biomarker", "serum/blood biomarker"],
        default="",
    )
    audit["nhanes_cycles"] = audit["cycle_list"]
    audit["n_cycles"] = audit["n_cycles_available"]
    audit["C_tag"] = np.select([audit["n_cycles_available"].fillna(0) <= 2, audit["n_cycles_available"].fillna(0) < 5], [0, 1], default=2)
    audit["E_tag"] = audit["entity_E_tag"]
    audit["N_tag"] = "pending_manual_review"
    audit["M_tag"] = audit["molecular_level"]
    audit["eligible_permissive"] = (audit["E_tag"].eq(1) & audit["X_tag"].eq(1) & audit["B_tag"].eq(1) & audit["D_tag"].fillna(0).ge(1) & audit["C_tag"].ge(1) & audit["T_tag"].fillna(0).ge(1))
    audit["eligible_moderate"] = (audit["E_tag"].eq(1) & audit["X_tag"].eq(1) & audit["B_tag"].eq(1) & audit["D_tag"].fillna(0).ge(1) & audit["C_tag"].ge(2) & audit["T_tag"].fillna(0).ge(1))
    audit["eligible_strict"] = (audit["E_tag"].eq(1) & audit["X_tag"].eq(1) & audit["B_tag"].eq(1) & audit["D_tag"].fillna(0).eq(2) & audit["C_tag"].eq(2) & audit["T_tag"].fillna(0).eq(2))
    audit["audit_status"] = np.where(audit["mapping_status"].eq("mapped"), "resolved_mapped", "resolved_no_supported_candidate_specific_analyte")
    audit["audit_evidence"] = audit.apply(lambda r: f"CDC NHANES catalog domains checked: {r['searched_domains'] or 'none matched'}; {r['mapping_source']}; local XPT registry variables checked for candidate aliases.", axis=1)
    audit["what_was_checked"] = audit.apply(lambda r: f"CTD identifiers/name; CDC NHANES laboratory catalog; {r['searched_domain_count']} relevant laboratory domains; downloaded XPT value/LOD/weight columns.", axis=1)
    audit["databases_checked"] = "CTD Phase 1 export; CDC NHANES continuous laboratory catalog; CDC NHANES environmental XPT files"
    audit["possible_proxy"] = audit["candidate_to_axis_relationship"].replace("none", "")
    audit["reason_unresolved"] = np.where(audit["mapping_status"].eq("mapped"), "", audit["mapping_source"])
    audit["manual_review_required"] = audit["entity_manual_review_required"].eq(True) | audit["mapping_status"].eq("unresolved_registry_gap")
    audit["manual_review_reason"] = np.select(
        [audit["entity_manual_review_required"].eq(True), audit["mapping_status"].eq("unresolved_registry_gap")],
        [audit["entity_E_reason"], audit["mapping_source"]],
        default="",
    )
    audit["priority_tier"] = np.select([
        audit["eligible_strict"], audit["eligible_moderate"], audit["eligible_permissive"], audit["M_tag"].eq("M2")
    ], ["Tier_A_strict", "Tier_A_moderate", "Tier_B_human_testable", "Tier_C_molecular_only"], default="not_human_testable")
    audit["final_disposition"] = np.where(audit["eligible_permissive"], "enter_systematic_human_screen", "retain_with_resolved_audit_no_primary_screen")
    audit["PRIORITIZATION_OUTCOME_BLINDED"] = True

    identity.to_csv(OUT_IDENTITY, index=False)
    mapping_output = mapping_df.merge(
        test_df[["ChemicalID", "biological_matrix", "nhanes_lab_component", "cycle_list", "n_cycles_available"]],
        on="ChemicalID",
        how="left",
        validate="one_to_one",
    )
    mapping_output["nhanes_variable"] = mapping_output["selected_primary_biomarker"]
    mapping_output["biomarker_type"] = np.select(
        [mapping_output["biological_matrix"].astype(str).str.contains("urine", na=False), mapping_output["biological_matrix"].astype(str).str.contains("serum|blood", case=False, regex=True, na=False)],
        ["urinary biomarker", "serum/blood biomarker"],
        default="",
    )
    mapping_output.to_csv(OUT_MAPPING, index=False)
    detect_df.to_csv(OUT_DETECT, index=False)
    test_df.to_csv(OUT_TEST, index=False)
    audit.to_csv(OUT_MATRIX, index=False)
    # The queue is reserved for genuine identity/registry exceptions, not
    # every candidate that fails an eligibility gate. Resolved negatives stay
    # in the full 267-row matrix with their evidence trail.
    audit.loc[audit["manual_review_required"]].to_csv(OUT_REVIEW, index=False)

    # Sequential flow counts are the actual attrition path, rather than
    # marginal counts that can be misread as a funnel.
    flow_rows = []
    flow_specs = [
        ("total_core_chemicals", pd.Series(True, index=audit.index), "all Phase 1 primary chemicals"),
        ("E_entity_valid", audit["E_tag"].eq(1), "E_tag == 1"),
        ("E_and_X_interpretable_exposure", audit["E_tag"].eq(1) & audit["X_tag"].eq(1), "E_tag == 1 and X_tag == 1"),
        ("E_X_and_B_biomarker_available", audit["E_tag"].eq(1) & audit["X_tag"].eq(1) & audit["B_tag"].eq(1), "E/X/B gates"),
        ("E_X_B_and_D_detectable", audit["E_tag"].eq(1) & audit["X_tag"].eq(1) & audit["B_tag"].eq(1) & audit["D_tag"].fillna(0).ge(1), "E/X/B and D >= 1"),
        ("E_X_B_D_and_C_coverage", audit["E_tag"].eq(1) & audit["X_tag"].eq(1) & audit["B_tag"].eq(1) & audit["D_tag"].fillna(0).ge(1) & audit["C_tag"].fillna(0).ge(1), "E/X/B/D and C >= 1"),
        ("E_X_B_D_C_and_T_testable", audit["eligible_permissive"].eq(True), "E=1 & X=1 & B=1 & D>=1 & C>=1 & T>=1"),
    ]
    for label, mask, definition in flow_specs:
        flow_rows.append({"stage": label, "n": int(mask.sum()), "definition": definition})
    flow_rows.append({"stage": "moderate_eligibility", "n": int(audit["eligible_moderate"].sum()), "definition": "E=1 & X=1 & B=1 & D>=1 & C>=2 & T>=1"})
    flow_rows.append({"stage": "strict_eligibility", "n": int(audit["eligible_strict"].sum()), "definition": "E=1 & X=1 & B=1 & D=2 & C=2 & T=2"})
    pd.DataFrame(flow_rows).to_csv(OUT_FLOW, index=False)

    forbidden = ["crc_or", "crc_ci", "crc_p", "crc_fdr", "loco_effect", "cycle_specific_effect", "human_or", "human_p"]
    eligibility_features = ["E_tag", "X_tag", "B_tag", "D_tag", "C_tag", "T_tag"]
    source_columns = set(audit.columns)
    forbidden_present = sorted(c for c in source_columns if any(f in c.lower() for f in forbidden))
    firewall = {
        "PRIORITIZATION_OUTCOME_BLINDED": True,
        "eligibility_features": eligibility_features,
        "forbidden_human_outcome_fields": forbidden,
        "forbidden_columns_present_in_actionability_matrix": forbidden_present,
        "candidate_specific_crc_effect_used": False,
        "note": "CRC case counts are used only as pre-outcome testability infrastructure; no human effect estimate enters eligibility.",
    }
    OUT_FIREWALL.write_text(json.dumps(firewall, indent=2, ensure_ascii=False), encoding="utf-8")

    counts = {"n_total": 267, "n_X": int(audit.X_tag.sum()), "n_B": int(audit.B_tag.sum()), "n_D": int(audit.D_tag.fillna(0).ge(1).sum()), "n_C": int(audit.C_tag.ge(1).sum()), "n_T": int(audit.T_tag.fillna(0).ge(1).sum()), "n_eligible_permissive": int(audit.eligible_permissive.sum()), "n_eligible_moderate": int(audit.eligible_moderate.sum()), "n_eligible_strict": int(audit.eligible_strict.sum()), "n_resolved": int((audit.audit_status != "unresolved").sum()), "n_unresolved": int((audit.audit_status == "unresolved").sum())}
    manifest = {"phase": "Complete 267-chemical Human Actionability Audit v2", "generated_at_utc": datetime.now(timezone.utc).isoformat(), "counts": counts, "environmental_xpt_file_count": int(len(registry[["cycle", "data_file"]].drop_duplicates())), "outputs": [str(p) for p in [OUT_IDENTITY, OUT_MAPPING, OUT_DETECT, OUT_TEST, OUT_MATRIX, OUT_FLOW, OUT_REVIEW, OUT_FIREWALL]], "outcome_blinded": True, "rule": "E=1 & X=1 & B=1 & D>=1 & C>=1 & T>=1", "note": "All 267 candidates receive an evidence-backed mapping decision; unresolved is reserved for genuine identity/mapping ambiguity, not absence from the generic queue."}
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(counts, ensure_ascii=False))


if __name__ == "__main__":
    main()
