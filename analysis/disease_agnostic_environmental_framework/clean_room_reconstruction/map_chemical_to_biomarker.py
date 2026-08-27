"""Clean-room chemical-to-analyte mapping and actionability gates.

The rules are technical exposure-infrastructure rules only.  No disease or
outcome field is accepted by this module.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd


METAL_TO_ANALYTE = {
    "barium": "URXUBA", "cadmium": "URXUCD", "cobalt": "URXUCO", "cesium": "URXUCS",
    "molybdenum": "URXUMO", "manganese": "URXUMN", "lead": "URXUPB", "platinum": "URXUPT",
    "antimony": "URXUSB", "tin": "URXUSN", "silver": "URXUSR", "thallium": "URXUTL",
    "tungsten": "URXUTU", "uranium": "URXUUR", "mercury": "URXUHG",
}
PHENOL_TO_ANALYTE = {
    "bisphenol a": "URXBPH", "bisphenol": "URXBPH", "benzophenone-3": "URXBP3",
    "butylparaben": "URXBUP", "ethylparaben": "URXEPB", "methylparaben": "URXMPB",
    "propylparaben": "URXPPB",
}
PAH_TO_ANALYTE = {
    "naphthalene": "URXP02", "methylnaphthalene": "URXP02", "dimethylnaphthalene": "URXP02",
    "fluorene": "URXP04", "phenanthrene": "URXP25", "methylphenanthrene": "URXP25",
    "pyrene": "URXP10", "benzo": "URXP10", "chrysene": "URXP10", "anthracene": "URXP25",
    "fluoranthene": "URXP10", "acenaphth": "URXP02", "indeno": "URXP10", "perylene": "URXP10",
    "picene": "URXP10", "retene": "URXP10", "corannulene": "URXP10",
}
PHTHALATE_TO_ANALYTE = {
    "butylbenzyl phthalate": ["URXMZP"], "mono-benzyl phthalate": ["URXMZP"],
    "monobutyl phthalate": ["URXMBP"], "monoethyl phthalate": ["URXMEP"],
    "mono-(2-ethylhexyl)phthalate": ["URXMHP"], "mono(2-ethyl-5-hydroxyhexyl)phthalate": ["URXMHH"],
    "mono(2-ethyl-5-oxohexyl)phthalate": ["URXMOH"], "2-ethyl-5-carboxypentyl phthalate": ["URXECP"],
    "mono-isobutyl phthalate": ["URXMIB"], "monoisononylphthalate": ["URXMNP"],
    "mono(carboxy-isooctyl)phthalate": ["URXCOP"], "diethylhexyl phthalate": ["URXMHH", "URXMOH", "URXECP"],
    "dinonylphthalate": ["URXMNP", "URXCOP"], "diisononyl phthalate": ["URXMNP", "URXCOP"],
}
PFAS_TO_ANALYTE = {
    "perfluorooctane sulfonic acid": ["LBXNFOS"], "perfluorobutanesulfonic acid": ["LBXPFBS"],
    "perfluorohexanesulfonic acid": ["LBXPFHS"], "perfluorohexanoic acid": [],
    "perfluorodecanoic acid": ["LBXPFDE"], "perfluorononanoic acid": ["LBXPFNA"],
    "perfluoro-n-nonanoic acid": ["LBXPFNA"], "perfluoroundecanoic acid": ["LBXPFUA"],
    "perfluorododecanoic acid": ["LBXPFDO"], "perfluoroheptanoic acid": ["LBXPFHP"],
    "perfluoro-n-heptanoic acid": ["LBXPFHP"], "2-(n-methyl-pfosa) acetate": ["LBXMPAH"],
}


def norm(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value or "").lower()).strip()


def contains_any(text: str, terms: list[str]) -> bool:
    return any(term in text for term in terms)


def map_candidate(row: pd.Series, registry: pd.DataFrame) -> dict[str, object]:
    name = norm(row["canonical_name"])
    raw = str(row["canonical_name"]).lower()
    chemical_class = str(row.get("chemical_class", "")).lower()
    candidates: list[str] = []
    confidence, relationship, axis, reason = "unresolved", "none", "", ""

    genx_like = any(token in raw for token in ("4,8-dioxa-3h-perfluorononanoic", "hexafluoropropylene oxide dimer acid", "hfpo-da", "hfpo da", "genx"))
    for key, values in PHTHALATE_TO_ANALYTE.items():
        if key in raw:
            candidates = values
            confidence = "high" if raw.startswith("mono") else "moderate"
            relationship = "direct urinary metabolite" if raw.startswith("mono") else "parent-to-validated-urinary-metabolite axis"
            axis = "DINP-related exposure axis" if any(x in raw for x in ("isononyl", "carboxy-isooctyl", "dinonyl")) else ("DEHP-related exposure axis" if "ethylhexyl" in raw else "phthalate exposure axis")
            reason = "CDC NHANES phthalate/plasticizer metabolite panel and codebook"
            break

    if not candidates and ("perfluoro" in raw or "pfas" in chemical_class):
        if genx_like:
            relationship, axis = "PFAS class searched; no candidate-specific NHANES GenX/HFPO-DA variable", "PFAS exposure axis"
            reason = "CDC NHANES PFAS panels searched; PFNA explicitly excluded as a non-equivalent analyte"
        else:
            for key, values in PFAS_TO_ANALYTE.items():
                if key in raw:
                    candidates, confidence, relationship, axis = values, "high", "direct serum PFAS analyte", "PFAS exposure axis"
                    reason = "CDC NHANES PFAS serum panel and codebook"
                    break
        if not candidates and not reason and "fluorotelomer" not in raw:
            relationship, axis = "PFAS class searched; no candidate-specific NHANES variable", "PFAS exposure axis"
            reason = "CDC NHANES PFAS panels searched; no specific analyte name match"

    if not candidates and ("bisphenol" in raw or "bisphenol" in chemical_class):
        for key, value in PHENOL_TO_ANALYTE.items():
            if key in raw:
                candidates, confidence, relationship, axis = [value], "high", "direct urinary analyte", "bisphenol exposure axis"
                reason = "CDC NHANES environmental phenols panel and codebook"
                break
        if not candidates:
            relationship, axis = "bisphenol family searched; no candidate-specific NHANES variable", "bisphenol exposure axis"
            reason = "CDC NHANES environmental phenols/parabens panels searched"

    if not candidates and ("pah" in chemical_class or contains_any(raw, list(PAH_TO_ANALYTE))):
        for key, value in PAH_TO_ANALYTE.items():
            if key in raw:
                candidates, confidence, relationship, axis = [value], "moderate", "parent PAH to validated urinary OH-PAH proxy", "PAH exposure axis"
                reason = "CDC NHANES urinary OH-PAH panel; parent-specificity is limited"
                break
        if not candidates and "pah" in chemical_class:
            candidates, confidence, relationship, axis = ["URXP10"], "low", "family-level urinary OH-PAH proxy", "PAH exposure axis"
            reason = "CDC NHANES urinary OH-PAH panel; only family-level proxy assigned"

    if not candidates and "heavy_metals" in chemical_class:
        for element, value in METAL_TO_ANALYTE.items():
            if element in raw:
                candidates, confidence, relationship, axis = [value], "moderate", "elemental urinary biomarker for parent metal/species", "metal exposure axis"
                reason = "CDC NHANES metals blood/urine panels and codebooks"
                break
        if not candidates:
            relationship, axis = "metal family searched; no mapped elemental analyte", "metal exposure axis"
            reason = "CDC NHANES Metals - Urine and blood metal panels searched"

    if not candidates and "pesticide" in chemical_class:
        relationship, axis = "candidate-specific pesticide biomarker not established from NHANES domains", "pesticide exposure axis"
        reason = "CDC NHANES pesticide, organophosphate, carbamate, and organochlorine panels searched"
    if not candidates and ("flame" in chemical_class or "cps " in raw):
        available = sorted(registry.loc[registry["laboratory_title"].str.lower().str.contains("flame retard"), "variable"].unique())
        if available:
            candidates, confidence, relationship, axis = available[:1], "low", "family-level flame-retardant metabolite proxy", "flame-retardant exposure axis"
            reason = "CDC NHANES flame-retardant metabolite panels searched; candidate-specific identity not verified"
        else:
            relationship, axis = "no flame-retardant analyte file with individual survey results", "flame-retardant exposure axis"
            reason = "CDC NHANES flame-retardant panels searched"
    if not candidates and ("voc" in chemical_class or "volatile" in raw or "fluorocarbon" in raw):
        relationship, axis = "candidate not matched to NHANES VOC biomarker variable", "VOC exposure axis"
        reason = "CDC NHANES VOC blood and urinary metabolite panels searched"

    available = registry.loc[registry["variable"].isin(candidates)].copy() if candidates else registry.iloc[0:0].copy()
    if candidates and available.empty:
        reason += "; proposed variable code absent from downloaded CDC XPT registry"
        confidence = "unresolved"
    if not available.empty:
        candidates = sorted(available["variable"].unique().tolist())
    return {
        "chemical_id": row["chemical_id"], "all_candidate_biomarkers": ";".join(candidates),
        "exposure_axis": axis, "mapping_type": relationship, "mapping_confidence": confidence,
        "mapping_source": reason, "mapping_status": "mapped" if not available.empty else "resolved_no_candidate_specific_analyte",
    }


def variable_summary(registry: pd.DataFrame) -> pd.DataFrame:
    if registry.empty:
        return pd.DataFrame()
    out = registry.groupby("variable", as_index=False).agg(
        n_cycles=("cycle", "nunique"), n_measured=("n_measured", "sum"), n_above_lod=("n_above_lod", "sum"),
        min_cycle_above_lod_pct=("above_lod_pct", "min"), median_cycle_above_lod_pct=("above_lod_pct", "median"),
        max_cycle_above_lod_pct=("above_lod_pct", "max"), cycle_list=("cycle", lambda x: ";".join(sorted(set(map(str, x))))),
        weight_variable=("weight_variable", lambda x: ";".join(sorted(set(map(str, x))))),
        matrix=("matrix", lambda x: ";".join(sorted(set(map(str, x))))), laboratory_title=("laboratory_title", lambda x: ";".join(sorted(set(map(str, x))))),
    )
    out["pooled_above_lod_pct"] = 100 * out["n_above_lod"] / out["n_measured"].replace(0, np.nan)
    out["D_tag"] = np.select([out["pooled_above_lod_pct"] < 50, out["pooled_above_lod_pct"] < 90], [0, 1], default=2)
    out["C_tag"] = np.select([out["n_cycles"] <= 2, out["n_cycles"] < 5], [0, 1], default=2)
    out["F_tag"] = out["weight_variable"].ne("").astype(int)
    return out
