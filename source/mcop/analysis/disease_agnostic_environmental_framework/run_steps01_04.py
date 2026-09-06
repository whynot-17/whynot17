"""Build and freeze the disease-agnostic environmental test space (Steps 1–4).

This runner deliberately uses only CTD chemical vocabulary/classification and
NHANES laboratory infrastructure. It does not import the CRC harmonized frame,
GeneCards files, CTD chemical–gene interactions, or any disease outcome field.
The existing audit module is imported only for its outcome-free NHANES registry
and deterministic chemical-to-analyte mapping functions.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "analysis" / "disease_agnostic_environmental_framework"
DEFAULT_CTD = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data" / "CTD_chemicals.tsv.gz"
DEFAULT_RULES = ROOT / "work" / "environmental_toxicology_crc_phase1" / "chemical_class_rules.json"
DEFAULT_DRUGCENTRAL = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data" / "drugcentral_structures.smiles.tsv"
DEFAULT_PAH_FORMULAS = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data" / "pubchem_pah_formula.json"
DEFAULT_CATALOG = ROOT / "outputs" / "environmental_crc_267_nhanes_environmental_lab_catalog.csv"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def as_bool(value: object) -> bool:
    return bool(value) if not pd.isna(value) else False


def classify_identity(row: pd.Series, chemical_class: str) -> dict[str, str]:
    name = clean(row.get("ChemicalName"))
    synonyms = ";".join(x for x in [clean(row.get("MESHSynonyms")), clean(row.get("CTDCuratedSynonyms"))] if x)
    status = "environmental chemical entity"
    low = name.lower()
    if "mono" in low or "metabolite" in low or "hydroxy" in low or "oxo" in low or "carboxy" in low:
        status = "possible metabolite or transformation product; parent relationship not inferred"
    elif chemical_class:
        status = "parent/classified environmental chemical entity"
    return {"synonym": synonyms, "parent_metabolite_status": status}


def build_s1(ctd: pd.DataFrame, classified: pd.DataFrame) -> pd.DataFrame:
    ctd = ctd.copy()
    ctd["ChemicalID"] = ctd["ChemicalID"].astype(str).str.replace(r"^MESH:", "", regex=True)
    frame = ctd.merge(
        classified[
            [
                "ChemicalID", "chemical_class", "classification_rule_prefixes", "is_core",
                "drug_like_exclusion", "drugcentral_match", "classification_exclusion_reason",
            ]
        ],
        on="ChemicalID", how="left", validate="one_to_one",
    )
    rows = []
    for row in frame.to_dict("records"):
        cls = clean(row.get("chemical_class"))
        ident = classify_identity(pd.Series(row), cls)
        eligible = as_bool(row.get("is_core", False))
        exclusion = "" if eligible else (clean(row.get("classification_exclusion_reason")) or "outside prespecified CTD environmental MeSH branches")
        rows.append(
            {
                "chemical_id": clean(row.get("ChemicalID")),
                "canonical_name": clean(row.get("ChemicalName")),
                "synonym": ident["synonym"],
                "CAS": clean(row.get("CasRN")),
                "PubChemCID": clean(row.get("PubChemCID")),
                "DTXSID": clean(row.get("DTXSID")),
                "InChIKey": clean(row.get("InChIKey")),
                "chemical_class": cls,
                "source": "CTD chemical vocabulary; prespecified CTD MeSH hierarchy rules",
                "parent_metabolite_status": ident["parent_metabolite_status"],
                "eligible_as_environmental_exposure": eligible,
                "exclusion_reason": exclusion,
                "classification_rule_prefixes": clean(row.get("classification_rule_prefixes")),
                "drug_like_exclusion": as_bool(row.get("drug_like_exclusion")),
                "drugcentral_match": as_bool(row.get("drugcentral_match")),
            }
        )
    out = pd.DataFrame(rows).drop_duplicates("chemical_id", keep="first")
    if out["chemical_id"].duplicated().any():
        raise AssertionError("S1 contains duplicate chemical IDs")
    return out


def variable_summary(registry: pd.DataFrame) -> pd.DataFrame:
    if registry.empty:
        return pd.DataFrame(columns=["variable", "n_cycles", "cycle_list", "pooled_above_lod_pct", "D_tag", "weight_variable", "matrix", "laboratory_title"])
    out = registry.groupby("variable", as_index=False).agg(
        n_cycles=("cycle", "nunique"),
        cycle_list=("cycle", lambda x: ";".join(sorted(set(map(str, x))))),
        n_measured=("n_measured", "sum"),
        n_above_lod=("n_above_lod", "sum"),
        min_cycle_above_lod_pct=("above_lod_pct", "min"),
        median_cycle_above_lod_pct=("above_lod_pct", "median"),
        max_cycle_above_lod_pct=("above_lod_pct", "max"),
        weight_variable=("weight_variable", lambda x: ";".join(sorted(set(map(str, x))))),
        matrix=("matrix", lambda x: ";".join(sorted(set(map(str, x))))),
        laboratory_title=("laboratory_title", lambda x: ";".join(sorted(set(map(str, x))))),
    )
    out["pooled_above_lod_pct"] = 100.0 * out["n_above_lod"] / out["n_measured"].replace(0, np.nan)
    out["D_tag"] = np.select([out["pooled_above_lod_pct"] < 50, out["pooled_above_lod_pct"] < 90], [0, 1], default=2)
    out["C_tag"] = np.select([out["n_cycles"] <= 2, out["n_cycles"] < 5], [0, 1], default=2)
    out["F_tag"] = out["weight_variable"].ne("").astype(int)
    return out


def build_s2_s3(universe: pd.DataFrame, registry: pd.DataFrame, mapper) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary = variable_summary(registry)
    summary_lookup = summary.set_index("variable", drop=False).to_dict("index") if not summary.empty else {}
    s2_rows: list[dict[str, object]] = []
    s3_rows: list[dict[str, object]] = []
    detect_rows: list[pd.DataFrame] = []

    for candidate in universe.to_dict("records"):
        candidate_series = pd.Series({"ChemicalID": candidate["chemical_id"], "ChemicalName": candidate["canonical_name"], "chemical_class": candidate["chemical_class"], "CasRN": candidate["CAS"], "PubChemCID": candidate["PubChemCID"], "DTXSID": candidate["DTXSID"]})
        mapping = mapper.mapping_for_candidate(candidate_series, registry)
        detail, choice = mapper.detectability_for_mapping(mapping, registry)
        selected = clean(choice.get("selected_primary_biomarker"))
        proposed = [x for x in clean(mapping.get("all_candidate_biomarkers")).split(";") if x]
        if not proposed:
            proposed = [""]
        for variable in proposed:
            info = summary_lookup.get(variable, {})
            mapped = bool(variable and variable in summary_lookup and clean(mapping.get("mapping_status")) == "mapped")
            s2_rows.append(
                {
                    "chemical_id": candidate["chemical_id"],
                    "chemical_name": candidate["canonical_name"],
                    "chemical_class": candidate["chemical_class"],
                    "human_biomarker": variable,
                    "exposure_axis": clean(mapping.get("exposure_axis")),
                    "mapping_type": clean(mapping.get("candidate_to_axis_relationship")),
                    "mapping_confidence": clean(mapping.get("mapping_confidence")),
                    "mapping_source": clean(mapping.get("mapping_source")),
                    "mapping_status": clean(mapping.get("mapping_status")),
                    "candidate_is_primary": bool(variable and variable == selected),
                    "NHANES_variable": variable,
                    "matrix": clean(info.get("matrix")),
                    "NHANES_lab_component": clean(info.get("laboratory_title")),
                    "cycle_list": clean(info.get("cycle_list")),
                    "n_cycles_available": int(info.get("n_cycles", 0) or 0),
                    "weight_variable": clean(info.get("weight_variable")),
                    "n_measured": int(info.get("n_measured", 0) or 0),
                    "n_above_lod": int(info.get("n_above_lod", 0) or 0),
                    "pooled_above_lod_pct": float(info.get("pooled_above_lod_pct", np.nan)) if info else np.nan,
                    "mapped": mapped,
                }
            )
            identity_gate = int(bool(candidate["chemical_id"] and candidate["canonical_name"] and mapped))
            exposure_gate = int(mapped and clean(mapping.get("mapping_confidence")) in {"high", "moderate", "low"})
            nhanes_gate = int(mapped)
            d_tag = int(info.get("D_tag", 0) or 0)
            c_tag = int(info.get("C_tag", 0) or 0)
            f_tag = int(info.get("F_tag", 0) or 0)
            actionable = bool(identity_gate and exposure_gate and nhanes_gate and d_tag >= 1 and c_tag >= 1 and f_tag == 1)
            if not mapped:
                reason = "no candidate-specific NHANES analyte with a usable weight in the downloaded registry"
            elif d_tag == 0:
                reason = "detectability gate failed: pooled above-LOD <50%"
            elif c_tag == 0:
                reason = "cycle-coverage gate failed: <=2 cycles"
            elif f_tag == 0:
                reason = "survey-design infrastructure gate failed: no laboratory weight"
            else:
                reason = "all fixed pre-disease actionability gates passed"
            s3_rows.append(
                {
                    "chemical_id": candidate["chemical_id"],
                    "chemical_name": candidate["canonical_name"],
                    "human_biomarker": variable,
                    "mapping_status": clean(mapping.get("mapping_status")),
                    "identity_gate_A": identity_gate,
                    "exposure_interpretability_gate_B": exposure_gate,
                    "NHANES_availability_gate_C": nhanes_gate,
                    "detectability_D_tag": d_tag,
                    "cycle_coverage_E_tag": c_tag,
                    "survey_design_gate_F": f_tag,
                    "pooled_above_lod_pct": float(info.get("pooled_above_lod_pct", np.nan)) if info else np.nan,
                    "n_cycles_available": int(info.get("n_cycles", 0) or 0),
                    "cycle_list": clean(info.get("cycle_list")),
                    "weight_variable": clean(info.get("weight_variable")),
                    "actionable_mapping": actionable,
                    "exclusion_reason": "" if actionable else reason,
                }
            )
        if not detail.empty:
            detail = detail.copy()
            detail["chemical_id"] = candidate["chemical_id"]
            detail["chemical_name"] = candidate["canonical_name"]
            detect_rows.append(detail)

    s2 = pd.DataFrame(s2_rows)
    s3 = pd.DataFrame(s3_rows)
    detects = pd.concat(detect_rows, ignore_index=True) if detect_rows else pd.DataFrame()
    return s2, s3, detects


def build_s4(s2: pd.DataFrame, s3: pd.DataFrame) -> pd.DataFrame:
    actionable = s3.loc[s3["actionable_mapping"] & s3["human_biomarker"].ne("")].copy()
    if actionable.empty:
        return pd.DataFrame(columns=["test_id", "biomarker", "variable", "matrix", "cycles", "weight", "mapping_count", "chemical_ids", "chemical_names", "exposure_axes"])
    joined = actionable.merge(
        s2[["chemical_id", "human_biomarker", "matrix", "NHANES_lab_component", "exposure_axis", "mapping_confidence"]],
        on=["chemical_id", "human_biomarker"], how="left", validate="one_to_many",
    )
    rows = []
    for variable, group in joined.groupby("human_biomarker", sort=True):
        rows.append(
            {
                "test_id": f"NHANES_{variable}",
                "biomarker": variable,
                "variable": variable,
                "matrix": ";".join(sorted(set(group["matrix"].dropna().astype(str)))),
                "cycles": ";".join(sorted(set(";".join(group["cycle_list"].fillna("")).split(";")) - {""})),
                "weight": ";".join(sorted(set(group["weight_variable"].dropna().astype(str)) - {""})),
                "mapping_count": int(len(group)),
                "chemical_ids": ";".join(sorted(set(group["chemical_id"]))),
                "chemical_names": ";".join(sorted(set(group["chemical_name"]))),
                "exposure_axes": ";".join(sorted(set(group["exposure_axis"].dropna().astype(str)) - {""})),
                "n_cycles": int(group["n_cycles_available"].max()),
                "pooled_above_lod_pct": float(group["pooled_above_lod_pct"].iloc[0]),
            }
        )
    return pd.DataFrame(rows).sort_values("test_id").reset_index(drop=True)


def write_audit(out: Path, s1: pd.DataFrame, s2: pd.DataFrame, s3: pd.DataFrame, s4: pd.DataFrame, source_meta: dict[str, object]) -> None:
    eligible = s1["eligible_as_environmental_exposure"].astype(bool)
    mapped = s2.loc[s2["mapped"] & s2["human_biomarker"].ne("")]
    actionable = s3.loc[s3["actionable_mapping"]]
    lines = [
        "# S1 universe audit",
        "",
        f"- Raw CTD vocabulary rows: **{len(s1):,}**",
        f"- Unique CTD chemical IDs: **{s1['chemical_id'].nunique():,}**",
        f"- Prespecified environmental universe after CTD MeSH classification: **{int(eligible.sum()):,}**",
        f"- Excluded from environmental universe: **{int((~eligible).sum()):,}**",
        "- Exclusion was determined by the frozen CTD MeSH hierarchy and deterministic drug-exclusion rules; no disease field was read.",
        "",
        "## Classification exclusions",
        "",
        "| Reason | N |",
        "| --- | ---: |",
    ]
    reasons = s1.loc[~eligible, "exclusion_reason"].value_counts(dropna=False)
    for reason, n in reasons.items():
        lines.append(f"| {str(reason).replace('|', '\\|')} | {int(n):,} |")
    lines += [
        "",
        "# PRE-DISEASE audit",
        "",
        f"- Environmental chemical universe N: **{int(eligible.sum()):,}**",
        f"- Chemical–biomarker mappings with a registry analyte: **{int(len(mapped)):,}**",
        f"- Actionable chemical–biomarker mappings: **{int(len(actionable)):,}**",
        f"- Unique NHANES biomarker tests: **{int(len(s4)):,}**",
        "- The mapping ledger was not collapsed before actionability assessment.",
        "- Step 4 collapse was performed only at the unique NHANES test level.",
        "- No CRC outcome, CRC association, GeneCards disease set, CTD chemical–gene interaction, or transcriptomic result was loaded.",
        "",
        "## Frozen unique test list",
        "",
        "| Test | Biomarker | Matrix | Cycles | Mapping count |",
        "| --- | --- | --- | --- | ---: |",
    ]
    for row in s4.itertuples(index=False):
        lines.append(f"| {row.test_id} | {row.biomarker} | {row.matrix} | {row.cycles} | {row.mapping_count} |")
    lines += [
        "",
        "## Source manifest",
        "",
        f"- CTD source: `{source_meta['ctd_path']}`",
        f"- CTD SHA-256: `{source_meta['ctd_sha256']}`",
        f"- NHANES catalog SHA-256: `{source_meta['catalog_sha256']}`",
        f"- Runner SHA-256: `{source_meta['runner_sha256']}`",
        "",
        "`DISEASE INFORMATION USED IN STEPS 1–4: NO`",
    ]
    (out / "audits" / "S1_universe_audit.md").write_text("\n".join(lines), encoding="utf-8")
    (out / "audits" / "PRE_DISEASE_AUDIT.md").write_text("\n".join(lines[lines.index("# PRE-DISEASE audit"):]), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ctd", type=Path, default=DEFAULT_CTD)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--drugcentral", type=Path, default=DEFAULT_DRUGCENTRAL)
    parser.add_argument("--pah-formulas", type=Path, default=DEFAULT_PAH_FORMULAS)
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for required in [args.ctd, args.rules, args.catalog]:
        if not required.exists():
            raise FileNotFoundError(required)

    phase1 = load_module(ROOT / "work" / "scripts" / "environmental_toxicology_crc_phase1.py", "ctd_outcome_free_classifier")
    mapper = load_module(ROOT / "work" / "scripts" / "environmental_crc_267_biomarker_audit_v2.py", "nhanes_outcome_free_mapper")
    actionability_rules_path = args.outdir / "step01_environmental_universe" / "actionability_rules.json"
    if not actionability_rules_path.exists():
        raise FileNotFoundError(actionability_rules_path)

    ctd = phase1.read_ctd_tsv(args.ctd)
    ctd["ChemicalID"] = ctd["ChemicalID"].astype(str).str.replace(r"^MESH:", "", regex=True)
    classified = phase1.classify_chemicals(ctd, args.rules, args.drugcentral if args.drugcentral.exists() else None, args.pah_formulas if args.pah_formulas.exists() else None)
    s1 = build_s1(ctd, classified)
    env = s1.loc[s1["eligible_as_environmental_exposure"]].copy()

    catalog = pd.read_csv(args.catalog, dtype=str, keep_default_na=False)
    catalog = catalog.loc[catalog["data_url"].notna() & catalog["data_url"].ne("")].copy()
    registry, _ = mapper.build_registry(catalog)
    s2, s3, detect = build_s2_s3(env, registry, mapper)
    s4 = build_s4(s2, s3)

    out1 = args.outdir / "step01_environmental_universe" / "environmental_universe.csv"
    out2 = args.outdir / "step02_biomarker_mapping" / "chemical_biomarker_mapping.csv"
    out3 = args.outdir / "step03_actionability" / "actionability_ledger.csv"
    out3x = args.outdir / "step03_actionability" / "exclusion_ledger.csv"
    out4 = args.outdir / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
    out_detect = args.outdir / "data_processed" / "detectability_registry_outcome_blinded.csv"
    for path in [out1, out2, out3, out3x, out4, out_detect]:
        path.parent.mkdir(parents=True, exist_ok=True)
    # Keep the declared S1 universe as the 2,042 eligible entities. The full
    # CTD classification ledger is retained separately for exclusion audit.
    env.to_csv(out1, index=False)
    s1.to_csv(out1.parent / "ctd_classification_ledger.csv", index=False)
    s2.to_csv(out2, index=False)
    s3.to_csv(out3, index=False)
    s3.loc[~s3["actionable_mapping"]].to_csv(out3x, index=False)
    s4.to_csv(out4, index=False)
    detect.to_csv(out_detect, index=False)

    rules_meta = json.loads(args.rules.read_text(encoding="utf-8"))
    actionability_rules = json.loads(actionability_rules_path.read_text(encoding="utf-8"))
    source_meta = {
        "ctd_path": str(args.ctd),
        "ctd_sha256": sha256(args.ctd),
        "rules_path": str(args.rules),
        "rules_sha256": sha256(args.rules),
        "actionability_rules_path": str(actionability_rules_path),
        "actionability_rules_sha256": sha256(actionability_rules_path),
        "catalog_path": str(args.catalog),
        "catalog_sha256": sha256(args.catalog),
        "runner_path": str(Path(__file__)),
        "runner_sha256": sha256(Path(__file__)),
        "drugcentral_path": str(args.drugcentral) if args.drugcentral.exists() else None,
        "drugcentral_sha256": sha256(args.drugcentral) if args.drugcentral.exists() else None,
        "pah_formulas_path": str(args.pah_formulas) if args.pah_formulas.exists() else None,
        "pah_formulas_sha256": sha256(args.pah_formulas) if args.pah_formulas.exists() else None,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "disease_information_used": False,
        "disease_fields_loaded": [],
        "gene_disease_information_used": False,
        "environmental_universe_n": int(env["chemical_id"].nunique()),
        "raw_ctd_rows_n": int(len(s1)),
        "actionable_mapping_n": int(len(s3.loc[s3["actionable_mapping"]])),
        "unique_test_n": int(len(s4)),
        "planned_downstream_fdr_denominator": int(len(s4)),
        "actionability_rules": rules_meta,
        "pre_disease_gate_rules": actionability_rules,
    }
    lock = {
        "lock_type": "PRE_DISEASE_TESTSET_LOCK",
        "lock_timestamp_utc": source_meta["generated_utc"],
        "git_commit": "resolved at commit time",
        "input_hashes": source_meta,
        "classification_rules": rules_meta,
        "actionability_rules": actionability_rules,
        "included_chemical_biomarker_mappings": int(len(s3.loc[s3["actionable_mapping"]])),
        "excluded_chemical_biomarker_mappings": int(len(s3.loc[~s3["actionable_mapping"]])),
        "unique_test_count": int(len(s4)),
        "biomarkers": s4["biomarker"].tolist(),
        "planned_fdr_denominator": int(len(s4)),
        "disease_information_used": False,
        "gene_disease_information_used": False,
        "historical_results_used_for_eligibility": False,
    }
    (args.outdir / "step04_testset_freeze" / "PRE_DISEASE_TESTSET_LOCK.json").write_text(json.dumps(lock, ensure_ascii=False, indent=2), encoding="utf-8")
    (args.outdir / "data_processed" / "run_manifest.json").write_text(json.dumps(source_meta, ensure_ascii=False, indent=2), encoding="utf-8")
    write_audit(args.outdir, s1, s2, s3, s4, source_meta)

    print(json.dumps({"environmental_chemical_universe_n": source_meta["environmental_universe_n"], "actionable_mapping_n": source_meta["actionable_mapping_n"], "unique_test_n": source_meta["unique_test_n"], "biomarkers": s4["biomarker"].tolist()}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
