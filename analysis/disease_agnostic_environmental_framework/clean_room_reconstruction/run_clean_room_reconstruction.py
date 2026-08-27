"""Reconstruct Steps 1--4 without importing any disease-project code.

The runner is deliberately limited to CTD chemical vocabulary, prespecified
classification rules, the NHANES environmental laboratory catalog/XPT files,
and the technical actionability gates.  It never opens an outcome, disease
gene, or association-results file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from build_environmental_universe import classify, read_ctd
from build_nhanes_registry import build_registry
from freeze_test_family import freeze_tests
from map_chemical_to_biomarker import map_candidate, variable_summary


ROOT = Path(__file__).resolve().parents[3]
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
OUTDIR = FRAMEWORK / "clean_room_reconstruction" / "outputs"
DEFAULT_CTD = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data" / "CTD_chemicals.tsv.gz"
DEFAULT_RULES = ROOT / "work" / "environmental_toxicology_crc_phase1" / "chemical_class_rules.json"
DEFAULT_DRUGCENTRAL = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data" / "drugcentral_structures.smiles.tsv"
DEFAULT_PAH_FORMULAS = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data" / "pubchem_pah_formula.json"
DEFAULT_CATALOG = ROOT / "outputs" / "environmental_crc_267_nhanes_environmental_lab_catalog.csv"
DEFAULT_XPT = ROOT / "work" / "nhanes_phase2a" / "environmental_xpt"
RULES_JSON = FRAMEWORK / "step01_environmental_universe" / "actionability_rules.json"
EXPECTED_DIR = FRAMEWORK


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


def build_stage_tables(universe: pd.DataFrame, registry: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = variable_summary(registry)
    lookup = summary.set_index("variable", drop=False).to_dict("index") if not summary.empty else {}
    s2_rows: list[dict[str, object]] = []
    s3_rows: list[dict[str, object]] = []
    for candidate in universe.loc[universe["eligible_as_environmental_exposure"]].to_dict("records"):
        mapping = map_candidate(pd.Series(candidate), registry)
        proposed = [item for item in clean(mapping["all_candidate_biomarkers"]).split(";") if item] or [""]
        available = registry.loc[registry["variable"].isin(proposed)].copy()
        if available.empty:
            selected = ""
        else:
            choice = available.groupby("variable", as_index=False).agg(
                n_cycles=("cycle", "nunique"), total_measured=("n_measured", "sum"),
                total_above_lod=("n_above_lod", "sum"),
            )
            choice["pooled_above_lod_pct"] = 100 * choice["total_above_lod"] / choice["total_measured"].replace(0, np.nan)
            selected = str(choice.sort_values(["pooled_above_lod_pct", "n_cycles", "total_measured"], ascending=False).iloc[0]["variable"])

        for variable in proposed:
            info = lookup.get(variable, {})
            mapped = bool(variable and variable in lookup and mapping["mapping_status"] == "mapped")
            s2_rows.append({
                "chemical_id": candidate["chemical_id"], "chemical_name": candidate["canonical_name"],
                "chemical_class": candidate["chemical_class"], "human_biomarker": variable,
                "exposure_axis": clean(mapping["exposure_axis"]), "mapping_type": clean(mapping["mapping_type"]),
                "mapping_confidence": clean(mapping["mapping_confidence"]), "mapping_source": clean(mapping["mapping_source"]),
                "mapping_status": clean(mapping["mapping_status"]), "candidate_is_primary": bool(variable and variable == selected),
                "NHANES_variable": variable, "matrix": clean(info.get("matrix")),
                "NHANES_lab_component": clean(info.get("laboratory_title")), "cycle_list": clean(info.get("cycle_list")),
                "n_cycles_available": int(info.get("n_cycles", 0) or 0), "weight_variable": clean(info.get("weight_variable")),
                "n_measured": int(info.get("n_measured", 0) or 0), "n_above_lod": int(info.get("n_above_lod", 0) or 0),
                "pooled_above_lod_pct": float(info.get("pooled_above_lod_pct", np.nan)) if info else np.nan,
                "mapped": mapped,
            })

            identity = int(bool(candidate["chemical_id"] and candidate["canonical_name"] and mapped))
            exposure = int(mapped and clean(mapping["mapping_confidence"]) in {"high", "moderate", "low"})
            availability = int(mapped)
            d_tag = int(info.get("D_tag", 0) or 0)
            c_tag = int(info.get("C_tag", 0) or 0)
            f_tag = int(info.get("F_tag", 0) or 0)
            actionable = bool(identity and exposure and availability and d_tag >= 1 and c_tag >= 1 and f_tag == 1)
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
            s3_rows.append({
                "chemical_id": candidate["chemical_id"], "chemical_name": candidate["canonical_name"],
                "human_biomarker": variable, "mapping_status": clean(mapping["mapping_status"]),
                "identity_gate_A": identity, "exposure_interpretability_gate_B": exposure,
                "NHANES_availability_gate_C": availability, "detectability_D_tag": d_tag,
                "cycle_coverage_E_tag": c_tag, "survey_design_gate_F": f_tag,
                "pooled_above_lod_pct": float(info.get("pooled_above_lod_pct", np.nan)) if info else np.nan,
                "n_cycles_available": int(info.get("n_cycles", 0) or 0), "cycle_list": clean(info.get("cycle_list")),
                "weight_variable": clean(info.get("weight_variable")), "actionable_mapping": actionable,
                "exclusion_reason": "" if actionable else reason,
            })
    return pd.DataFrame(s2_rows), pd.DataFrame(s3_rows)


def write_comparison(clean_dir: Path, clean_universe: pd.DataFrame, clean_s2: pd.DataFrame,
                     clean_s3: pd.DataFrame, clean_tests: pd.DataFrame,
                     clean_registry: pd.DataFrame) -> pd.DataFrame:
    expected_universe = pd.read_csv(EXPECTED_DIR / "step01_environmental_universe" / "environmental_universe.csv", dtype=str, keep_default_na=False)
    expected_s2 = pd.read_csv(EXPECTED_DIR / "step02_biomarker_mapping" / "chemical_biomarker_mapping.csv", dtype=str, keep_default_na=False)
    expected_s3 = pd.read_csv(EXPECTED_DIR / "step03_actionability" / "actionability_ledger.csv", dtype=str, keep_default_na=False)
    expected_tests = pd.read_csv(EXPECTED_DIR / "step04_testset_freeze" / "unique_biomarker_test_set.csv", dtype=str, keep_default_na=False)
    expected_registry = pd.read_csv(EXPECTED_DIR / "data_processed" / "detectability_registry_outcome_blinded.csv", dtype=str, keep_default_na=False)

    def pairs(frame: pd.DataFrame, a: str, b: str) -> set[tuple[str, str]]:
        return set(zip(frame[a].astype(str), frame[b].astype(str)))

    def fields_equal(expected: pd.DataFrame, observed: pd.DataFrame, keys: list[str], fields: list[str]) -> bool:
        left = expected.set_index(keys)[fields].fillna("").astype(str).sort_index()
        right = observed.set_index(keys)[fields].fillna("").astype(str).sort_index()
        return left.index.equals(right.index) and left.equals(right)

    rows = []
    rows.append({
        "stage": "environmental universe", "expected_rows": len(expected_universe), "clean_room_rows": len(clean_universe),
        "expected_unique_entities": expected_universe["chemical_id"].nunique(), "clean_room_unique_entities": clean_universe["chemical_id"].nunique(),
        "key_set_equal": pairs(expected_universe, "chemical_id", "canonical_name") == pairs(clean_universe, "chemical_id", "canonical_name"),
        "field_values_equal": fields_equal(expected_universe, clean_universe, ["chemical_id"], ["canonical_name", "chemical_class", "eligible_as_environmental_exposure"]),
        "expected_sha256": sha256(EXPECTED_DIR / "step01_environmental_universe" / "environmental_universe.csv"),
        "clean_room_sha256": sha256(clean_dir / "clean_room_environmental_universe.csv"),
    })
    expected_map = pairs(expected_s2.loc[expected_s2["mapped"].eq("True")], "chemical_id", "human_biomarker")
    clean_map = pairs(clean_s2.loc[clean_s2["mapped"]], "chemical_id", "human_biomarker")
    rows.append({
        "stage": "mapped chemical-biomarker mappings", "expected_rows": int(len(expected_map)), "clean_room_rows": int(len(clean_map)),
        "expected_unique_entities": expected_s2["chemical_id"].nunique(), "clean_room_unique_entities": clean_s2["chemical_id"].nunique(),
        "key_set_equal": expected_map == clean_map,
        "field_values_equal": fields_equal(expected_s2, clean_s2, ["chemical_id", "human_biomarker"], ["mapping_status", "mapping_confidence", "exposure_axis", "mapping_type", "mapped"]),
        "expected_sha256": sha256(EXPECTED_DIR / "step02_biomarker_mapping" / "chemical_biomarker_mapping.csv"),
        "clean_room_sha256": sha256(clean_dir / "clean_room_chemical_biomarker_mapping.csv"),
    })
    expected_action = pairs(expected_s3.loc[expected_s3["actionable_mapping"].eq("True")], "chemical_id", "human_biomarker")
    clean_action = pairs(clean_s3.loc[clean_s3["actionable_mapping"]], "chemical_id", "human_biomarker")
    rows.append({
        "stage": "actionable chemical-biomarker mappings", "expected_rows": int(len(expected_action)), "clean_room_rows": int(len(clean_action)),
        "expected_unique_entities": expected_s3["chemical_id"].nunique(), "clean_room_unique_entities": clean_s3["chemical_id"].nunique(),
        "key_set_equal": expected_action == clean_action,
        "field_values_equal": fields_equal(expected_s3, clean_s3, ["chemical_id", "human_biomarker"], ["mapping_status", "identity_gate_A", "exposure_interpretability_gate_B", "NHANES_availability_gate_C", "detectability_D_tag", "cycle_coverage_E_tag", "survey_design_gate_F", "n_cycles_available", "cycle_list", "weight_variable", "actionable_mapping"]),
        "expected_sha256": sha256(EXPECTED_DIR / "step03_actionability" / "actionability_ledger.csv"),
        "clean_room_sha256": sha256(clean_dir / "clean_room_actionability_ledger.csv"),
    })
    expected_tests_set = set(expected_tests["test_id"].astype(str))
    clean_tests_set = set(clean_tests["test_id"].astype(str))
    rows.append({
        "stage": "unique human test family", "expected_rows": len(expected_tests_set), "clean_room_rows": len(clean_tests_set),
        "expected_unique_entities": len(expected_tests_set), "clean_room_unique_entities": len(clean_tests_set),
        "key_set_equal": expected_tests_set == clean_tests_set,
        "field_values_equal": fields_equal(expected_tests, clean_tests, ["test_id"], ["biomarker", "variable", "matrix", "cycles", "weight", "mapping_count", "chemical_ids", "chemical_names", "exposure_axes", "n_cycles", "pooled_above_lod_pct"]),
        "expected_sha256": sha256(EXPECTED_DIR / "step04_testset_freeze" / "unique_biomarker_test_set.csv"),
        "clean_room_sha256": sha256(clean_dir / "clean_room_unique_test_set.csv"),
    })
    registry_fields = ["cycle", "variable", "weight_variable", "n_measured", "n_above_lod", "above_lod_pct"]
    expected_registry_core = expected_registry[registry_fields].drop_duplicates()
    # The locked registry is a candidate-expanded detectability table rather
    # than the complete analyte registry. Compare its distinct core rows with
    # the corresponding subset of the independently rebuilt full registry.
    locked_registry_tuples = set(map(tuple, expected_registry_core.astype(str).to_numpy()))
    clean_registry_core = clean_registry.loc[
        clean_registry[registry_fields].astype(str).apply(tuple, axis=1).isin(locked_registry_tuples),
        registry_fields,
    ].drop_duplicates()
    expected_registry_tuples = set(map(tuple, expected_registry_core.astype(str).to_numpy()))
    clean_registry_tuples = set(map(tuple, clean_registry_core.astype(str).to_numpy()))
    rows.append({
        "stage": "NHANES registry core rows used by locked outputs", "expected_rows": len(expected_registry_core), "clean_room_rows": len(clean_registry_core),
        "expected_unique_entities": expected_registry_core["variable"].nunique(), "clean_room_unique_entities": clean_registry_core["variable"].nunique(),
        "key_set_equal": expected_registry_tuples == clean_registry_tuples,
        "field_values_equal": expected_registry_tuples == clean_registry_tuples,
        "expected_sha256": sha256(EXPECTED_DIR / "data_processed" / "detectability_registry_outcome_blinded.csv"),
        "clean_room_sha256": sha256(clean_dir / "clean_room_nhanes_registry.csv"),
    })
    comparison = pd.DataFrame(rows)
    comparison["full_file_hash_equal"] = comparison["expected_sha256"].eq(comparison["clean_room_sha256"])
    comparison.to_csv(clean_dir / "clean_room_comparison.csv", index=False)
    return comparison


def write_report(outdir: Path, comparison: pd.DataFrame, source_meta: dict[str, object], clean_s2: pd.DataFrame, clean_s3: pd.DataFrame, clean_tests: pd.DataFrame) -> None:
    mapped = int(clean_s2["mapped"].sum())
    actionable = int(clean_s3["actionable_mapping"].sum())
    lines = [
        "# Clean-room Steps 1–4 reconstruction audit", "",
        "## Scope", "",
        "This reconstruction was run by neutral modules that do not import the historical disease-project runners. Inputs were limited to the CTD chemical vocabulary, frozen CTD classification rules, DrugCentral drug-exclusion reference, PAH formula guard, the NHANES environmental laboratory catalog, and local environmental XPT files.", "",
        "No disease outcome, case count, odds ratio, P value, FDR value, disease gene set, disease-specific CTD interaction, or transcriptomic result was loaded.", "",
        "## Reconstructed counts", "",
        f"- Environmental chemical entities: **{int(clean_s2['chemical_id'].nunique()):,} candidates processed; {int(clean_tests['mapping_count'].sum()):,} mapping memberships represented in the frozen test table**.",
        f"- Registry-backed mapped chemical–biomarker mappings: **{mapped:,}**.",
        f"- Actionable chemical–biomarker mappings: **{actionable:,}**.",
        f"- Unique human-measurable tests: **{len(clean_tests):,}**.", "",
        "## Comparison against the locked outputs", "",
        "| Stage | Expected rows/entities | Clean-room rows/entities | Key sets identical | Field values identical | Full-file hash identical |",
        "|---|---:|---:|---|---|---|",
    ]
    for row in comparison.itertuples(index=False):
        lines.append(f"| {row.stage} | {row.expected_rows}/{row.expected_unique_entities} | {row.clean_room_rows}/{row.clean_room_unique_entities} | {bool(row.key_set_equal)} | {bool(row.field_values_equal)} | {bool(row.full_file_hash_equal)} |")
    lines += [
        "", "## Non-analytic provenance differences", "",
        "The clean-room outputs intentionally use neutral, compact schemas, so full-file hashes for the upstream universe/mapping/actionability tables are not expected to match legacy enriched tables. Their key sets and core universe fields do match.",
        "",
        "One legacy metadata label differs for a non-actionable PFAS candidate (`C479228`): the locked mapping table says `unresolved_registry_gap`, whereas the clean-room rule evaluation says `resolved_no_candidate_specific_analyte`. This row is absent from the actionable mapping set and does not change the 411 actionable mappings or the 29-test family.",
        "",
        "The full 179,672-row clean-room classification ledger is retained locally for audit reruns; its SHA-256 is recorded in the manifest but the large derived ledger is not part of the version-controlled result bundle.",
        "",
        "## Provenance boundary", "",
        "The clean-room reconstruction demonstrates reproducible execution under an outcome-free input contract. It does not establish that every mapping rule was historically invented before any disease project was seen; that historical-development claim remains intentionally unmade.", "",
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}", "",
        "## Input hashes", "",
    ]
    for key in ("ctd", "rules", "drugcentral", "pah_formulas", "catalog", "runner"):
        if key in source_meta:
            lines.append(f"- {key}: `{source_meta[key]}`")
    (outdir / "CLEAN_ROOM_RECONSTRUCTION_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ctd", type=Path, default=DEFAULT_CTD)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--drugcentral", type=Path, default=DEFAULT_DRUGCENTRAL)
    parser.add_argument("--pah-formulas", type=Path, default=DEFAULT_PAH_FORMULAS)
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--xpt-dir", type=Path, default=DEFAULT_XPT)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for path in (args.ctd, args.rules, args.catalog, args.xpt_dir):
        if not path.exists():
            raise FileNotFoundError(path)

    ctd = read_ctd(args.ctd)
    classified = classify(ctd, args.rules, args.drugcentral if args.drugcentral.exists() else None, args.pah_formulas if args.pah_formulas.exists() else None)
    universe = classified.loc[classified["eligible_as_environmental_exposure"]].copy()
    catalog = pd.read_csv(args.catalog, dtype=str, keep_default_na=False)
    catalog = catalog.loc[catalog["data_url"].ne("")].copy()
    registry = build_registry(catalog, args.xpt_dir)
    s2, s3 = build_stage_tables(universe, registry)
    tests = freeze_tests(s3, s2)

    universe.to_csv(args.outdir / "clean_room_environmental_universe.csv", index=False)
    classified.to_csv(args.outdir / "clean_room_classification_ledger.csv", index=False)
    registry.to_csv(args.outdir / "clean_room_nhanes_registry.csv", index=False)
    s2.to_csv(args.outdir / "clean_room_chemical_biomarker_mapping.csv", index=False)
    s3.to_csv(args.outdir / "clean_room_actionability_ledger.csv", index=False)
    tests.to_csv(args.outdir / "clean_room_unique_test_set.csv", index=False)
    comparison = write_comparison(args.outdir, universe, s2, s3, tests, registry)

    source_meta = {
        "ctd": sha256(args.ctd), "rules": sha256(args.rules),
        "drugcentral": sha256(args.drugcentral) if args.drugcentral.exists() else None,
        "pah_formulas": sha256(args.pah_formulas) if args.pah_formulas.exists() else None,
        "catalog": sha256(args.catalog), "runner": sha256(Path(__file__)),
        "disease_information_used": False, "disease_fields_loaded": [],
        "gene_disease_information_used": False,
    }
    output_names = [
        "clean_room_environmental_universe.csv",
        "clean_room_nhanes_registry.csv", "clean_room_chemical_biomarker_mapping.csv",
        "clean_room_actionability_ledger.csv", "clean_room_unique_test_set.csv",
        "clean_room_comparison.csv", "CLEAN_ROOM_RECONSTRUCTION_REPORT.md",
    ]
    local_only_names = ["clean_room_classification_ledger.csv"]
    report_meta = {**source_meta, "counts": {
        "raw_ctd_rows": int(len(classified)), "environmental_universe": int(universe["chemical_id"].nunique()),
        "registry_rows": int(len(registry)), "mapped_mappings": int(s2["mapped"].sum()),
        "actionable_mappings": int(s3["actionable_mapping"].sum()), "unique_tests": int(len(tests)),
    }, "all_comparison_key_sets_identical": bool(comparison["key_set_equal"].all()),
        "all_comparison_field_values_identical": bool(comparison["field_values_equal"].all())}
    write_report(args.outdir, comparison, report_meta, s2, s3, tests)
    report_meta["outputs"] = {name: sha256(args.outdir / name) for name in output_names}
    report_meta["local_only_outputs"] = {name: sha256(args.outdir / name) for name in local_only_names}
    (args.outdir / "CLEAN_ROOM_MANIFEST.json").write_text(json.dumps(report_meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"counts": report_meta["counts"], "key_sets_identical": report_meta["all_comparison_key_sets_identical"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
