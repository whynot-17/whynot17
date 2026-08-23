"""Outcome-blinded actionability followed by the unified 267-candidate human screen.

The actionability matrix is created before this script is called.  Every
permissive eligible chemical is represented in the member list, while the
model is fit once per unique exposure-axis/primary-NHANES-analyte combination
to avoid counting the same biomarker as multiple independent hypotheses.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = ROOT / "outputs"
XPT_DIR = ROOT / "work" / "nhanes_phase2a" / "environmental_xpt"
MATRIX = OUTPUTS / "environmental_crc_267_actionability_matrix_v2.csv"
DETECT = OUTPUTS / "environmental_crc_267_detectability_by_cycle.csv"
OUT_RESULTS = OUTPUTS / "environmental_crc_systematic_human_screen_v2.csv"
OUT_FDR = OUTPUTS / "environmental_crc_systematic_human_screen_fdr_v2.csv"
OUT_REPORT = OUTPUTS / "environmental_crc_systematic_human_screen_report_v2.md"
OUT_MANIFEST = OUTPUTS / "environmental_crc_systematic_human_screen_manifest_v2.json"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh_fdr(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg q-values over all finite tested axes."""
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = p.notna() & np.isfinite(p)
    if not valid.any():
        return out
    order = p.loc[valid].sort_values().index
    ranked = p.loc[order].to_numpy(float)
    m = len(ranked)
    q = np.minimum(1.0, ranked * m / np.arange(1, m + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out.loc[order] = q
    return out


def load_harmonized(audit_module) -> pd.DataFrame:
    frame = audit_module.load_harmonized()
    if len(frame) == 0:
        raise RuntimeError("The harmonized NHANES CRC frame is empty")
    frame["SEQN"] = pd.to_numeric(frame["SEQN"], errors="coerce").astype("Int64")
    return frame


def read_environmental_axis(axis_rows: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    """Read one selected analyte with one weighted XPT record per cycle."""
    rows = axis_rows.drop_duplicates(["cycle", "data_file", "variable"]).copy()
    usable = rows.loc[rows["weight_variable"].notna() & rows["weight_variable"].ne("")].copy()
    cycles = sorted(usable["cycle"].dropna().astype(str).unique().tolist())
    if not cycles:
        return pd.DataFrame(), {"status": "not_estimable", "reason": "no weighted analyte file"}
    pieces = []
    source_rows = []
    for row in usable.itertuples(index=False):
        path = Path(row.local_xpt)
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_sas(path, format="xport", encoding="latin1")
        needed = ["SEQN", row.variable, row.weight_variable]
        if not set(needed).issubset(frame.columns):
            continue
        part = frame[needed].copy()
        part["SEQN"] = pd.to_numeric(part["SEQN"], errors="coerce").astype("Int64")
        part["exposure_raw"] = pd.to_numeric(part[row.variable], errors="coerce")
        part["survey_weight"] = pd.to_numeric(part[row.weight_variable], errors="coerce")
        # CDC's 1999-2000/2001-2002 phthalate 4-year weight covers twice the
        # width of a single 2-year cycle. Match the frozen harmonizer: apply
        # the 2x conversion only when the file actually carries WTSPH4YR.
        weight_multiplier = 2.0 if str(row.cycle) in {"1999-2000", "2001-2002"} and str(row.weight_variable) == "WTSPH4YR" else 1.0
        part["survey_weight"] = part["survey_weight"] * weight_multiplier
        part["cycle"] = str(row.cycle)
        part = part[["SEQN", "cycle", "exposure_raw", "survey_weight"]]
        pieces.append(part)
        source_rows.append({"cycle": row.cycle, "data_file": row.data_file, "variable": row.variable, "weight_variable": row.weight_variable, "weight_multiplier": weight_multiplier})
    if not pieces:
        return pd.DataFrame(), {"status": "not_estimable", "reason": "no readable weighted analyte file"}
    exposure = pd.concat(pieces, ignore_index=True)
    # Supplemental duplicates have already been removed by data_file; protect
    # the participant-level model from any residual duplicate SEQN/cycle.
    exposure = exposure.sort_values(["cycle", "SEQN"]).drop_duplicates(["cycle", "SEQN"], keep="first")
    exposure["axis_log2"] = np.log2(exposure["exposure_raw"].where(exposure["exposure_raw"] > 0))
    exposure["pooled_weight"] = exposure["survey_weight"] / len(cycles)
    return exposure, {"status": "ok", "cycles": cycles, "source_rows": source_rows, "n_raw": len(exposure)}


def fit_axis(axis: pd.Series, axis_rows: pd.DataFrame, harmonized: pd.DataFrame, model) -> dict[str, object]:
    variable = str(axis["selected_primary_biomarker"])
    exposure, source = read_environmental_axis(axis_rows)
    base = {
        "axis_key": axis["axis_key"],
        "exposure_axis": axis["exposure_axis"],
        "primary_biomarker": variable,
        "biomarker_type": axis.get("biomarker_type", ""),
        "biological_matrix": axis.get("biological_matrix", ""),
        "eligible_chemical_count": int(axis["eligible_chemical_count"]),
        "eligible_chemical_ids": axis["eligible_chemical_ids"],
        "eligible_chemical_names": axis["eligible_chemical_names"],
        "mapping_confidence": axis.get("mapping_confidence", ""),
        "n_cycles_available": int(axis.get("n_cycles_available", 0)),
        "cycle_list": axis.get("cycle_list", ""),
        "source_registry_status": source.get("status"),
        "source_registry_n": source.get("n_raw", 0),
    }
    if exposure.empty:
        return {**base, "status": "not_estimable", "reason": source.get("reason", "empty exposure frame")}
    frame = exposure.merge(harmonized, on=["SEQN", "cycle"], how="inner", validate="one_to_one", suffixes=("", "_outcome"))
    # Use the frozen CRC-vs-cancer-free population construction. The current
    # model gets the exposure variable through a generic axis_log2 column.
    urine = "urine" in str(axis.get("biological_matrix", "")).lower()
    continuous = ["axis_log2", "age", "bmi", "pir"]
    if urine:
        continuous.append("creatinine_log2")
    categorical = ["sex", "race", "smoking"]
    population = model.population_frames(frame)["CRC_vs_cancer_free"]
    fit = model.fit_survey_logistic(
        population,
        continuous,
        categorical,
        exposure_name="axis_log2",
        levels=model.LEVELS,
    )
    clean = {k: v for k, v in fit.items() if k not in {"coefficients", "covariance"}}
    # Avoid case-insensitive duplicate CSV headers (PowerShell and some
    # spreadsheet readers treat Control_N and control_n as the same field).
    clean = {
        ("fit_N" if k == "N" else "fit_crc_cases" if k == "CRC_N" else "fit_control_n" if k == "Control_N" else k): v
        for k, v in clean.items()
    }
    complete_req = ["outcome", "axis_log2", "age", "bmi", "pir", *categorical, "pooled_weight", "psu", "strata"]
    if urine:
        complete_req.append("creatinine_log2")
    complete = population.dropna(subset=complete_req)
    complete = complete[complete["pooled_weight"].gt(0)]
    return {
        **base,
        **clean,
        "analytic_n": int(len(complete)),
        "crc_cases": int(complete["outcome"].sum()) if len(complete) else 0,
        "control_n": int(len(complete) - complete["outcome"].sum()) if len(complete) else 0,
        "model_specification": "CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR" + (" + log2(creatinine)" if urine else ""),
        "weight_rule": f"selected analyte survey weight / {len(source.get('cycles', []))} cycles",
    }


def main() -> None:
    matrix = pd.read_csv(MATRIX, low_memory=False)
    detect = pd.read_csv(DETECT, low_memory=False)
    if len(matrix) != 267 or matrix["ChemicalID"].nunique() != 267:
        raise ValueError("The systematic screen requires the complete 267-row actionability matrix")
    truth = matrix["eligible_permissive"].astype(str).str.lower().eq("true")
    eligible = matrix.loc[truth].copy()
    # One test per actual NHANES biomarker axis. Parent chemicals that resolve
    # to the same analyte (even if their descriptive exposure_axis labels
    # differ) are members of the same hypothesis and share one model/FDR slot.
    eligible["axis_key"] = eligible["biological_matrix"].fillna("") + "|" + eligible["selected_primary_biomarker"].fillna("")
    axis_meta = eligible.groupby("axis_key", as_index=False).agg(
        exposure_axis=("exposure_axis", lambda s: ";".join(sorted(set(s.astype(str))))),
        selected_primary_biomarker=("selected_primary_biomarker", "first"),
        biomarker_type=("biomarker_type", "first"),
        biological_matrix=("biological_matrix", "first"),
        mapping_confidence=("mapping_confidence", "first"),
        eligible_chemical_count=("ChemicalID", "size"),
        eligible_chemical_ids=("ChemicalID", lambda s: ";".join(s.astype(str))),
        eligible_chemical_names=("ChemicalName", lambda s: ";".join(s.astype(str))),
        n_cycles_available=("n_cycles_available", "first"),
        cycle_list=("cycle_list", "first"),
    )
    audit_module = load_module(ROOT / "work" / "scripts" / "environmental_crc_267_biomarker_audit_v2.py", "screen_audit")
    model = load_module(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "screen_model")
    harmonized = load_harmonized(audit_module)

    rows = []
    for _, axis in axis_meta.iterrows():
        selected = axis["selected_primary_biomarker"]
        axis_rows = detect.loc[detect["variable"].eq(selected)].copy()
        rows.append(fit_axis(axis, axis_rows, harmonized, model))
    results = pd.DataFrame(rows)
    results["BH_FDR"] = bh_fdr(results.get("P", pd.Series(np.nan, index=results.index)))
    results["tested_in_unified_screen"] = True
    results["outcome_blinded_for_eligibility"] = True
    results.to_csv(OUT_RESULTS, index=False)
    fdr = results.sort_values(["BH_FDR", "P"], na_position="last").reset_index(drop=True)
    fdr["screen_rank"] = np.arange(1, len(fdr) + 1)
    fdr.to_csv(OUT_FDR, index=False)

    finite = results["P"].notna() & np.isfinite(pd.to_numeric(results["P"], errors="coerce"))
    mcop = results.loc[results["primary_biomarker"].eq("URXCOP")]
    lines = [
        "# Systematic NHANES human screen v2",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Frozen scope",
        "",
        f"The complete actionability audit covers 267 original chemicals. {len(eligible)} chemical rows satisfy the permissive outcome-blinded rule and collapse to {len(axis_meta)} unique exposure-axis/primary-analyte tests. Every one of these {len(axis_meta)} axes was passed through the same complex-survey logistic model; BH-FDR was calculated over all finite P values from this unified screen.",
        "",
        "Model: `CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR`; urinary analytes additionally include `log2(creatinine)`. Serum/blood analytes do not receive creatinine adjustment.",
        "",
        f"Finite fitted P values: {int(finite.sum())}/{len(results)}. Survey weights were taken from the selected analyte's own NHANES file and divided by the number of included cycles for that analyte. The existing validated Taylor-style PSU/strata implementation was used.",
        "",
        "## Results",
        "",
        fdr[[c for c in ["screen_rank", "exposure_axis", "primary_biomarker", "eligible_chemical_count", "analytic_n", "crc_cases", "OR", "CI_low", "CI_high", "P", "BH_FDR", "status"] if c in fdr.columns]].to_markdown(index=False),
        "",
        "## MCOP and MiNP interpretation",
        "",
        "MCOP (`URXCOP`) is included as one of the prespecified eligible axes. MiNP (`URXMNP`) remains a separate DINP molecular nomination but is not a primary human-screen axis because direct MiNP detectability fails the frozen D>=1 gate; it was not removed or merged into MCOP.",
        "",
        "This screen is an epidemiologic association scan, not causal evidence. Candidate selection/eligibility was frozen before reading these OR/P/FDR results.",
    ]
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest = {
        "analysis": "Complete 267-chemical systematic NHANES human screen v2",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_actionability_matrix": str(MATRIX),
        "input_actionability_matrix_sha256": sha256(MATRIX),
        "n_original_chemicals": int(len(matrix)),
        "n_permissive_eligible_chemical_rows": int(len(eligible)),
        "n_unique_exposure_axes_tested": int(len(axis_meta)),
        "n_finite_p_values": int(finite.sum()),
        "bh_fdr_scope": "all finite P values from all unique eligible exposure axes",
        "model": "validated Python Taylor-style complex-survey logistic; common model with creatinine only for urine",
        "outcome_blinded_for_eligibility": True,
        "outputs": [str(OUT_RESULTS), str(OUT_FDR), str(OUT_REPORT)],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"n_chemical_rows": len(eligible), "n_axes": len(axis_meta), "n_finite_p": int(finite.sum())}, ensure_ascii=False))


if __name__ == "__main__":
    main()
