"""Systematic human screen for the frozen outcome-blinded eligible axes.

This script is intentionally downstream of
``environmental_crc_267_actionability_prioritization.py``.  Eligibility is
read from the actionability matrix before any human association statistic is
used.  The present local audit has one eligible direct axis (MCOP), so the
across-axis BH-FDR is transparently a one-test adjustment rather than being
presented as a broad screen that has not yet been run for all 267 chemicals.
"""

from __future__ import annotations

import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = ROOT / "outputs"
MATRIX = OUTPUTS / "environmental_crc_267_actionability_matrix.csv"
PHASE2H = ROOT / "work" / "scripts" / "mcop_crc_phase2h_survey_audit.py"

OUT_RAW = OUTPUTS / "environmental_crc_systematic_human_screen.csv"
OUT_FDR = OUTPUTS / "environmental_crc_systematic_human_screen_fdr.csv"
OUT_REPORT = OUTPUTS / "environmental_crc_systematic_human_screen_report.md"
OUT_MANIFEST = OUTPUTS / "environmental_crc_systematic_human_screen_manifest.json"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bh_fdr(p_values: pd.Series) -> pd.Series:
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(float)
    result = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    if not valid.any():
        return pd.Series(result, index=p_values.index)
    order = np.argsort(values[valid])
    sorted_p = values[valid][order]
    m = len(sorted_p)
    adjusted = np.minimum.accumulate((sorted_p * m / np.arange(1, m + 1))[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    positions = np.flatnonzero(valid)[order]
    result[positions] = adjusted
    return pd.Series(result, index=p_values.index)


def build_data_with_pickle_fallback(phase2h, phase2):
    """Read the frozen frame, rebuilding it from local XPTs only if needed.

    The repository pickle was produced under a different pandas StringDtype
    implementation. Rebuilding the same harmonized frame from the versioned
    local XPT inputs is a compatibility repair, not an analytic change.
    """
    try:
        return phase2h.build_seven_cycle_frame(
            phase2,
            phase2h.DEFAULT_HARMONIZED,
            phase2h.DEFAULT_DATA_DIR,
        )
    except (NotImplementedError, TypeError, ValueError) as error:
        if not isinstance(error, NotImplementedError):
            raise
        mbzp = load_module(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "systematic_mbzp_harmonizer")
        frames = []
        input_manifest = []
        for idx, spec in enumerate(mbzp.CYCLES):
            frame, manifest_rows = mbzp.load_cycle(spec, idx)
            frames.append(frame)
            input_manifest.extend(manifest_rows)
        harmonized = pd.concat(frames, ignore_index=True)
        harmonized = harmonized[(harmonized["age"] >= 20) & harmonized["cancer_outcome_available"]].copy()
        mcop, mcop_manifest = phase2.read_mcop(phase2h.DEFAULT_DATA_DIR)
        data = phase2.build_frame(harmonized, mcop)
        data = data[data["cycle"].isin(phase2h.MCOP_CYCLES)].copy()
        data["pooled_weight"] = pd.to_numeric(data["phthalate_weight_base"], errors="coerce") / phase2h.N_CYCLES
        data["pooled_weight_rule"] = f"cycle-specific phthalate subsample weight / {phase2h.N_CYCLES}"
        data["weight_rule_audit"] = np.where(
            data["weight_source"].eq(data["cycle"].map(phase2h.EXPECTED_WEIGHT_SOURCES)),
            "pass",
            "FAIL",
        )
        if data["weight_rule_audit"].ne("pass").any():
            raise RuntimeError("Weight-source audit failed during pickle compatibility rebuild")
        return data, {"harmonized_rebuilt_from_local_xpt": True, "source_manifest_rows": len(input_manifest), "mcop_manifest": mcop_manifest}


def main() -> None:
    matrix = pd.read_csv(MATRIX, low_memory=False)
    eligible = matrix.loc[matrix["eligible_permissive"].astype(str).str.lower().eq("true")].copy()
    phase2h = load_module(PHASE2H, "systematic_mcop_phase2h")
    phase2, model, _paper, _stability = phase2h.load_components()

    # The frame is rebuilt from the same seven-cycle source and frozen primary
    # model used by the Phase 2H audit. This is not a read-through of the old
    # MCOP result file.
    data, input_manifest = build_data_with_pickle_fallback(phase2h, phase2)
    populations = model.population_frames(data)
    primary = populations["CRC_vs_cancer_free"].copy()

    rows: list[dict[str, object]] = []
    for _, candidate in eligible.iterrows():
        # Candidate-specific analyte harmonization is kept explicit.  A future
        # candidate can only enter here after a matching harmonized biomarker
        # variable and the same frozen model have been audited.
        if candidate["ChemicalID"] != "C573544" or candidate["human_biomarker"] != "MCOP":
            rows.append(
                {
                    "ChemicalID": candidate["ChemicalID"],
                    "ChemicalName": candidate["ChemicalName"],
                    "human_biomarker": candidate["human_biomarker"],
                    "status": "eligible_but_analyte_pipeline_not_registered",
                    "outcome_blinded_eligibility": True,
                    "selection_used_human_crc_effect": False,
                }
            )
            continue

        fit = phase2h.fit_primary(primary, model, "MCOP_primary_uniform_7_cycle_screen")
        rows.append(
            {
                "ChemicalID": candidate["ChemicalID"],
                "ChemicalName": candidate["ChemicalName"],
                "human_biomarker": candidate["human_biomarker"],
                "exposure_axis_name": candidate["exposure_axis_name"],
                "phase1_molecular_level": candidate["molecular_level"],
                "phase1_unfiltered_rank": candidate["phase1_unfiltered_rank"],
                "actionability_rule": "permissive",
                "status": fit.get("status"),
                "N": fit.get("N"),
                "CRC_N": fit.get("CRC_N"),
                "Control_N": fit.get("Control_N"),
                "beta": fit.get("beta"),
                "SE": fit.get("SE"),
                "OR": fit.get("OR"),
                "CI_low": fit.get("CI_low"),
                "CI_high": fit.get("CI_high"),
                "P": fit.get("P"),
                "design_df": fit.get("design_df"),
                "PSU_N": fit.get("PSU_N"),
                "strata_N": fit.get("strata_N"),
                "message": fit.get("message"),
                "cycles": ";".join(phase2h.MCOP_CYCLES),
                "weight_rule": "cycle-specific phthalate subsample weight / 7",
                "model_specification": "CRC_vs_cancer_free ~ log2(MCOP) + age + sex + race + BMI + smoking + PIR + log2(creatinine)",
                "outcome_blinded_eligibility": True,
                "selection_used_human_crc_effect": False,
            }
        )

    raw = pd.DataFrame(rows)
    raw.to_csv(OUT_RAW, index=False)
    raw["human_screen_bh_fdr"] = bh_fdr(raw["P"] if "P" in raw else pd.Series(dtype=float))
    raw.to_csv(OUT_FDR, index=False)

    tested = raw.loc[raw["status"].eq("ok")].copy()
    primary_row = tested.iloc[0] if len(tested) else None
    if primary_row is not None:
        result_sentence = (
            f"MCOP was the only currently eligible direct axis and yielded "
            f"OR={float(primary_row['OR']):.3f} (95% CI "
            f"{float(primary_row['CI_low']):.3f}–{float(primary_row['CI_high']):.3f}), "
            f"P={float(primary_row['P']):.4g}; across-axis BH-FDR is "
            f"{float(primary_row['human_screen_bh_fdr']):.4g} because one axis was tested."
        )
    else:
        result_sentence = "No currently eligible axis produced an estimable human model."

    report = f"""# Systematic human screen after 267-chemical actionability filtering

## Scope

This run uses only candidates that passed the prespecified, CRC-outcome-blinded
permissive actionability rule. The actionability matrix was frozen before this
human association was fitted. No human OR, CI, P value, LOCO result, or
cycle-specific effect entered candidate eligibility.

## Current result

{result_sentence}

The current local biomonitoring audit supports **one** direct human-testable
axis: MCOP (ChemicalID C573544). Therefore the BH-FDR reported here is a
one-test correction, not evidence that all 267 chemicals have already received
an equivalent NHANES analysis.

## Model

- Seven NHANES cycles: 2005–06 through 2017–18.
- Primary comparison: CRC versus cancer-free controls.
- Exposure: log2(MCOP), per doubling.
- Covariates: age, sex, race, BMI, smoking, PIR, and log2 urinary creatinine.
- Survey design: cycle-specific phthalate subsample weight divided by 7,
  cycle-specific strata and PSU.
- The frame is rebuilt from the harmonized source and fit through the validated
  Python Taylor-style complex-survey implementation.

## Identity firewall

MCOP is analyzed as the direct candidate/urinary analyte. It is not relabeled as
MiNP, and MiNP's molecular nomination is not used to retroactively select MCOP.
The paper's direct-discovery statement is therefore:

`267 core chemicals → outcome-blinded actionability → MCOP retained → MCOP human screen`

## Limitation

The full 267-chemical biomonitoring/actionability queue is not yet complete:
265 chemicals remain manual-review unknowns. This output supports the direct
MCOP provenance audit and its downstream human test, but a claim that MCOP was
the unique winner of a fully completed 267-axis epidemiologic screen must wait
until the remaining candidates are annotated and run under the same model.

Generated: {datetime.now(timezone.utc).isoformat()}
"""
    OUT_REPORT.write_text(report, encoding="utf-8")

    manifest = {
        "phase": "Phase 2I systematic human screen",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "actionability_matrix": str(MATRIX),
        "n_eligible_permissive": int(len(eligible)),
        "n_tested_estimable": int(len(tested)),
        "candidate_ids_tested": tested["ChemicalID"].tolist() if len(tested) else [],
        "outcome_blinded_eligibility": True,
        "selection_used_human_crc_effect": False,
        "input_manifest": input_manifest,
        "outputs": [str(OUT_RAW), str(OUT_FDR), str(OUT_REPORT), str(OUT_MANIFEST)],
        "warning": "Only one axis is currently eligible; 265 chemical actionability audits remain pending.",
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps({"eligible": len(eligible), "tested": len(tested), "result": result_sentence, "outputs": manifest["outputs"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
