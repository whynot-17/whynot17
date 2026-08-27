"""Collapse clean-room actionability results to unique human tests."""

from __future__ import annotations

import pandas as pd


def freeze_tests(actionability: pd.DataFrame, mappings: pd.DataFrame) -> pd.DataFrame:
    actionable = actionability.loc[actionability["actionable_mapping"] & actionability["human_biomarker"].ne("")].copy()
    if actionable.empty:
        return pd.DataFrame(columns=["test_id", "biomarker", "variable", "matrix", "cycles", "weight", "mapping_count", "chemical_ids", "chemical_names", "exposure_axes"])
    joined = actionable.merge(
        # cycle_list and weight_variable are actionability-stage fields.  Do
        # not merge duplicate copies from the mapping table, which would
        # create suffixed columns and obscure the frozen test definition.
        mappings[["chemical_id", "human_biomarker", "matrix", "exposure_axis"]],
        on=["chemical_id", "human_biomarker"], how="left", validate="one_to_many",
    )
    rows = []
    for variable, group in joined.groupby("human_biomarker", sort=True):
        cycles = ";".join(sorted(set(";".join(group["cycle_list"].fillna("")).split(";")) - {""}))
        rows.append({
            "test_id": f"NHANES_{variable}", "biomarker": variable, "variable": variable,
            "matrix": ";".join(sorted(set(group["matrix"].dropna().astype(str)))), "cycles": cycles,
            "weight": ";".join(sorted(set(group["weight_variable"].dropna().astype(str)) - {""})),
            "mapping_count": int(len(group)), "chemical_ids": ";".join(sorted(set(group["chemical_id"]))),
            "chemical_names": ";".join(sorted(set(group["chemical_name"]))),
            "exposure_axes": ";".join(sorted(set(group["exposure_axis"].dropna().astype(str)) - {""})),
            "n_cycles": int(group["n_cycles_available"].max()), "pooled_above_lod_pct": float(group["pooled_above_lod_pct"].iloc[0]),
        })
    return pd.DataFrame(rows).sort_values("test_id").reset_index(drop=True)
