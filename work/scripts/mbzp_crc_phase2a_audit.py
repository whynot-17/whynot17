"""Phase 2A: NHANES MBzP/CRC sample audit only.

This script deliberately stops before any regression or association model.
Raw NHANES XPT files are local/ignored; the cycle-specific file map below
prevents accidental use of the wrong phthalate file for early cycles.
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "work" / "nhanes_phase2a" / "data"
OUTPUT_DIR = REPO_ROOT / "outputs"
CRC_TYPE_CODES = {16: "colon", 31: "rectal"}

CYCLE_SPECS = [
    # Local names are normalized to cycle + component; the manifest preserves them.
    ("1999-2000", "PHPYPA", "MCQ", "DEMO"),
    ("2001-2002", "PHPYPA", "MCQ", "DEMO"),  # source components PHPYPA_B/MCQ_B/DEMO_B
    ("2003-2004", "L24PH", "MCQ", "DEMO"),  # source components L24PH_C/MCQ_C/DEMO_C
    ("2005-2006", "PHTHTE", "MCQ", "DEMO"),
    ("2007-2008", "PHTHTE", "MCQ", "DEMO"),
    ("2009-2010", "PHTHTE", "MCQ", "DEMO"),
    ("2011-2012", "PHTHTE", "MCQ", "DEMO"),
    ("2013-2014", "PHTHTE", "MCQ", "DEMO"),
    ("2015-2016", "PHTHTE", "MCQ", "DEMO"),
    ("2017-2018", "PHTHTE", "MCQ", "DEMO"),
]


def read_xpt(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing local NHANES file: {path}")
    return pd.read_sas(path, format="xport", encoding="latin1")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def crc_flags(mcq: pd.DataFrame) -> pd.DataFrame:
    required = {"SEQN", "MCQ220"}
    missing = required - set(mcq.columns)
    if missing:
        raise ValueError(f"MCQ file missing required columns: {sorted(missing)}")
    type_cols = [c for c in mcq.columns if c.upper() in {"MCQ230A", "MCQ230B", "MCQ230C", "MCQ230D"}]
    if not type_cols:
        raise ValueError("MCQ file has no MCQ230A-D cancer type columns")
    type_frame = mcq[type_cols].apply(pd.to_numeric, errors="coerce")
    colon = type_frame.eq(16).any(axis=1)
    rectal = type_frame.eq(31).any(axis=1)
    return pd.DataFrame(
        {
            "SEQN": mcq["SEQN"],
            "cancer_outcome_available": mcq["MCQ220"].notna(),
            "cancer_reported": pd.to_numeric(mcq["MCQ220"], errors="coerce").eq(1),
            "cancer_free": pd.to_numeric(mcq["MCQ220"], errors="coerce").eq(2),
            "colon": colon,
            "rectal": rectal,
            "crc": colon | rectal,
            "both_colon_rectal": colon & rectal,
        }
    )


def first_weight_column(lab: pd.DataFrame) -> str | None:
    for column in ("WTSPH2YR", "WTSB2YR", "WTSA2YR", "WTSPH4YR"):
        if column in lab.columns:
            return column
    return None


def audit_cycle(cycle: str, lab_prefix: str, mcq_prefix: str, demo_prefix: str):
    lab_path = DATA_DIR / f"{cycle}_{lab_prefix}.XPT"
    mcq_path = DATA_DIR / f"{cycle}_{mcq_prefix}.XPT"
    demo_path = DATA_DIR / f"{cycle}_{demo_prefix}.XPT"
    lab = read_xpt(lab_path)
    mcq = read_xpt(mcq_path)
    demo = read_xpt(demo_path)
    if "SEQN" not in lab.columns or "URXMZP" not in lab.columns:
        raise ValueError(f"{lab_path.name} lacks SEQN/URXMZP")
    if "SEQN" not in demo.columns:
        raise ValueError(f"{demo_path.name} lacks SEQN")
    flags = crc_flags(mcq)
    merged = lab[["SEQN", "URXMZP"]].merge(flags, on="SEQN", how="inner", validate="one_to_one")
    mbzp = merged["URXMZP"].notna()
    mbzp_outcome = mbzp & merged["cancer_outcome_available"]
    result = {
        "cycle": cycle,
        "nhanes_total_n": int(len(demo)),
        "mbzp_available_n": int(mbzp.sum()),
        "cancer_outcome_available_n": int(merged["cancer_outcome_available"].sum()),
        "mbzp_crc_outcome_n": int(mbzp_outcome.sum()),
        "crc_cases_n": int((mbzp_outcome & merged["crc"]).sum()),
        "colon_cases_n": int((mbzp_outcome & merged["colon"]).sum()),
        "rectal_cases_n": int((mbzp_outcome & merged["rectal"]).sum()),
        "both_colon_rectal_cases_n": int((mbzp_outcome & merged["both_colon_rectal"]).sum()),
        "controls_non_crc_n": int((mbzp_outcome & ~merged["crc"]).sum()),
        "cancer_free_controls_n": int((mbzp_outcome & merged["cancer_free"]).sum()),
        "lab_rows_n": int(len(lab)),
        "mcq_rows_n": int(len(mcq)),
        "demo_rows_n": int(len(demo)),
        "mbzp_weight_column": first_weight_column(lab) or "",
        "lab_file": lab_path.name,
        "mcq_file": mcq_path.name,
        "demo_file": demo_path.name,
        "crc_definition": "MCQ230A-D contains 16 (Colon) or 31 (Rectum); combined CRC is the union",
    }
    manifest_files = []
    for path in (lab_path, mcq_path, demo_path):
        manifest_files.append({"file": path.name, "sha256": file_sha256(path), "bytes": path.stat().st_size})
    return result, manifest_files


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cycle_rows = []
    source_files = []
    for spec in CYCLE_SPECS:
        row, files = audit_cycle(*spec)
        cycle_rows.append(row)
        source_files.extend(files)
    cycle_df = pd.DataFrame(cycle_rows)
    cycle_df.to_csv(OUTPUT_DIR / "mbzp_crc_phase2_cycle_audit.csv", index=False)

    sum_columns = [
        "nhanes_total_n", "mbzp_available_n", "cancer_outcome_available_n", "mbzp_crc_outcome_n",
        "crc_cases_n", "colon_cases_n", "rectal_cases_n", "both_colon_rectal_cases_n",
        "controls_non_crc_n", "cancer_free_controls_n",
    ]
    totals = cycle_df[sum_columns].sum().astype(int).to_dict()
    sample_rows = [
        {"metric": "nhanes_total_n", "n": totals["nhanes_total_n"], "definition": "Sum of DEMO participant rows across the ten 2-year cycles; unweighted."},
        {"metric": "mbzp_available_n", "n": totals["mbzp_available_n"], "definition": "URXMZP is non-missing in the phthalate laboratory file; unweighted."},
        {"metric": "cancer_outcome_available_n", "n": totals["cancer_outcome_available_n"], "definition": "MCQ220 is non-missing after lab–MCQ SEQN intersection; unweighted."},
        {"metric": "mbzp_crc_outcome_n", "n": totals["mbzp_crc_outcome_n"], "definition": "Both URXMZP and MCQ220 are non-missing; final Phase 2A audit denominator; unweighted."},
        {"metric": "crc_cases_n", "n": totals["crc_cases_n"], "definition": "Combined CRC cases among MBzP + outcome participants; MCQ230A-D code 16 or 31; unweighted."},
        {"metric": "colon_cases_n", "n": totals["colon_cases_n"], "definition": "At least one MCQ230A-D code 16 (Colon); unweighted."},
        {"metric": "rectal_cases_n", "n": totals["rectal_cases_n"], "definition": "At least one MCQ230A-D code 31 (Rectum); unweighted."},
        {"metric": "both_colon_rectal_cases_n", "n": totals["both_colon_rectal_cases_n"], "definition": "Participants with both a colon and rectal cancer type code; unweighted."},
        {"metric": "controls_non_crc_n", "n": totals["controls_non_crc_n"], "definition": "MBzP + outcome participants without a combined CRC type code; unweighted."},
        {"metric": "cancer_free_controls_n", "n": totals["cancer_free_controls_n"], "definition": "MBzP + MCQ220=2 (reported no cancer); unweighted."},
        {"metric": "model_1_complete_case_n", "n": "", "definition": "Deferred: Phase 2A audit only; no covariate model run."},
        {"metric": "model_2_complete_case_n", "n": "", "definition": "Deferred: Phase 2A audit only; no covariate model run."},
        {"metric": "model_3_complete_case_n", "n": "", "definition": "Deferred: Phase 2A audit only; no covariate model run."},
    ]
    pd.DataFrame(sample_rows).to_csv(OUTPUT_DIR / "mbzp_crc_phase2_sample_audit.csv", index=False)

    run_time = datetime.now(timezone.utc).isoformat()
    manifest = {
        "analysis": "MBzP-CRC Phase 2A NHANES sample audit",
        "run_timestamp_utc": run_time,
        "audit_only": True,
        "regression_performed": False,
        "cycles": [spec[0] for spec in CYCLE_SPECS],
        "data_source": "CDC NHANES public XPT files",
        "source_file_map": source_files,
        "laboratory_measure": "URXMZP (mono-benzyl phthalate; ng/mL)",
        "outcome_variable": "MCQ220 cancer history availability",
        "crc_type_codes": CRC_TYPE_CODES,
        "crc_definition": "Union of MCQ230A-D code 16 (Colon) and code 31 (Rectum)",
        "counts_are": "unweighted participant counts",
        "totals": totals,
        "weight_columns_recorded": cycle_df[["cycle", "mbzp_weight_column"]].to_dict(orient="records"),
        "software": {"python": platform.python_version(), "pandas": pd.__version__},
        "next_gate": "Do not run Phase 2B regression until the event-count audit is accepted.",
    }
    (OUTPUT_DIR / "mbzp_crc_phase2_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    report = f"""# MBzP–CRC Phase 2A NHANES sample audit

Run timestamp (UTC): `{run_time}`  
Scope: sample/event audit only; no regression was performed.

## Primary audit results

- MBzP + CRC outcome available: **{totals['mbzp_crc_outcome_n']:,}** participants.
- Combined CRC cases: **{totals['crc_cases_n']:,}** unweighted cases.
- Colon cases: {totals['colon_cases_n']:,}; rectal cases: {totals['rectal_cases_n']:,}; both type codes: {totals['both_colon_rectal_cases_n']:,}.

The combined CRC count is the union of MCQ230A-D code 16 (Colon) and code 31 (Rectum), so colon and rectal subtype counts are not added to obtain the primary event count.

## Cycle-level CRC cases

| Cycle | MBzP + outcome N | CRC cases |
|---|---:|---:|
""" + "\n".join(f"| {r['cycle']} | {r['mbzp_crc_outcome_n']:,} | {r['crc_cases_n']:,} |" for r in cycle_rows) + """

All values are unweighted counts. NHANES files and sampling weights are retained locally and are not committed; the manifest records the exact file map, hashes, variables, and weight columns. Phase 2B models remain deferred pending acceptance of this audit.
"""
    (OUTPUT_DIR / "mbzp_crc_phase2a_audit_report.md").write_text(report, encoding="utf-8")

    print(cycle_df[["cycle", "mbzp_crc_outcome_n", "crc_cases_n"]].to_string(index=False))
    print("TOTAL", totals)


if __name__ == "__main__":
    main()
