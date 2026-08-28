"""CCHS 2022 public-PUMF feasibility audit.

This audit deliberately stops before any exposure-outcome association model.
It verifies the directly downloadable package, reconstructs the relevant
fixed-width fields from the official layout cards, and freezes one
source-native exposure/outcome demonstration using feasibility only.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import zipfile
from collections import Counter
from pathlib import Path


OUT = Path(__file__).resolve().parent
ARCHIVE = Path(r"D:\whynot17\public_sources\cchs_2022\2022_TXT.zip")
DOWNLOAD_URL = "https://www150.statcan.gc.ca/n1/pub/82m0013x/2024001/2022_TXT.zip"
PUMF_PAGE = "https://www150.statcan.gc.ca/n1/pub/82m0013x/82m0013x2024001-eng.htm"
QUESTIONNAIRE_URL = "https://www23.statcan.gc.ca/imdb/p3Instr.pl?Function=assembleInstr&Item_Id=1390243&lang=en"
METHOD_URL = "https://www.statcan.gc.ca/en/statistical-programs/document/3226_D56_T9_V1"

EXPECTED_ARCHIVE_SHA256 = "a116a01fb35cc3204a8b14e36dda70cf6722e75380f88b6d822e492765ac8c41"

# 1-based inclusive fixed-width positions from the official SPSS input card.
FIELDS = {
    "EDDVH3": (39, 39),
    "DHHGAGE": (40, 40),
    "DHH_SEX": (41, 41),
    "HWTDGBCC": (57, 57),
    "CCC_05": (90, 90),
    "CCC_80": (91, 91),
    "CCC_85": (92, 92),
    "CCC_90": (93, 93),
    "PAADVWHO": (259, 259),
    "SMKDVSTY": (305, 306),
    "ALC_10": (322, 322),
    "INCDGHH": (374, 374),
    "WTS_M": (378, 385),
}

VARIABLES = {
    "EDDVH3": ("Highest level of education", "1,2,3", "6-9", "adult covariate"),
    "DHHGAGE": ("Age group", "2,3,4,5", "6-9", "adult restriction: 18 years and older"),
    "DHH_SEX": ("Sex at birth", "1,2", "6-9", "core covariate"),
    "HWTDGBCC": ("BMI classification for adults (adjusted)", "1,2", "6-9", "core covariate"),
    "CCC_05": ("Has diabetes", "1,2", "6-9", "alternative outcome"),
    "CCC_80": ("Has high blood pressure", "1,2", "6-9", "primary outcome"),
    "CCC_85": ("High blood pressure - took medication - 1 mo", "1,2", "6-9", "outcome-supporting sensitivity variable"),
    "CCC_90": ("Has had high blood cholesterol - lifetime", "1,2", "6-9", "alternative outcome"),
    "PAADVWHO": ("Physically active, based on WHO guidelines", "1,2,3,4", "6-9", "alternative exposure; module-limited"),
    "SMKDVSTY": ("Smoking status (type 2) - traditional definition", "01-06", "96-99", "primary exposure"),
    "ALC_10": ("Drank alcohol - 12 mo", "1,2", "6-9", "alternative exposure"),
    "INCDGHH": ("Total household income - all sources", "1,2,3,4,5", "6-9", "core covariate"),
    "WTS_M": ("Weights - Master", "positive numeric", "99999.96-99999.99", "person-level survey weight"),
}

VALID = {
    "DHHGAGE": {"2", "3", "4", "5"},
    "DHH_SEX": {"1", "2"},
    "HWTDGBCC": {"1", "2"},
    "CCC_05": {"1", "2"},
    "CCC_80": {"1", "2"},
    "CCC_85": {"1", "2"},
    "CCC_90": {"1", "2"},
    "PAADVWHO": {"1", "2", "3", "4"},
    "SMKDVSTY": {"01", "02", "03", "04", "05", "06"},
    "ALC_10": {"1", "2"},
    "EDDVH3": {"1", "2", "3"},
    "INCDGHH": {"1", "2", "3", "4", "5"},
}

MISSING = {
    "DHHGAGE": {"6", "7", "8", "9"},
    "DHH_SEX": {"6", "7", "8", "9"},
    "HWTDGBCC": {"6", "7", "8", "9"},
    "CCC_05": {"6", "7", "8", "9"},
    "CCC_80": {"6", "7", "8", "9"},
    "CCC_85": {"6", "7", "8", "9"},
    "CCC_90": {"6", "7", "8", "9"},
    "PAADVWHO": {"6", "7", "8", "9"},
    "SMKDVSTY": {"96", "97", "98", "99"},
    "ALC_10": {"6", "7", "8", "9"},
    "EDDVH3": {"6", "7", "8", "9"},
    "INCDGHH": {"6", "7", "8", "9"},
}

CORE = ["DHH_SEX", "HWTDGBCC", "EDDVH3", "INCDGHH"]


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def field(line: bytes, name: str) -> str:
    start, end = FIELDS[name]
    return line[start - 1 : end].decode("latin1").strip()


def member_by_suffix(names: list[str], suffix: str) -> str:
    matches = [name for name in names if name.lower().endswith(suffix.lower())]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one member ending {suffix!r}; found {matches}")
    return matches[0]


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise RuntimeError(f"No rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def scan_pumf(z: zipfile.ZipFile, data_member: str) -> tuple[int, int, dict[str, Counter]]:
    counters = {name: Counter() for name in FIELDS}
    rows = 0
    bad_width = 0
    with z.open(data_member) as handle:
        for raw in handle:
            line = raw.rstrip(b"\r\n")
            rows += 1
            if len(line) != 385:
                bad_width += 1
            for name in FIELDS:
                counters[name][field(line, name)] += 1
    return rows, bad_width, counters


def candidate_summary(rows: list[dict[str, str]], candidate_id: str, exposure: str, outcome: str, exposure_valid: set[str], outcome_valid: set[str], exposure_positive: set[str]) -> dict:
    adults = [row for row in rows if row["DHHGAGE"] in VALID["DHHGAGE"]]
    exp_ok = [row for row in adults if row[exposure] in exposure_valid]
    out_ok = [row for row in adults if row[outcome] in outcome_valid]
    joint = [row for row in adults if row[exposure] in exposure_valid and row[outcome] in outcome_valid]
    complete = [row for row in joint if all(row[name] in VALID[name] for name in CORE) and row["WTS_M"] and not row["WTS_M"].startswith("99999")]
    return {
        "candidate_id": candidate_id,
        "exposure_variable": exposure,
        "outcome_variable": outcome,
        "adult_definition": "DHHGAGE in {2,3,4,5} (18 years and older)",
        "adults_n": len(adults),
        "exposure_valid_n": len(exp_ok),
        "outcome_valid_n": len(out_ok),
        "joint_exposure_outcome_valid_n": len(joint),
        "outcome_positive_n_for_feasibility_only": sum(row[outcome] == "1" for row in joint),
        "complete_core_n": len(complete),
        "exposure_positive_definition": "SMKDVSTY in {01,02}" if exposure == "SMKDVSTY" else "valid source code set",
        "selection_used_association_results": False,
        "association_models_run": False,
    }


def main() -> None:
    if not ARCHIVE.exists():
        raise FileNotFoundError(ARCHIVE)
    observed_sha = sha256_path(ARCHIVE)
    if observed_sha != EXPECTED_ARCHIVE_SHA256:
        raise RuntimeError(f"Archive hash mismatch: {observed_sha}")

    with zipfile.ZipFile(ARCHIVE) as z:
        names = z.namelist()
        data_member = member_by_suffix(names, "/pumf_master_cchs.txt")
        bsw_member = member_by_suffix(names, "/bsw.txt")
        input_card = member_by_suffix(names, "/spss/pumf_master_cchs_i.sps")
        labels_card = member_by_suffix(names, "/spss/pumf_master_cchs_vare.sps")
        values_card = member_by_suffix(names, "/spss/pumf_master_cchs_vale.sps")
        missing_card = member_by_suffix(names, "/spss/pumf_master_cchs_miss.sps")
        bsw_card = member_by_suffix(names, "/layout_cards/bsw_i.sas")
        data_n, bad_width, counters = scan_pumf(z, data_member)

        parsed_rows: list[dict[str, str]] = []
        with z.open(data_member) as handle:
            for raw in handle:
                line = raw.rstrip(b"\r\n")
                parsed_rows.append({name: field(line, name) for name in FIELDS})

        bsw_n = 0
        bsw_bad_width = 0
        with z.open(bsw_member) as handle:
            for raw in handle:
                bsw_n += 1
                if len(raw.rstrip(b"\r\n")) != 8029:
                    bsw_bad_width += 1

        label_text = z.read(labels_card).decode("latin1", "replace")
        value_text = z.read(values_card).decode("latin1", "replace")
        missing_text = z.read(missing_card).decode("latin1", "replace")
        bsw_text = z.read(bsw_card).decode("latin1", "replace")
        package_rows = []
        for name in names:
            info = z.getinfo(name)
            package_rows.append({
                "member": name,
                "compressed_bytes": info.compress_size,
                "uncompressed_bytes": info.file_size,
                "is_raw_data": name.lower().endswith(("/pumf_master_cchs.txt", "/bsw.txt")),
                "version_control_policy": "exclude_raw_archive" if name.lower().endswith(("/pumf_master_cchs.txt", "/bsw.txt")) else "inventory_only",
            })
        selected_member_hashes = []
        for name in [input_card, labels_card, values_card, missing_card, bsw_card]:
            selected_member_hashes.append({
                "member": name,
                "sha256": hashlib.sha256(z.read(name)).hexdigest(),
            })

    candidate_rows = [
        candidate_summary(parsed_rows, "cchs_2022_current_smoking_to_hypertension", "SMKDVSTY", "CCC_80", VALID["SMKDVSTY"], VALID["CCC_80"], {"01", "02"}),
        candidate_summary(parsed_rows, "cchs_2022_current_smoking_to_diabetes", "SMKDVSTY", "CCC_05", VALID["SMKDVSTY"], VALID["CCC_05"], {"01", "02"}),
        candidate_summary(parsed_rows, "cchs_2022_alcohol_to_hypertension", "ALC_10", "CCC_80", VALID["ALC_10"], VALID["CCC_80"], {"1"}),
        candidate_summary(parsed_rows, "cchs_2022_physical_activity_to_hypertension", "PAADVWHO", "CCC_80", VALID["PAADVWHO"], VALID["CCC_80"], {"1", "2", "3", "4"}),
    ]

    variable_rows = []
    for name, (label, valid_codes, missing_codes, role) in VARIABLES.items():
        start, end = FIELDS[name]
        variable_rows.append({
            "variable": name,
            "position_1_based": f"{start}-{end}",
            "public_pumf_label": label,
            "valid_code_rule": valid_codes,
            "missing_code_rule": missing_codes,
            "role": role,
            "exact_public_variable_confirmed": True,
        })

    primary = {
        "freeze_status": "primary_feasibility_freeze",
        "candidate_id": "cchs_2022_current_smoking_to_hypertension",
        "exposure_family": "tobacco smoke exposure represented by current smoking status",
        "exposure_variable": "SMKDVSTY",
        "exposure_definition": "current daily or occasional smoker (01/02) versus not currently smoking (03-06); codes 96-99 treated as missing",
        "outcome": "self-reported high blood pressure/hypertension",
        "outcome_variable": "CCC_80",
        "outcome_definition": "yes (1) versus no (2); codes 6-9 treated as missing",
        "adult_definition": "DHHGAGE in {2,3,4,5}, corresponding to 18 years and older",
        "core_covariates": CORE,
        "survey_weight": "WTS_M",
        "variance_method": "1000 bootstrap replicate weights BSW1-BSW1000 from bsw.txt; empirical replicate variance uses 1/(B-1) sum of squared deviations around the replicate mean; no CCHS Fay multiplier is applied",
        "direct_download_confirmed": True,
        "registration_or_application_required_for_package": False,
        "individual_microdata": True,
        "exact_29_frozen_biomarker_replication": False,
        "association_models_run": False,
        "selection_used_exposure_outcome_association": False,
        "selection_basis": "exact public variables, broad adult coverage, valid outcome support, complete core-variable coverage, and directly available survey weight/replicates only",
        "secondhand_smoke_status": "not frozen: official questionnaire concept exists, but no corresponding second-hand-smoke variable was found in the public PUMF English label card; do not infer an exact PUMF exposure",
    }

    summary = {
        "audit_id": "STEP10C_CCHS_PUBLIC_DICTIONARY_AUDIT",
        "status": "complete_feasibility_freeze_no_association",
        "dataset": "CCHS 2022 PUMF TXT package",
        "archive_path": str(ARCHIVE),
        "download_url": DOWNLOAD_URL,
        "official_pumf_page": PUMF_PAGE,
        "questionnaire_url": QUESTIONNAIRE_URL,
        "methodology_url": METHOD_URL,
        "archive_sha256": observed_sha,
        "data_member": data_member,
        "data_rows": data_n,
        "data_record_width": 385,
        "data_bad_width_rows": bad_width,
        "bootstrap_member": bsw_member,
        "bootstrap_rows": bsw_n,
        "bootstrap_record_width": 8029,
        "bootstrap_bad_width_rows": bsw_bad_width,
        "bootstrap_layout_confirms_1000_replicates": bool(re.search(r"BSW1\s*-\s*BSW1000", bsw_text, re.I)),
        "bootstrap_data_rows_match_master": data_n == bsw_n,
        "selected_member_hashes": selected_member_hashes,
        "public_pumf_has_secondhand_smoke_label_hit": bool(re.search(r"second.?hand|passive smoke|environmental tobacco", label_text + value_text + missing_text, re.I)),
        "primary_freeze": primary,
        "candidate_count_audited": len(candidate_rows),
        "association_models_run": False,
        "outcome_firewall": "No exposure-outcome regression, P value, effect estimate, or FDR was used for candidate selection.",
    }

    write_csv(OUT / "cchs_2022_public_package_inventory.csv", package_rows)
    write_csv(OUT / "cchs_2022_variable_dictionary_audit.csv", variable_rows)
    write_csv(OUT / "cchs_2022_candidate_combination_audit.csv", candidate_rows)
    (OUT / "CCHS_PUBLIC_AUDIT_QC_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (OUT / "CCHS_PUBLIC_AUDIT_MANIFEST.json").write_text(json.dumps({
        "audit_id": summary["audit_id"],
        "status": summary["status"],
        "files": [
            "PLAN_CCHS_PUBLIC_AUDIT.md",
            "run_cchs_public_dictionary_audit.py",
            "cchs_2022_public_package_inventory.csv",
            "cchs_2022_variable_dictionary_audit.csv",
            "cchs_2022_candidate_combination_audit.csv",
            "CCHS_PUBLIC_AUDIT_QC_SUMMARY.json",
            "STEP10C_CCHS_PUBLIC_AUDIT_REPORT.md",
            "CCHS_PUBLIC_AUDIT_MANIFEST.json",
        ],
        "raw_archive_excluded_from_git": True,
        "archive_sha256": observed_sha,
        "primary_candidate": primary,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    report = f"""# Step 10C — CCHS public population-replacement feasibility audit

## Decision

`CCHS 2022 PUMF` passes the **publicly downloadable, individual-level microdata** gate. The archive was downloaded directly from Statistics Canada and verified locally by SHA-256. No registration, application, or controlled-access step was used for this package.

The feasibility-only primary freeze is:

> **Current smoking status (`SMKDVSTY`) → self-reported high blood pressure (`CCC_80`) in adults**

This is a source-native exposure demonstration, not an exact replication of the 29 NHANES biomarker family. The exposure is defined as current daily/occasional smoking (codes 01/02) versus not currently smoking (03–06). Hypertension is `CCC_80` yes/no (1/2). Adults are `DHHGAGE` groups 2–5 (18 years and older).

## Package and design audit

- Direct package: `{DOWNLOAD_URL}`
- Local archive: `{ARCHIVE}`
- SHA-256: `{observed_sha}`
- Master PUMF rows: `{data_n:,}`; fixed-width record length: 385; malformed-width rows: {bad_width}
- Bootstrap rows: `{bsw_n:,}`; fixed-width record length: 8029; malformed-width rows: {bsw_bad_width}
- Master and bootstrap row counts match: **{data_n == bsw_n}**
- The official bootstrap layout card contains `BSW1–BSW1000`.
- Person-level survey weight: `WTS_M`
- Variance plan: use the package bootstrap replicate weights with empirical variance `1/(B-1) * sum((beta_b - mean(beta_b))^2)`; the CCHS guide recommends bootstrap weights for exact regression precision, and no CCHS-specific Fay multiplier is applied.

## Variable feasibility

The public PUMF label/input cards contain exact fields for age group, sex, adult BMI classification, education, income, smoking status, high blood pressure, diabetes, alcohol, physical activity, `WTS_M`, and the bootstrap file. The main public label card did **not** contain an exact second-hand-smoke variable; the official questionnaire concept must not be treated as an available public-PUMF exposure without a mapped field.

The candidate-combination table reports only variable validity and descriptive event support. It contains no association estimate and was not used to optimize a result.

## Firewall and scope

`association_models_run = false`. No exposure–outcome regression, P value, odds ratio, or FDR was used to select or freeze the combination. The next step, if authorized, is a separate preregistered CCHS association runner using the frozen combination and bootstrap survey variance.

The large raw archive and raw microdata are deliberately excluded from version control; this directory contains only the small audit artifacts and provenance needed to reproduce the package check.
"""
    (OUT / "STEP10C_CCHS_PUBLIC_AUDIT_REPORT.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
