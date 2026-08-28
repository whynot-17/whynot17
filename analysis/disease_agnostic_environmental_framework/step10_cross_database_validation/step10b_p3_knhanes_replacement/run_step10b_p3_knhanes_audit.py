#!/usr/bin/env python
"""Outcome-free KNHANES population-replacement feasibility audit.

The official KNHANES raw-data catalogue is queried for metadata only.  No
microdata file is downloaded, no user registration is attempted, and no
exposure--diabetes association is estimated.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
FROZEN = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
RUN_UTC = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")

PORTAL_URL = "https://www.data.go.kr/tcs/dss/selectFileDataDetailView.do?publicDataPk=15076556"
KNHANES_MAIN_URL = "https://knhanes.kdca.go.kr/knhanes/main.do"
KNHANES_EN_URL = "https://knhanes.kdca.go.kr/knhanes/eng/main.do"
RAW_PAGE_URL = "https://knhanes.kdca.go.kr/knhanes/rawDataDwnld/rawDataDwnld.do"
RAW_API_URL = "https://knhanes.kdca.go.kr/knhanes/rawDataDwnld/findRawdtaList.json"
RAW_REDIRECT_URL = "https://knhanes.kdca.go.kr/knhanes/postSendPage.do?url=/rawDataDwnld/rawDataDwnld.do&postparam=%7B%22menuId%22:%2210031001%22%7D"

SOURCES = {
    "KDCA_OFFICIAL_EN": {
        "title": "KDCA National Health and Nutrition Survey overview",
        "url": KNHANES_EN_URL,
        "source_type": "official_program_page",
    },
    "KDCA_OFFICIAL_KO": {
        "title": "KNHANES official portal",
        "url": KNHANES_MAIN_URL,
        "source_type": "official_data_portal",
    },
    "RAW_DATA_CATALOGUE": {
        "title": "KNHANES raw-data catalogue page",
        "url": RAW_PAGE_URL,
        "source_type": "official_raw_data_catalogue",
    },
    "DATA_GO_KR": {
        "title": "KDCA raw-data user guide record on Korea public data portal",
        "url": PORTAL_URL,
        "source_type": "official_public_data_record",
    },
    "ARSENIC_DIABETES": {
        "title": "KNHANES 2008–2009 urinary arsenic and diabetes feasibility precedent",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3678002/",
        "source_type": "peer_reviewed_feasibility_precedent",
    },
    "BLOOD_METALS": {
        "title": "KNHANES blood lead, cadmium and mercury exposure precedent",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5743444/",
        "source_type": "peer_reviewed_exposure_precedent",
    },
    "DIABETES_TRENDS": {
        "title": "KNHANES diabetes definition and trend precedent",
        "url": "https://diabetesjournals.org/care/article/32/11/2016/25967/Prevalence-and-Management-of-Diabetes-in-Korean",
        "source_type": "peer_reviewed_outcome_precedent",
    },
}

KNHANES_YEARS = [str(year) for year in range(2007, 2025)]


# KNHANES papers and the official portal document a different environmental
# sampling programme from KoNEHS.  The only frozen-test near matches retained
# here are the blood-metal family precedents; they are explicitly matrix
# mismatches for the frozen urine tests.  Frozen serum PFAS and urinary
# phthalate/PAH/metal tests are not upgraded without direct KNHANES evidence.
EXPOSURE_EVIDENCE = {
    "LBXPFDE": dict(domain="PFAS", analyte="PFDA/PFDeA", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES serum PFDA/PFDeA evidence confirmed in the audited sources."),
    "LBXPFHS": dict(domain="PFAS", analyte="PFHxS", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES serum PFHxS evidence confirmed in the audited sources."),
    "LBXPFNA": dict(domain="PFAS", analyte="PFNA", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES serum PFNA evidence confirmed in the audited sources."),
    "URXBPH": dict(domain="bisphenol", analyte="BPA", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary BPA evidence confirmed in the audited sources."),
    "URXCOP": dict(domain="phthalate", analyte="MCOP/MCiOP", match="not_confirmed_public", years=[], basis="", note="MCOP exposure evidence was found for KoNEHS, not KNHANES; it is not transferred across surveys."),
    "URXECP": dict(domain="phthalate", analyte="MECPP", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary MECPP evidence confirmed in the audited sources."),
    "URXMBP": dict(domain="phthalate", analyte="MnBP", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary MnBP evidence confirmed in the audited sources."),
    "URXMEP": dict(domain="phthalate", analyte="MEP", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary MEP evidence confirmed in the audited sources."),
    "URXMHH": dict(domain="phthalate", analyte="MEHHP", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary MEHHP evidence confirmed in the audited sources."),
    "URXMHP": dict(domain="phthalate", analyte="MEHP", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary MEHP evidence confirmed in the audited sources."),
    "URXMIB": dict(domain="phthalate", analyte="MiBP", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary MiBP evidence confirmed in the audited sources."),
    "URXMOH": dict(domain="phthalate", analyte="MEOHP", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary MEOHP evidence confirmed in the audited sources."),
    "URXMZP": dict(domain="phthalate", analyte="MBzP", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary MBzP evidence confirmed in the audited sources."),
    "URXP02": dict(domain="PAH", analyte="naphthalene metabolite family", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary PAH family evidence confirmed in the audited sources."),
    "URXP04": dict(domain="PAH", analyte="fluorene metabolite family", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary PAH family evidence confirmed in the audited sources."),
    "URXP10": dict(domain="PAH", analyte="pyrene metabolite family", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary PAH family evidence confirmed in the audited sources."),
    "URXP25": dict(domain="PAH", analyte="phenanthrene metabolite family", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary PAH family evidence confirmed in the audited sources."),
    "URXUBA": dict(domain="metal", analyte="urine barium", match="not_confirmed_public", years=[], basis="", note="The audited KNHANES metal precedent did not confirm frozen urine barium."),
    "URXUCD": dict(domain="metal", analyte="urine cadmium", match="family_public_matrix_mismatch", years=["2008", "2009"], basis="BLOOD_METALS", note="KNHANES public precedent documents blood cadmium, not the frozen urine-cadmium test."),
    "URXUCO": dict(domain="metal", analyte="urine cobalt", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary cobalt evidence confirmed in the audited sources."),
    "URXUCS": dict(domain="metal", analyte="urine cesium", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary cesium evidence confirmed in the audited sources."),
    "URXUMO": dict(domain="metal", analyte="urine molybdenum", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary molybdenum evidence confirmed in the audited sources."),
    "URXUPB": dict(domain="metal", analyte="urine lead", match="family_public_matrix_mismatch", years=["2008", "2009"], basis="BLOOD_METALS", note="KNHANES public precedent documents blood lead, not the frozen urine-lead test."),
    "URXUSB": dict(domain="metal", analyte="urine antimony", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary antimony evidence confirmed in the audited sources."),
    "URXUSN": dict(domain="metal", analyte="urine tin", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary tin evidence confirmed in the audited sources."),
    "URXUSR": dict(domain="metal", analyte="urine silver", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary silver evidence confirmed in the audited sources."),
    "URXUTL": dict(domain="metal", analyte="urine thallium", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary thallium evidence confirmed in the audited sources."),
    "URXUTU": dict(domain="metal", analyte="urine tungsten", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary tungsten evidence confirmed in the audited sources."),
    "URXUUR": dict(domain="metal", analyte="urine uranium", match="not_confirmed_public", years=[], basis="", note="No exact KNHANES urinary uranium evidence confirmed in the audited sources."),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fetch_bytes(url: str, method: str = "GET", body: bytes | None = None, headers: dict[str, str] | None = None, limit: int = 2_000_000) -> dict[str, object]:
    entry: dict[str, object] = {"url": url, "method": method, "retrieval_status": "not_attempted"}
    try:
        req_headers = {"User-Agent": "whynot17-KNHANES-audit/1.0"}
        if headers:
            req_headers.update(headers)
        req = Request(url, data=body, headers=req_headers, method=method)
        with urlopen(req, timeout=12) as response:  # nosec B310: frozen public URLs
            data = response.read(limit + 1)
            entry.update({
                "http_status": getattr(response, "status", None),
                "retrieval_status": "ok" if len(data) <= limit else "truncated",
                "content_bytes_captured": min(len(data), limit),
                "sha256_captured": hashlib.sha256(data[:limit]).hexdigest(),
                "content_type": response.headers.get("Content-Type", ""),
            })
    except (HTTPError, URLError, TimeoutError, OSError) as exc:
        entry.update({"retrieval_status": "failed_best_effort", "error_type": type(exc).__name__, "error": str(exc)[:300]})
    return entry


def fetch_source_snapshot() -> dict[str, object]:
    snapshot: dict[str, object] = {"retrieved_utc": RUN_UTC, "sources": {}}
    for source_id, metadata in SOURCES.items():
        entry = dict(metadata)
        entry.update({k: v for k, v in fetch_bytes(metadata["url"]).items() if k not in {"url"}})
        snapshot["sources"][source_id] = entry
    return snapshot


def fetch_raw_catalog() -> tuple[list[dict[str, object]], dict[str, object]]:
    all_rows: list[dict[str, object]] = []
    probe: dict[str, object] = {
        "catalog_api_url": RAW_API_URL,
        "raw_page_url": RAW_PAGE_URL,
        "official_portal_url": PORTAL_URL,
        "api_method": "POST",
        "metadata_only": True,
        "microdata_downloaded": False,
        "user_registration_attempted": False,
        "file_download_status": "registration_required_and_not_completed_in_audit",
        "api_status": "not_attempted",
    }
    first_body = json.dumps({"bizBgngYr": "", "bizEndYr": "", "rawdtaDmnReqList": ["5", "6", "7", "8"], "pageNo": 1}).encode("utf-8")
    try:
        req = Request(RAW_API_URL, data=first_body, headers={"User-Agent": "whynot17-KNHANES-audit/1.0", "Content-Type": "application/json;charset=UTF-8"}, method="POST")
        with urlopen(req, timeout=20) as response:  # nosec B310: official public metadata API
            payload = json.loads(response.read(5_000_000).decode("utf-8"))
        data = payload.get("data", {})
        probe.update({"api_status": payload.get("status", "unknown"), "total_count": data.get("totCnt"), "page_count": data.get("totPageCnt")})
        all_rows.extend(data.get("rawdtaList", []))
        total_pages = int(data.get("totPageCnt") or 1)
        for page in range(2, total_pages + 1):
            body = json.dumps({"bizBgngYr": "", "bizEndYr": "", "rawdtaDmnReqList": ["5", "6", "7", "8"], "pageNo": page}).encode("utf-8")
            req = Request(RAW_API_URL, data=body, headers={"User-Agent": "whynot17-KNHANES-audit/1.0", "Content-Type": "application/json;charset=UTF-8"}, method="POST")
            with urlopen(req, timeout=20) as response:  # nosec B310: official public metadata API
                page_payload = json.loads(response.read(5_000_000).decode("utf-8"))
            all_rows.extend(page_payload.get("data", {}).get("rawdtaList", []))
    except (HTTPError, URLError, TimeoutError, OSError, ValueError, KeyError) as exc:
        probe.update({"api_status": "failed_best_effort", "error_type": type(exc).__name__, "error": str(exc)[:300]})
    return all_rows, probe


def build_crosswalk(frozen: list[dict[str, str]]) -> list[dict[str, object]]:
    rows = []
    for item in frozen:
        variable = item["variable"]
        evidence = EXPOSURE_EVIDENCE[variable]
        rows.append({
            "test_id": item["test_id"],
            "frozen_biomarker": item["biomarker"],
            "nhanes_variable": variable,
            "nhanes_matrix": item["matrix"],
            "frozen_nhanes_cycles": item["cycles"],
            "domain": evidence["domain"],
            "knhanes_candidate_analyte": evidence["analyte"],
            "knhanes_match_class": evidence["match"],
            "knhanes_years_with_public_evidence": ";".join(evidence["years"]),
            "evidence_basis": evidence["basis"] or "none_in_audited_public_sources",
            "evidence_sources": evidence["basis"] or "",
            "exact_analyte_confirmed": "no",
            "exact_matrix_confirmed": "no",
            "raw_catalog_record_for_exact_test": "not_found",
            "diabetes_joint_extract_confirmed": "no",
            "individual_access_status": "registration_required_and_not_completed_in_audit",
            "notes": evidence["note"],
        })
    return rows


def build_year_readiness(crosswalk: list[dict[str, object]], catalog: list[dict[str, object]]) -> list[dict[str, object]]:
    env_records = [row for row in catalog if str(row.get("clsf", "")).find("환경") >= 0 or str(row.get("tpcCn", "")).find("환경") >= 0 or str(row.get("clsf", "")).find("유해") >= 0]
    env_summary = ";".join(sorted({f"{row.get('bizYr','')}:{row.get('tpcCn','')}" for row in env_records}))
    rows = []
    for item in crosswalk:
        for year in KNHANES_YEARS:
            evidence_year = year in str(item["knhanes_years_with_public_evidence"]).split(";")
            rows.append({
                "test_id": item["test_id"],
                "nhanes_variable": item["nhanes_variable"],
                "knhanes_year": year,
                "public_exact_or_family_evidence": "yes" if evidence_year else "no_or_not_confirmed",
                "match_class": item["knhanes_match_class"] if evidence_year else "not_confirmed_for_year",
                "exact_variable_confirmed": "no",
                "exact_matrix_confirmed": "no",
                "raw_catalog_exact_record": "not_found",
                "environmental_catalog_metadata_context": env_summary,
                "joint_diabetes_confirmed": "no",
                "weights_design_confirmed": "no",
            })
    return rows


def build_outcome_audit() -> list[dict[str, object]]:
    specs = [
        ("raw_data_download_path", "access", "public_registration_path", "Official KDCA/data.go.kr pages expose raw-data guidance and a raw-data catalogue; the portal presents registration/consent before file download."),
        ("microdata_file_download", "access", "not_executed", "No file was downloaded and no personal information was entered in this audit."),
        ("diabetes_health_statistic", "outcome", "documented", "KDCA official survey overview lists diabetes among health-examination/chronic-disease content."),
        ("fasting_glucose", "outcome", "documented_by_published_precedent", "KNHANES 2008–2009 diabetes/arsenic paper used fasting plasma glucose in the analytic workflow."),
        ("HbA1c", "outcome", "documented_by_published_precedent", "KNHANES diabetes studies use glycemic measures; exact current field/version requires codebook confirmation."),
        ("physician_diagnosis_medication", "outcome", "documented_by_published_precedent", "Published KNHANES diabetes analyses use diagnosis and/or diabetes treatment information."),
        ("type1_exclusion", "outcome", "conditional", "A reproducible T2D-specific restriction must be confirmed from the current questionnaire/medication fields."),
        ("age_sex", "covariate", "documented", "Official survey scope and published studies include demographic variables."),
        ("BMI_anthropometry", "covariate", "documented", "Health examination content includes obesity/anthropometric measures."),
        ("smoking_alcohol_activity", "covariate", "documented", "Official survey overview lists smoking, drinking and physical activity questionnaire domains."),
        ("SES_education_income", "covariate", "documented_by_published_precedent", "Published exposure–diabetes analyses use education/income or residence-related covariates."),
        ("urinary_creatinine", "laboratory_qc", "conditional", "Urine-specific creatinine and laboratory flags must be confirmed for any exact exposure file; not relevant to blood-only near matches."),
        ("survey_weights", "survey_design", "documented_methodology_names_pending", "KNHANES is a probability sample and published analyses use survey procedures; exact current weight variable names require the codebook."),
        ("strata_psu", "survey_design", "documented_methodology_names_pending", "Complex stratified/multistage design is documented in published KNHANES analyses; exact public-use fields require codebook confirmation."),
        ("environmental_lab_panel", "exposure", "partial", "The public catalogue exposes a 2021 HN_IAQ environmental/indoor-air-related module, but it does not establish coverage of the frozen 29 biomarkers."),
        ("published_external_precedent", "overall_feasibility", "yes_but_not_replication", "KNHANES publications demonstrate linkage of biomarker, diabetes/glycemia, covariates and survey design; no published estimate is imported into this audit."),
    ]
    rows = []
    for key, domain, status, note in specs:
        rows.append({
            "audit_item": key,
            "domain": domain,
            "status": status,
            "evidence_sources": "KDCA_OFFICIAL_EN;KDCA_OFFICIAL_KO;RAW_DATA_CATALOGUE;DATA_GO_KR;ARSENIC_DIABETES;BLOOD_METALS;DIABETES_TRENDS",
            "required_next_action": "After official registration, download the relevant raw database and verify exact fields, matrix, missingness and design variables." if status not in {"documented", "documented_by_published_precedent", "yes_but_not_replication", "partial"} else "Retain as feasibility evidence; verify exact coding in the downloaded codebook.",
            "notes": note,
        })
    return rows


def main() -> int:
    frozen = read_csv(FROZEN)
    if len(frozen) != 29:
        raise ValueError(f"Expected 29 frozen tests, got {len(frozen)}")
    missing = sorted(set(row["variable"] for row in frozen) - set(EXPOSURE_EVIDENCE))
    if missing:
        raise ValueError(f"Missing exposure mapping(s): {missing}")

    catalog, access_probe = fetch_raw_catalog()
    crosswalk = build_crosswalk(frozen)
    year_rows = build_year_readiness(crosswalk, catalog)
    outcome_rows = build_outcome_audit()
    source_snapshot = fetch_source_snapshot()

    catalog_rows = []
    for row in catalog:
        catalog_rows.append({
            "bizYr": row.get("bizYr", ""),
            "domain": row.get("dmn", ""),
            "domain_name": row.get("dmnNm", ""),
            "survey_area": row.get("clsf", ""),
            "database_topic": row.get("tpcCn", ""),
            "last_modified": row.get("mdfcnDt", ""),
            "sas_record_id": row.get("sasRawdtaSeq", ""),
            "spss_record_id": row.get("spssRawdtaSeq", ""),
            "catalog_metadata_only": "yes",
        })
    if not catalog_rows:
        catalog_rows = [{"bizYr": "", "domain": "", "domain_name": "", "survey_area": "", "database_topic": "", "last_modified": "", "sas_record_id": "", "spss_record_id": "", "catalog_metadata_only": "yes", "catalog_status": "no_records_returned"}]

    write_csv(OUT / "knhanes_29_test_crosswalk.csv", crosswalk)
    write_csv(OUT / "knhanes_year_readiness.csv", year_rows)
    write_csv(OUT / "knhanes_outcome_covariate_design_audit.csv", outcome_rows)
    write_csv(OUT / "knhanes_raw_data_catalog.csv", catalog_rows)
    (OUT / "KNHANES_ACCESS_PROBE.json").write_text(json.dumps(access_probe, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (OUT / "KNHANES_SOURCE_SNAPSHOT.json").write_text(json.dumps(source_snapshot, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    exact = sum(row["knhanes_match_class"] == "exact_public" for row in crosswalk)
    mismatch = sum(row["knhanes_match_class"] == "family_public_matrix_mismatch" for row in crosswalk)
    not_confirmed = sum(row["knhanes_match_class"] == "not_confirmed_public" for row in crosswalk)
    qc = {
        "generated_utc": RUN_UTC,
        "status": "catalog_access_confirmed_exact_29_test_crosswalk_not_confirmed",
        "decision": "not_currently_actionable_for_frozen_29_exact_replication",
        "outcome_free_crosswalk": True,
        "association_models_run": False,
        "published_effect_estimates_imported": False,
        "microdata_downloaded": False,
        "user_registration_attempted": False,
        "frozen_test_count_expected": 29,
        "frozen_test_count_observed": len(frozen),
        "exact_public_match_count": exact,
        "family_public_matrix_mismatch_count": mismatch,
        "not_confirmed_public_count": not_confirmed,
        "exact_match_fraction": round(exact / 29, 6),
        "raw_catalog_record_count": len(catalog),
        "raw_catalog_api_status": access_probe.get("api_status"),
        "raw_catalog_total_count": access_probe.get("total_count"),
        "raw_catalog_page_count": access_probe.get("page_count"),
        "diabetes_outcome_feasibility": "documented_by_official_scope_and_published_precedent_exact_t2d_coding_pending",
        "survey_design_feasibility": "documented_by_published_precedent_exact_weight_strata_psu_names_pending",
        "access_status": "public_registration_path_confirmed_file_download_not_completed",
        "limitations": [
            "The two related metal precedents are blood-matrix evidence and do not satisfy frozen urine tests.",
            "No exact frozen biomarker/matrix match was confirmed in the audited KNHANES sources.",
            "The raw-data catalogue is accessible as metadata, but the file download requires the official user-information/consent workflow.",
            "A published diabetes analysis is feasibility precedent, not this project's external replication.",
        ],
    }
    (OUT / "KNHANES_QC_SUMMARY.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    report = f"""# Step 10B-P3 — KNHANES executable population-replacement audit

Generated: `{RUN_UTC}`  \\
Status: **{qc['status']}**  \\
Decision: **{qc['decision']}**

## Executive result

KNHANES has a real official raw-data access route, but the present audit does **not** establish an exact match to any of the 29 frozen biomarker tests. The crosswalk contains **29/29 tests**, with **{exact}/29 exact public matches**, **{mismatch}/29 related-family blood-matrix mismatches**, and **{not_confirmed}/29 not confirmed in audited public sources**.

This is not a failed data-access test. The official KDCA portal and Korea public-data record expose a raw-data catalogue and a registration/consent workflow. The metadata API returned a catalogue with `{len(catalog)}` records (the API-reported total was `{access_probe.get('total_count')}` when available). However, this audit did not enter personal information or download a raw file, so file-level access is recorded as `public_registration_path_confirmed_file_download_not_completed`.

## Exact crosswalk result

The audited KNHANES literature documents blood lead, blood cadmium and blood mercury, and urinary total arsenic in specific 2008–2009 analyses. These are valuable environmental-health precedents, but the frozen panel contains serum PFAS and urine phthalate/PAH/metal tests. Blood lead/cadmium were therefore retained as explicit **matrix mismatches**, not upgraded to exact matches. KoNEHS phthalate/PFAS/MCOP evidence was not transferred to KNHANES.

No exact KNHANES evidence was confirmed for frozen serum PFDA/PFDeA, PFHxS, PFNA; urinary BPA, phthalate metabolites, PAH families; or the frozen urinary metals. The row-level evidence and notes are in `knhanes_29_test_crosswalk.csv`; annual metadata status is in `knhanes_year_readiness.csv`.

## Outcome, covariate and design feasibility

The official KDCA survey overview lists diabetes among chronic-disease/health-examination content and includes smoking, drinking, physical activity, obesity and other health domains. Published KNHANES studies demonstrate use of fasting glucose/glycemia, diabetes ascertainment, covariates and complex probability-sample methods. Exact current field names, type-1 handling, missingness, weights, strata and PSU variables remain codebook-level items to verify after download.

Thus the outcome/design side is **feasible by precedent**, but it cannot rescue a missing exact exposure matrix. The official raw catalogue also lists an environmental/indoor-air-related module in 2021; that catalogue metadata does not demonstrate coverage of the frozen 29 tests.

## What this means for population replacement

KNHANES is preferable to a purely hypothetical source because its official raw-data route is demonstrable and registration-based access is available. Nevertheless, under the frozen 29-test rule it is currently **not an exact biomarker population-replacement dataset**. It can be promoted only if the authorized download/codebook reveals additional environmental laboratory files containing one or more frozen analytes with the required matrix and a same-person diabetes frame.

## Frozen exclusions

- No association model was run.
- No published PFAS, arsenic, metal or diabetes effect estimate was adopted as our replication.
- No candidate was selected because it had a prior Korean paper.
- No KoNEHS measurement was treated as KNHANES evidence.
- No raw microdata were downloaded.

## Files

- `knhanes_29_test_crosswalk.csv`
- `knhanes_year_readiness.csv`
- `knhanes_outcome_covariate_design_audit.csv`
- `knhanes_raw_data_catalog.csv`
- `KNHANES_ACCESS_PROBE.json`
- `KNHANES_SOURCE_SNAPSHOT.json`
- `KNHANES_QC_SUMMARY.json`
- `STEP10B_P3_KNHANES_MANIFEST.json`
"""
    (OUT / "STEP10B_P3_KNHANES_AUDIT_REPORT.md").write_text(report, encoding="utf-8")

    output_names = [
        "PLAN.md",
        "run_step10b_p3_knhanes_audit.py",
        "knhanes_29_test_crosswalk.csv",
        "knhanes_year_readiness.csv",
        "knhanes_outcome_covariate_design_audit.csv",
        "knhanes_raw_data_catalog.csv",
        "KNHANES_ACCESS_PROBE.json",
        "KNHANES_SOURCE_SNAPSHOT.json",
        "KNHANES_QC_SUMMARY.json",
        "STEP10B_P3_KNHANES_AUDIT_REPORT.md",
    ]
    manifest = {
        "generated_utc": RUN_UTC,
        "script": "run_step10b_p3_knhanes_audit.py",
        "frozen_input": str(FROZEN.relative_to(ROOT)).replace("\\", "/"),
        "frozen_test_count": len(frozen),
        "outcome_free": True,
        "association_models_run": False,
        "published_effect_estimates_imported": False,
        "microdata_downloaded": False,
        "user_registration_attempted": False,
        "official_catalog_api": RAW_API_URL,
        "outputs": [],
    }
    for name in output_names:
        data = (OUT / name).read_bytes()
        manifest["outputs"].append({"path": name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
    (OUT / "STEP10B_P3_KNHANES_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
