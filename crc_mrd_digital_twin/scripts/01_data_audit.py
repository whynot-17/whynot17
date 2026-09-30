"""Phase 1 audit of Chen et al. CRC postoperative serial ctDNA supplements.

Reads original files without treating an unobserved draw as a negative result.
This script writes audit summaries only. It does not create a model-ready cohort
or fit survival models.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from docx import Document

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data_raw"
TABLES = ROOT / "outputs" / "tables"
QC = ROOT / "outputs" / "qc"
REPORTS = ROOT / "reports"
for folder in (TABLES, QC, REPORTS):
    folder.mkdir(parents=True, exist_ok=True)


def write_csv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(TABLES / name, index=False, encoding="utf-8-sig")


def read_xlsx(name: str, header: int) -> pd.DataFrame:
    return pd.read_excel(RAW / name, header=header, engine="openpyxl").dropna(how="all")


def read_baseline() -> pd.DataFrame:
    table = Document(RAW / "13045_2021_1089_MOESM1_ESM.docx").tables[0]
    headers = [cell.text.strip() for cell in table.rows[1].cells]
    rows = [[cell.text.strip() or pd.NA for cell in row.cells] for row in table.rows[2:]]
    frame = pd.DataFrame(rows, columns=headers).rename(columns={
        "Patient No.": "Patient_No", "PTLa": "PTL", "Hist_typeb": "Hist_type",
        "Hist_gradec": "Hist_grade", "LV_invasiond": "LV_invasion",
        "N_invasione": "N_invasion", "ACTf": "ACT", "R_timeg": "R_time",
        "R_treatmenth": "R_treatment",
    })
    frame["Age"] = pd.to_numeric(frame["Age"], errors="coerce")
    frame["R_time"] = pd.to_numeric(frame["R_time"], errors="coerce")
    frame["recurrence_event"] = frame["R_time"].notna()
    frame["ACT_yes"] = frame["ACT"].astype(str).str.lower().eq("yes")
    return frame


def quantiles(values: pd.Series) -> dict:
    x = pd.to_numeric(values, errors="coerce").dropna()
    if not len(x):
        return {"n": 0, "median": None, "q1": None, "q3": None, "min": None, "max": None}
    return {"n": int(len(x)), "median": float(x.median()), "q1": float(x.quantile(.25)),
            "q3": float(x.quantile(.75)), "min": float(x.min()), "max": float(x.max())}


def verify_downloads() -> tuple[int, int]:
    manifest = json.loads((RAW / "figshare_manifest.json").read_text(encoding="utf-8-sig"))
    verified = 0
    for item in manifest:
        path = RAW / item["filename"]
        digest = hashlib.md5(path.read_bytes()).hexdigest()
        if digest != item["md5"].lower():
            raise ValueError(f"MD5 mismatch: {path.name}")
        verified += 1
    return len(manifest), verified


n_manifest, n_verified = verify_downloads()
base = read_baseline()
variants = read_xlsx("13045_2021_1089_MOESM6_ESM.xlsx", 12)
ctdna = read_xlsx("13045_2021_1089_MOESM7_ESM.xlsx", 5)
cea = read_xlsx("13045_2021_1089_MOESM9_ESM.xlsx", 4)

ctdna["Time_postop"] = pd.to_numeric(ctdna["Time_postop"], errors="coerce")
ctdna["sample_rank"] = pd.to_numeric(ctdna["Sample_point"].astype(str).str.extract(r"P(\d+)")[0])
ctdna["positive"] = ctdna["ctDNA_status"].astype(str).str.lower().eq("positive")
cea["Time_postop"] = pd.to_numeric(cea["Time_postop"], errors="coerce")
cea["CEA_value"] = pd.to_numeric(cea["CEA_value"], errors="coerce")
cea_pre = cea[cea["Time_postop"] <= 0]
cea_pre_ids = set(cea_pre["Patient_No"])
patient_ids = set(base["Patient_No"])
recurrence_map = base.set_index("Patient_No")["R_time"].to_dict()
cea_ids = set(cea["Patient_No"])
ctdna["R_time"] = ctdna["Patient_No"].map(recurrence_map)
ctdna["pre_event"] = ctdna["R_time"].isna() | (ctdna["Time_postop"] < ctdna["R_time"])
ctdna["at_or_after_recurrence"] = ctdna["R_time"].notna() & (ctdna["Time_postop"] >= ctdna["R_time"])

# Patient-level draw counts are an audit output, not a prediction table.
draws = ctdna.groupby("Patient_No").agg(
    n_ctDNA_draws=("ctDNA_status", "size"),
    n_postoperative_draws=("sample_rank", lambda s: int((s >= 2).sum())),
    n_positive_draws=("positive", "sum"),
).reset_index()
draws = base[["Patient_No", "recurrence_event", "R_time", "ACT"]].merge(draws, how="left", on="Patient_No")
write_csv(draws, "patient_draw_counts.csv")
draws["outcome_group"] = np.where(draws["recurrence_event"], "Documented recurrence", "No documented recurrence")
draw_group_rows = []
for label, group in draws.groupby("outcome_group"):
    draw_group_rows.append({"outcome_group": label, **quantiles(group["n_ctDNA_draws"])})
write_csv(pd.DataFrame(draw_group_rows), "draw_counts_by_recurrence.csv")

# Counts, status rates and observed timing by protocol point; both all records and
# pre-recurrence records are shown so post-event samples remain visible as leakage risk.
points = []
for point, group in ctdna.groupby("Sample_point"):
    pre = group[group["pre_event"]]
    points.append({
        "sample_point": point,
        "n_samples": len(group),
        "n_patients": group["Patient_No"].nunique(),
        "n_positive": int(group["positive"].sum()),
        "positive_rate": float(group["positive"].mean()),
        "median_time_months": float(group["Time_postop"].median()),
        "q1_time_months": float(group["Time_postop"].quantile(.25)),
        "q3_time_months": float(group["Time_postop"].quantile(.75)),
        "min_time_months": float(group["Time_postop"].min()),
        "max_time_months": float(group["Time_postop"].max()),
        "n_pre_recurrence_or_no_event": len(pre),
        "n_positive_pre_recurrence_or_no_event": int(pre["positive"].sum()),
        "n_at_or_after_recurrence": int(group["at_or_after_recurrence"].sum()),
    })
point_summary = pd.DataFrame(points)
point_summary["rank"] = pd.to_numeric(point_summary["sample_point"].str.extract(r"P(\d+)")[0])
point_summary = point_summary.sort_values("rank").drop(columns="rank")
write_csv(point_summary, "samplepoint_summary.csv")

# Adjacent observed transitions. Missing scheduled visits never enter as negatives.
def transition_table(frame: pd.DataFrame, scope: str) -> list[dict]:
    rows = []
    for _, group in frame.sort_values(["Patient_No", "Time_postop", "sample_rank"]).groupby("Patient_No"):
        statuses = group["ctDNA_status"].astype(str).tolist()
        rows.extend((a, b) for a, b in zip(statuses, statuses[1:]))
    total = len(rows)
    result = []
    for a, b, label in [
        ("Negative", "Negative", "negative_to_negative"),
        ("Negative", "Positive", "negative_to_positive"),
        ("Positive", "Negative", "positive_to_negative"),
        ("Positive", "Positive", "positive_to_positive"),
    ]:
        n = sum(x == a and y == b for x, y in rows)
        result.append({"scope": scope, "transition": label, "n_pairs": n,
                       "denominator_pairs": total, "proportion": n / total if total else None})
    return result

transitions = transition_table(ctdna, "all observed consecutive draws")
transitions += transition_table(
    ctdna[(ctdna["sample_rank"] >= 2) & ctdna["pre_event"]],
    "postoperative consecutive draws before recurrence",
)
write_csv(pd.DataFrame(transitions), "MRD_transition_frequencies.csv")


def trajectory(group: pd.DataFrame) -> str:
    if group.empty:
        return "No postoperative observation"
    statuses = group.sort_values(["Time_postop", "sample_rank"])["ctDNA_status"].astype(str).tolist()
    if len(statuses) == 1:
        return "Single postoperative observation: " + statuses[0]
    first, last = statuses[0], statuses[-1]
    n_positive = sum(x == "Positive" for x in statuses)
    if n_positive == 0:
        return "Persistent negative"
    if n_positive == len(statuses):
        return "Persistent positive"
    if first == "Positive" and last == "Negative":
        return "Positive to negative at last observed test"
    if first == "Negative" and last == "Positive":
        return "Negative to positive at last observed test"
    if first == last == "Negative":
        return "Intermittent positive; negative endpoints"
    return "Intermittent negative; positive endpoints"


trajectory_rows = []
for pid in sorted(patient_ids):
    group = ctdna[(ctdna["Patient_No"] == pid) & (ctdna["sample_rank"] >= 2) & ctdna["pre_event"]]
    label = base.loc[base["Patient_No"] == pid, "ACT"].iloc[0]
    trajectory_rows.append({"Patient_No": pid, "ACT": label, "trajectory": trajectory(group)})
trajectory_summary = pd.DataFrame(trajectory_rows).groupby(["ACT", "trajectory"]).size().rename("n_patients").reset_index()
write_csv(trajectory_summary, "treatment_by_mrd_trajectory.csv")

# Last MRD strictly before recurrence and earliest postoperative positive lead-time proxy.
recurrence_rows = []
for _, patient in base[base["recurrence_event"]].sort_values("Patient_No").iterrows():
    pid, event_time = patient["Patient_No"], float(patient["R_time"])
    history = ctdna[(ctdna["Patient_No"] == pid) & (ctdna["Time_postop"] < event_time)].sort_values(["Time_postop", "sample_rank"])
    post = history[history["sample_rank"] >= 2]
    positive = post[post["positive"]]
    last = history.iloc[-1] if len(history) else None
    first_positive = positive.iloc[0] if len(positive) else None
    recurrence_rows.append({
        "Patient_No": pid, "recurrence_time_months": event_time,
        "n_samples_before_recurrence": len(history),
        "last_pre_recurrence_sample_point": last["Sample_point"] if last is not None else None,
        "last_pre_recurrence_sample_time_months": float(last["Time_postop"]) if last is not None else None,
        "months_from_last_sample_to_recurrence": event_time - float(last["Time_postop"]) if last is not None else None,
        "n_postoperative_positive_samples_before_recurrence": len(positive),
        "first_postoperative_positive_time_months": float(first_positive["Time_postop"]) if first_positive is not None else None,
        "lead_time_proxy_months": event_time - float(first_positive["Time_postop"]) if first_positive is not None else None,
    })
rec_timing = pd.DataFrame(recurrence_rows)
write_csv(rec_timing, "recurrence_sample_timing_audit.csv")

# Timepoint coverage and samples at/after recurrence by outcome group.
coverage = []
for event, label in [(False, "No documented recurrence"), (True, "Documented recurrence")]:
    ids = set(base.loc[base["recurrence_event"] == event, "Patient_No"])
    for point in [f"P{i}" for i in range(1, 10)]:
        group = ctdna[(ctdna["Patient_No"].isin(ids)) & (ctdna["Sample_point"] == point)]
        coverage.append({"outcome_group": label, "sample_point": point,
                         "patients_with_draw": int(group["Patient_No"].nunique()),
                         "patients_in_group": len(ids),
                         "coverage": group["Patient_No"].nunique() / len(ids) if ids else None,
                         "draws_at_or_after_recurrence": int(group["at_or_after_recurrence"].sum())})
write_csv(pd.DataFrame(coverage), "sample_coverage_by_recurrence.csv")

# Within-patient calendar irregularity among postoperative samples.
interval_rows = []
for pid, group in ctdna[ctdna["sample_rank"] >= 2].sort_values(["Patient_No", "Time_postop", "sample_rank"]).groupby("Patient_No"):
    times = group["Time_postop"].dropna().astype(float).to_numpy()
    gaps = np.diff(times)
    if len(gaps) >= 2 and float(np.mean(gaps)) > 0:
        interval_rows.append({"Patient_No": pid, "n_postoperative_draws": len(times),
                              "interval_cv": float(np.std(gaps, ddof=1) / np.mean(gaps)),
                              "n_intervals": len(gaps)})
interval_cv = pd.DataFrame(interval_rows)
write_csv(interval_cv, "patient_sampling_regularity.csv")

# Field-level missingness. Recurrence-time blanks are structural: event time is not
# supplied for people without documented recurrence, and cannot be used as censor time.
missing = []
def add_missing(source: str, variable: str, series: pd.Series, interpretation: str = "") -> None:
    n = len(series)
    n_miss = int(series.isna().sum())
    missing.append({"source": source, "variable": variable, "n_rows": n, "n_missing": n_miss,
                    "missing_percent": 100 * n_miss / n if n else None, "interpretation": interpretation})

for field in ["Age", "Sex", "PTL", "Stage", "Hist_type", "Hist_grade", "LV_invasion", "N_invasion", "MSI", "ACT"]:
    add_missing("Table S1 patient-level", field, base[field])
add_missing("Table S1 patient-level", "R_time", base["R_time"],
            "Structural blank for patients without documented recurrence; not a last-follow-up time.")
add_missing("Table S1 patient-level", "R_treatment", base["R_treatment"],
            "Expected blank when no recurrence is documented.")
for field in variants.columns:
    add_missing("Dataset S1 primary-tumor variants", str(field), variants[field])
for field in ["Patient_No", "Time_postop", "Sample_point", "ctDNA_status"]:
    add_missing("Dataset S2 ctDNA samples", field, ctdna[field])
for field in ["Patient_No", "Time_postop", "CEA_value"]:
    add_missing("Dataset S3 CEA samples", field, cea[field])
missing.append({
    "source": "Dataset S3 patient availability", "variable": "Preoperative CEA available",
    "n_rows": len(patient_ids), "n_missing": len(patient_ids - cea_pre_ids),
    "missing_percent": 100 * len(patient_ids - cea_pre_ids) / len(patient_ids),
    "interpretation": "Nine patients lack a preoperative CEA; one has no CEA result at any time.",
})
missing.append({
    "source": "Dataset S3 patient availability", "variable": "Any CEA available",
    "n_rows": len(patient_ids), "n_missing": len(patient_ids - cea_ids),
    "missing_percent": 100 * len(patient_ids - cea_ids) / len(patient_ids),
    "interpretation": "Eight of the nine patients without preoperative CEA have a later CEA measurement.",
})
write_csv(pd.DataFrame(missing), "missingness_matrix.csv")

# Variant and CEA summaries.
variants_per_patient = variants.groupby("Patient_No").size()
write_csv(pd.DataFrame([{
    "n_variant_rows": len(variants),
    "n_patients_with_primary_tumor_variants": variants["Patient_No"].nunique(),
    "median_variants_per_patient": float(variants_per_patient.median()),
    "q1_variants_per_patient": float(variants_per_patient.quantile(.25)),
    "q3_variants_per_patient": float(variants_per_patient.quantile(.75)),
    "min_variants_per_patient": int(variants_per_patient.min()),
    "max_variants_per_patient": int(variants_per_patient.max()),
    "source_material": "Primary tumor tissue; not serial plasma variant data",
}]), "tumor_variant_summary.csv")

cea_above5 = int(cea_pre["CEA_value"].gt(5.0).sum())
write_csv(pd.DataFrame([{
    "n_cea_measurements": len(cea),
    "n_patients_with_cea": cea["Patient_No"].nunique(),
    "n_patients_without_any_cea": len(patient_ids - cea_ids),
    "n_patients_with_preoperative_cea": len(cea_pre_ids),
    "n_patients_without_preoperative_cea": len(patient_ids - cea_pre_ids),
    "n_preoperative_measurements_time_le_0": len(cea_pre),
    "n_preoperative_cea_above_5_ng_ml": cea_above5,
    "preoperative_rate_above_5": float(cea_pre["CEA_value"].gt(5.0).mean()),
    "threshold_source": "Supplementary Methods normal range 0.00-5.00 ng/mL; elevated derived as >5.00 for this audit",
}]), "cea_summary.csv")
cea_time_groups = []
for label, group in [
    ("All CEA measurements", cea),
    ("Preoperative: Time_postop <= 0", cea_pre),
    ("Postoperative: Time_postop > 0", cea[cea["Time_postop"] > 0]),
]:
    cea_time_groups.append({
        "period": label, "n_measurements": len(group),
        "n_patients": group["Patient_No"].nunique(),
        "median_time_months": float(group["Time_postop"].median()),
        "q1_time_months": float(group["Time_postop"].quantile(.25)),
        "q3_time_months": float(group["Time_postop"].quantile(.75)),
        "median_cea_ng_ml": float(group["CEA_value"].median()),
        "q1_cea_ng_ml": float(group["CEA_value"].quantile(.25)),
        "q3_cea_ng_ml": float(group["CEA_value"].quantile(.75)),
        "n_above_5_ng_ml": int(group["CEA_value"].gt(5).sum()),
        "rate_above_5_ng_ml": float(group["CEA_value"].gt(5).mean()),
    })
write_csv(pd.DataFrame(cea_time_groups), "cea_time_summary.csv")

# Count-level replication that does not require censoring.
first_postop = ctdna[ctdna["sample_rank"] == 2].sort_values(["Patient_No", "Time_postop"]).drop_duplicates("Patient_No")
first_postop = first_postop.merge(base[["Patient_No", "recurrence_event"]], on="Patient_No", how="left")
positive_events = int(((first_postop["ctDNA_status"] == "Positive") & first_postop["recurrence_event"]).sum())
negative_events = int(((first_postop["ctDNA_status"] == "Negative") & first_postop["recurrence_event"]).sum())
write_csv(pd.crosstab(first_postop["ctDNA_status"], first_postop["recurrence_event"]).reset_index(),
          "postoperative_ctdna_by_recurrence.csv")

events = base[base["recurrence_event"]]
draw_counts = ctdna.groupby("Patient_No").size().reindex(sorted(patient_ids), fill_value=0)
post_draw_counts = ctdna[ctdna["sample_rank"] >= 2].groupby("Patient_No").size().reindex(sorted(patient_ids), fill_value=0)
lead = rec_timing["lead_time_proxy_months"].dropna()
cv_stats = quantiles(interval_cv["interval_cv"]) if len(interval_cv) else quantiles(pd.Series(dtype=float))
summary = {
    "patients": int(base["Patient_No"].nunique()),
    "recurrence_events_from_R_time": int(base["recurrence_event"].sum()),
    "stage_II_events": int(((base["Stage"] == "II") & base["recurrence_event"]).sum()),
    "stage_III_events": int(((base["Stage"] == "III") & base["recurrence_event"]).sum()),
    "ctdna_samples": len(ctdna), "ctdna_patients": ctdna["Patient_No"].nunique(),
    "ctdna_draws_per_patient": quantiles(draw_counts),
    "postoperative_draws_per_patient": quantiles(post_draw_counts),
    "samplepoint_counts": {str(k): int(v) for k, v in ctdna["Sample_point"].value_counts().items()},
    "P1_positive": int(ctdna.loc[ctdna["sample_rank"] == 1, "positive"].sum()),
    "P1_total": int((ctdna["sample_rank"] == 1).sum()),
    "P2_positive": int(ctdna.loc[ctdna["sample_rank"] == 2, "positive"].sum()),
    "P2_total": int((ctdna["sample_rank"] == 2).sum()),
    "P2_positive_recurrences": positive_events, "P2_negative_recurrences": negative_events,
    "cea_samples": len(cea), "cea_patients": cea["Patient_No"].nunique(),
    "patients_without_cea": len(patient_ids - cea_ids),
    "preop_cea_n": len(cea_pre), "preop_cea_above5": cea_above5,
    "primary_tumor_variant_rows": len(variants), "variant_patients": variants["Patient_No"].nunique(),
    "draws_at_or_after_recurrence": int(ctdna["at_or_after_recurrence"].sum()),
    "lead_time_proxy": quantiles(lead), "lead_time_proxy_n": int(len(lead)),
    "recurrence_time_months": quantiles(events["R_time"]),
    "individual_censor_times_available": 0, "death_status_or_time_available": 0,
    "sample_ids_available": 0, "ACT_yes": int(base["ACT_yes"].sum()),
    "ACT_no": int(base["ACT"].astype(str).str.lower().eq("no").sum()),
    "manifest_files": n_manifest, "verified_files": n_verified,
    "postop_interval_cv": cv_stats,
}
(QC / "audit_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

# Machine-readable dictionary, including specifically requested fields not supplied.
dictionary_rows = [
("patient_id","Patient_No","Table S1 / Datasets S1-S3","Pseudonymous participant key","Available","No unique sample ID."),
("age","Age","Table S1","Age; unit not stated in column header, paper reports years","Available",""),
("sex","Sex","Table S1","Male/Female","Available",""),
("tumor_site","PTL","Table S1","Right/left primary tumor location","Available","No exact subsite or separate colon/rectum field."),
("stage","Stage","Table S1","Pathological stage II/III","Available",""),
("T_stage","—","—","Separate T category","Unavailable","Do not infer."),
("N_stage","—","—","Separate N category","Unavailable","Do not infer."),
("pathology_risk_factors","Hist_type; Hist_grade; LV_invasion; N_invasion; MSI","Table S1","Histology, grade, lymphovascular and nerve invasion, MSI","Available","Other pathology risk factors are not supplied."),
("baseline_CEA","CEA_value at Time_postop <= 0","Dataset S3","Preoperative CEA; ng/mL per Supplementary Methods","Available","9 patients lack a preoperative result; one patient lacks any CEA result."),
("surgery_date","—","—","Absolute date of operation","Unavailable","Only relative months in biomarker tables."),
("adjuvant_chemotherapy","ACT","Table S1","Ever received ACT, yes/no","Available",""),
("chemotherapy_regimen","—","—","Regimen / cycles","Unavailable",""),
("chemotherapy_timing","—","—","Start/end and blood-draw relation","Unavailable","Cannot identify post-ACT sample patient by patient."),
("recurrence_event","R_time nonblank","Table S1","First documented recurrence","Available (derived)","Does not include separately documented CRC death."),
("recurrence_time","R_time","Table S1","Months from surgery to first documented recurrence","Available for events only","Blank for no documented recurrence; not a censor time."),
("death","—","—","Death status, date and cause","Unavailable",""),
("follow_up_duration","—","—","Individual last-follow-up / censor time","Unavailable","Paper reports cohort median 27.4 months, not individual censor times."),
("sample_id","—","—","Unique blood draw ID","Unavailable",""),
("blood_draw_time","Time_postop","Datasets S2/S3","Relative time from surgery, months","Available","P1 is pre-op but stored at 0; use Sample_point for meaning."),
("sample_point","Sample_point","Dataset S2","Protocol point P1-P9","Available",""),
("ctDNA_status","ctDNA_status","Dataset S2","Binary Positive/Negative","Available","Missing draw is not negative."),
("ctDNA_quantitative_value","—","—","Plasma ctDNA concentration","Unavailable",""),
("plasma_VAF","—","—","ctDNA VAF by sample","Unavailable","Dataset S1 DP/AD are primary-tumor reads."),
("plasma_variant_count","—","—","Number of tracked/detected variants per plasma sample","Unavailable",""),
("tumor_mutations","Patient_No; Gene; Chr_start; Chr_end; Ref; Alt; Variant_Classification; AAChange; Recurrent_in_COSMIC; DP; AD","Dataset S1","Primary-tumor mutation records/read counts","Available","No serial plasma mutation-level calls."),
("assay_QC","—","—","Per-sample QC/pass metrics","Unavailable","Methods are documented, sample QC fields are not."),
("CEA_measurement_time","Time_postop","Dataset S3","Relative time from surgery, months","Available",""),
("CEA_value","CEA_value","Dataset S3","Serum CEA, ng/mL per Methods","Available",""),
("CEA_threshold","—","Supplementary Methods","Normal range 0.00-5.00 ng/mL","Available as method definition","Audit derives elevated as >5.00."),
]

dictionary_zh = {
    'patient_id': ('伪匿名患者标识', '可用', '未提供唯一采血样本编号。'),
    'age': ('年龄；列名未注明单位，论文报告为岁', '可用', ''),
    'sex': ('男性/女性', '可用', ''),
    'tumor_site': ('原发肿瘤部位：右侧/左侧', '可用', '未提供更精确的部位或结肠/直肠单独分类。'),
    'stage': ('病理分期 II/III 期', '可用', ''),
    'T_stage': ('单独的 T 分期', '不可用', '不可推断。'),
    'N_stage': ('单独的 N 分期', '不可用', '不可推断。'),
    'pathology_risk_factors': ('组织学类型、分级、淋巴血管侵犯、神经侵犯及 MSI', '可用', '未提供其他病理风险因素。'),
    'baseline_CEA': ('术前 CEA；按补充方法，单位为 ng/mL', '可用', '9 名患者缺少术前结果；其中 1 名患者完全没有 CEA 结果。'),
    'surgery_date': ('手术绝对日期', '不可用', '生物标志物表仅提供相对月份。'),
    'adjuvant_chemotherapy': ('是否曾接受辅助化疗（是/否）', '可用', ''),
    'chemotherapy_regimen': ('化疗方案/周期数', '不可用', ''),
    'chemotherapy_timing': ('化疗起止时间及其与采血的先后关系', '不可用', '无法逐患者识别辅助化疗后的采血。'),
    'recurrence_event': ('首次记录的复发事件', '可用（派生）', '不包括单独记录的结直肠癌死亡。'),
    'recurrence_time': ('从手术到首次记录复发的月数', '仅事件患者可用', '未记录复发者为空；该值不是删失时间。'),
    'death': ('死亡状态、日期及死因', '不可用', ''),
    'follow_up_duration': ('逐患者末次随访/删失时间', '不可用', '论文报告队列随访中位数 27.4 个月，但未提供逐患者删失时间。'),
    'sample_id': ('唯一采血编号', '不可用', ''),
    'blood_draw_time': ('相对手术时间（月）', '可用', 'P1 为术前采血，但记录时间为 0；解释时须结合 Sample_point。'),
    'sample_point': ('方案规定采样点 P1-P9', '可用', ''),
    'ctDNA_status': ('二元状态：阳性/阴性', '可用', '未观测到采血不等于 MRD 阴性。'),
    'ctDNA_quantitative_value': ('血浆 ctDNA 浓度', '不可用', ''),
    'plasma_VAF': ('逐样本 ctDNA 变异等位基因频率', '不可用', 'Dataset S1 的 DP/AD 是原发肿瘤测序读数。'),
    'plasma_variant_count': ('每份血浆样本中追踪/检出的变异数', '不可用', ''),
    'tumor_mutations': ('原发肿瘤突变记录/测序读数', '可用', '未提供连续血浆样本的突变级别结果。'),
    'assay_QC': ('逐样本检测质量控制/通过指标', '不可用', '方法有所说明，但没有样本级 QC 字段。'),
    'CEA_measurement_time': ('相对手术时间（月）', '可用', ''),
    'CEA_value': ('血清 CEA；按方法部分，单位为 ng/mL', '可用', ''),
    'CEA_threshold': ('正常范围 0.00-5.00 ng/mL', '方法中有定义', '本次审计将 >5.00 ng/mL 派生为升高。'),
}
dictionary_rows = [(concept, source_field, source_file, *dictionary_zh[concept])
                   for concept, source_field, source_file, _, _, _ in dictionary_rows]

write_csv(pd.DataFrame(dictionary_rows, columns=["requested_concept","source_field","source_file","definition","availability","audit_note"]),
          "data_dictionary.csv")

# Checks intentionally separate exact reproduction from crude event-count checks.
checks = [
("Cohort size","240",str(len(patient_ids)),"MATCH" if len(patient_ids)==240 else "MISMATCH"),
("Serial ctDNA rows","1290",str(len(ctdna)),"MATCH" if len(ctdna)==1290 else "MISMATCH"),
("Recurrence events","32",str(int(base["recurrence_event"].sum())),"MATCH" if int(base["recurrence_event"].sum())==32 else "MISMATCH"),
("Pre-op ctDNA positive","154/240",f"{summary['P1_positive']}/{summary['P1_total']}","MATCH" if (summary['P1_positive'],summary['P1_total'])==(154,240) else "MISMATCH"),
("Day 3-7 ctDNA positive","20/240",f"{summary['P2_positive']}/{summary['P2_total']}","MATCH" if (summary['P2_positive'],summary['P2_total'])==(20,240) else "MISMATCH"),
("P2-positive cases with recurrence","12/20",f"{positive_events}/20","MATCH" if positive_events==12 else "MISMATCH"),
("P2-negative cases with recurrence","20/220",f"{negative_events}/220","MATCH" if negative_events==20 else "MISMATCH"),
("Pre-op CEA >5 ng/mL","91/231",f"{cea_above5}/{len(cea_pre)}","MATCH" if (cea_above5,len(cea_pre))==(91,231) else "MISMATCH"),
("Postoperative ctDNA vs RFS HR/KM","HR 10.98 (95% CI 5.31-22.72)","Not estimable","FAIL: censor/death times unavailable"),
("Post-ACT ctDNA vs RFS","HR 12.76 (95% CI 5.39-30.19)","Cannot identify first post-ACT sample","FAIL: individual ACT timing unavailable"),
("Surveillance ctDNA vs RFS","HR 32.02 (95% CI 10.79-95.08)","Cannot reconstruct analysis set","FAIL: censor and treatment timing unavailable"),
("ctDNA imaging lead time","Mean 5.01 months","Proxy only; see timing audit","FAIL: post-definitive-treatment status unavailable"),
]
write_csv(pd.DataFrame(checks, columns=["endpoint","published_result","raw_data_result","status"]), "replication_checks.csv")

# Audit report.
draw_q = quantiles(draw_counts)
rec_q = quantiles(events["R_time"])
cv_q = quantiles(interval_cv["interval_cv"])
point_rates = "; ".join(
    f"{row.sample_point}: {row.n_positive}/{row.n_samples} ({100*row.positive_rate:.1f}%)"
    for row in point_summary.itertuples()
)
transition_lines = "\n".join(
    f"| {row['scope']} | {row['transition'].replace('_', ' ')} | {row['n_pairs']} / {row['denominator_pairs']} |"
    for row in transitions
)
trajectory_lines = "\n".join(
    f"| {row.ACT} | {row.trajectory} | {row.n_patients} |"
    for row in trajectory_summary.sort_values(["ACT", "trajectory"]).itertuples()
)
draw_group_text = "; ".join(
    f"{row.outcome_group}: median {row.median:.0f} (IQR {row.q1:.1f}-{row.q3:.1f}; n={int(row.n)})"
    for row in pd.DataFrame(draw_group_rows).itertuples()
)
last_gap = quantiles(rec_timing["months_from_last_sample_to_recurrence"])
lead_mean = float(lead.mean()) if len(lead) else float("nan")
irregular_n = int((interval_cv["interval_cv"] > .2).sum())
audit_text = f"""# Data audit report

## Scope and provenance

This is Phase 0-1 source retrieval and audit for the prospective stage II-III CRC cohort reported by [Chen et al.](https://doi.org/10.1186/s13045-021-01089-z). The full set of 18 supplementary records linked to the article DOI was downloaded from Springer Nature Figshare. All {n_verified}/{n_manifest} downloaded files passed the publisher-provided MD5 check. File IDs, DOIs, URLs, byte counts, and checksums are stored in data_raw/figshare_manifest.csv.

The patient-level sources are [Additional file 1](https://springernature.figshare.com/articles/journal_contribution/Additional_file_1_of_Postoperative_circulating_tumor_DNA_as_markers_of_recurrence_risk_in_stages_II_to_III_colorectal_cancer/14608859) (Table S1, DOCX), [Additional file 6](https://springernature.figshare.com/articles/dataset/Additional_file_6_of_Postoperative_circulating_tumor_DNA_as_markers_of_recurrence_risk_in_stages_II_to_III_colorectal_cancer/14608874) (Dataset S1, primary-tumor variants), [Additional file 7](https://springernature.figshare.com/articles/dataset/Additional_file_7_of_Postoperative_circulating_tumor_DNA_as_markers_of_recurrence_risk_in_stages_II_to_III_colorectal_cancer/14608877) (Dataset S2, serial ctDNA), and [Additional file 9](https://springernature.figshare.com/articles/dataset/Additional_file_9_of_Postoperative_circulating_tumor_DNA_as_markers_of_recurrence_risk_in_stages_II_to_III_colorectal_cancer/14608883) (Dataset S3, serial CEA). Other supplementary items contain methods, the gene panel, figures, and aggregate results. No absent patient-level data were imputed.

## Cohort and sample reconciliation

- Unique patients in Table S1: **{len(patient_ids)}**.
- Nonblank recurrence times: **{int(base['recurrence_event'].sum())}** ({summary['stage_II_events']} stage II; {summary['stage_III_events']} stage III).
- Serial ctDNA observations: **{len(ctdna)}** for {ctdna['Patient_No'].nunique()} patients.
- Draws per patient including P1: median **{draw_q['median']:.0f}**, IQR **{draw_q['q1']:.1f}-{draw_q['q3']:.1f}**, range **{draw_q['min']:.0f}-{draw_q['max']:.0f}**.
- CEA observations: **{len(cea)}** across {cea['Patient_No'].nunique()} patients. There are {len(cea_pre)} preoperative results across {len(cea_pre_ids)} patients, so {len(patient_ids-cea_pre_ids)} patients lack preoperative CEA. One patient has no CEA result at any time.
- Primary-tumor variant records: **{len(variants)}** across {variants['Patient_No'].nunique()} patients. Per-patient count median **{variants_per_patient.median():.1f}** (IQR {variants_per_patient.quantile(.25):.1f}-{variants_per_patient.quantile(.75):.1f}; range {variants_per_patient.min()}-{variants_per_patient.max()}). These are tumor-tissue records, not plasma calls.

See outputs/tables/patient_draw_counts.csv, samplepoint_summary.csv, tumor_variant_summary.csv, and replication_checks.csv.

## Data dictionary and field availability

The field-level machine-readable dictionary is outputs/tables/data_dictionary.csv. Table S1 contains Patient_No, Age, Sex, PTL, Stage, Hist_type, Hist_grade, LV_invasion, N_invasion, MSI, ACT, R_time, and R_treatment. R_time is months from surgery to first documented recurrence. A blank does not encode last follow-up.

Dataset S2 contains only Patient_No, Time_postop (months), Sample_point (P1-P9), and binary ctDNA_status. There is no unique sample ID, quantitative plasma ctDNA value, plasma VAF, plasma mutation count, or per-sample assay QC. Dataset S1 contains primary-tumor mutation-level data with DP and AD; those read counts must not be interpreted as plasma ctDNA values. Dataset S3 contains Patient_No, Time_postop (months), and CEA_value. Supplementary Methods define the normal CEA range as 0.00-5.00 ng/mL; this audit derives elevated as >5.00 ng/mL.

The separate T stage and N stage, absolute surgery date, ACT regimen and timing, death status/date/cause, and per-patient censor/last-follow-up time are unavailable. ACT is only yes/no. Sample times are in decimal months rather than dates or days. P1 is labeled preoperative but stored at time 0, so interpret it using Sample_point as well as numeric time.

## Longitudinal structure and MRD

Counts, positivity, and actual time distributions by observed sample point are in outputs/tables/samplepoint_summary.csv. Positive rates over all records are: {point_rates}. The table also reports counts strictly before known recurrence and rows at/after recurrence. At/after-recurrence ctDNA observations are post-outcome and cannot enter a future-recurrence model.

P1 is preoperative and P2 is the early postoperative draw. Retain observed decimal-month time rather than assigning exact dates to protocol points. Transition counts are in outputs/tables/MRD_transition_frequencies.csv. Their denominators include adjacent observed draws only.

| Observed-history scope | Transition | Pairs |
|---|---|---:|
{transition_lines}

An unobserved draw is never treated as MRD-negative. Treatment-by-trajectory counts use the ACT yes/no field and, among event cases, include only samples strictly before recurrence:

| ACT | Observed postoperative trajectory | Patients |
|---|---|---:|
{trajectory_lines}

Postoperative collection intervals are irregular. Among patients with at least three postoperative samples, the within-patient coefficient of variation of consecutive intervals has median {cv_q['median']:.2f} (IQR {cv_q['q1']:.2f}-{cv_q['q3']:.2f}); {irregular_n}/{len(interval_cv)} exceed the descriptive CV 0.20 marker. This is not a clinical threshold. Timepoint coverage and draw-count summaries by recurrence group are in outputs/tables/sample_coverage_by_recurrence.csv and draw_counts_by_recurrence.csv. {draw_group_text}.

## Recurrence, follow-up, and informative missingness

For event cases, R_time has median {rec_q['median']:.2f} months (IQR {rec_q['q1']:.2f}-{rec_q['q3']:.2f}). This is event time among recurrences, not cohort follow-up. No patient-level censor time is supplied for patients without documented recurrence, and no death field is present. The article reports median cohort follow-up of 27.4 months (95% CI 26.2-28.5); the supplement does not permit independent reconstruction of recurrence-versus-nonrecurrence follow-up duration.

For recurrence cases, outputs/tables/recurrence_sample_timing_audit.csv lists the last MRD observation strictly before recurrence and an exploratory interval from earliest postoperative positive MRD to recurrence. The last pre-recurrence draw is a median {last_gap['median']:.2f} months before the event (IQR {last_gap['q1']:.2f}-{last_gap['q3']:.2f}). A positive-to-recurrence proxy is available for {len(lead)} event patients (mean {lead_mean:.2f}, median {lead.median():.2f} months). It is not the paper's lead-time definition because ACT dates are absent and the paper uses ctDNA detection after definitive treatment. Samples at or after recurrence are excluded from this proxy; {int(ctdna['at_or_after_recurrence'].sum())} raw observations are at/after recurrence and must not enter prediction histories.

Observation intensity can be informative. Later visit coverage differs across recurrence groups, the planned schedule is not fully observed, and censoring/follow-up duration for non-events is unavailable. These files cannot distinguish missed visits, early loss to follow-up, treatment, death, withdrawal, or administrative censoring. Missing visits were not coded as negative; no observation-process model was attempted.

## Phase 0-1 conclusion

The public data support a patient-level clinical table, binary serial ctDNA with relative sample times, primary-tumor variants, and numeric serial CEA. They do not support all requested variables. Quantitative plasma MRD, per-sample QC, detailed ACT timing, individual censor times, and death data are unavailable.
"""
(REPORTS / "data_audit.md").write_text(audit_text, encoding="utf-8")
(REPORTS / "data_audit.md").write_text((ROOT / "scripts" / "templates" / "data_audit_zh.md").read_text(encoding="utf-8"), encoding="utf-8")

print(json.dumps({
    "patients": len(patient_ids),
    "recurrences": int(base["recurrence_event"].sum()),
    "ctdna_rows": len(ctdna),
    "draws_per_patient": quantiles(draw_counts),
    "P1_positive": f"{summary['P1_positive']}/{summary['P1_total']}",
    "P2_positive": f"{summary['P2_positive']}/{summary['P2_total']}",
    "P2_positive_recurrence": f"{positive_events}/20",
    "P2_negative_recurrence": f"{negative_events}/220",
    "CEA_preop_above5": f"{cea_above5}/{len(cea_pre)}",
    "lead_time_proxy": quantiles(lead),
}, indent=2))
