"""Phase 2 count-level reproduction and explicit RFS replication feasibility report."""
from __future__ import annotations
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
QC = ROOT / "outputs" / "qc"
REPORTS = ROOT / "reports"
REPORTS.mkdir(parents=True, exist_ok=True)

summary = json.loads((QC / "audit_summary.json").read_text(encoding="utf-8"))
checks = pd.read_csv(TABLES / "replication_checks.csv", encoding="utf-8-sig")
timing = pd.read_csv(TABLES / "recurrence_sample_timing_audit.csv", encoding="utf-8-sig")
lead = timing["lead_time_proxy_months"].dropna()
lead_mean = float(lead.mean()) if len(lead) else float("nan")
lead_median = float(lead.median()) if len(lead) else float("nan")

check_rows = "\n".join(
    f"| {row.endpoint} | {row.raw_data_result} | {row.published_result} | {row.status} |"
    for row in checks.itertuples()
)

report = f"""# Original-paper replication report

## Status

**Phase 2 count-level checks match the paper, but the specified RFS replications cannot be reconstructed from the public patient-level files. The stop rule is triggered, so no Phase 3+ modeling was started.**

## Checks

| Endpoint | Result from released files | Published result | Status |
|---|---|---|---|
{check_rows}

Raw file counts confirm {summary['patients']} patients, {summary['ctdna_samples']} serial ctDNA rows, and {summary['recurrence_events_from_R_time']} recurrence times ({summary['stage_II_events']} stage II; {summary['stage_III_events']} stage III). At P2, {summary['P2_positive']} of {summary['P2_total']} patients are ctDNA-positive. Documented recurrence occurred in {summary['P2_positive_recurrences']}/20 P2-positive patients and {summary['P2_negative_recurrences']}/220 P2-negative patients. These last two are crude event-status counts, not survival estimates.

## Why RFS cannot be reproduced

1. **Postoperative MRD versus RFS:** The paper defines RFS from surgery to radiological recurrence or CRC death, censored at last follow-up or non-CRC death. Table S1 supplies recurrence time for event cases, but no death status/date and no individual censor time for the other {summary['patients']-summary['recurrence_events_from_R_time']} patients. Kaplan-Meier curves, 2-year RFS, log-rank tests, and HR 10.98 (95% CI 5.31-22.72) cannot be estimated without inventing censoring times.
2. **First post-ACT MRD versus RFS:** ACT is only recorded yes/no. Regimen and start/end dates are absent, so the first post-ACT blood draw cannot be selected patient by patient. The published HR 12.76 (95% CI 5.39-30.19) is not reproducible.
3. **Surveillance MRD versus RFS:** The paper's surveillance subset requires post-definitive-treatment sampling and at least 24 months of follow-up or relapse. Individual ACT completion dates and non-event censor times are missing, so neither eligibility nor HR 32.02 (95% CI 10.79-95.08) can be reconstructed.
4. **ctDNA lead time:** There are {len(lead)} recurrence patients with a positive postoperative sample strictly before recorded recurrence. The earliest postoperative-positive to recurrence proxy has mean {lead_mean:.2f} and median {lead_median:.2f} months. This is not equivalent to the published 5.01-month lead time because individual ACT timing is unavailable and the publication uses post-definitive-treatment ctDNA detection.

## Decision

The raw counts support partial verification of cohort composition and crude recurrence status. They do not support the paper's RFS analyses. The pipeline stops at Phase 2 as requested. No fitted survival model, digital twin, or patient-specific future risk estimate is available or claimed.

To resume the planned modeling phases, obtain an authorized patient-level file with last-follow-up/censor time, death status/cause/time, ACT regimen and dates, and post-definitive-treatment sample status. Otherwise, revise the protocol before starting later phases.

## Source

Chen et al., *Postoperative circulating tumor DNA as markers of recurrence risk in stages II to III colorectal cancer*, Journal of Hematology & Oncology (2021), [DOI 10.1186/s13045-021-01089-z](https://doi.org/10.1186/s13045-021-01089-z). The paper reports 240 evaluable patients, 1,290 serial plasma samples, 32 radiological recurrences, and a mean ctDNA lead time of 5.01 months. The complete supplementary package and MD5 manifest are archived in data_raw/.
"""
(REPORTS / "replication_report.md").write_text(report, encoding="utf-8")
(REPORTS / "replication_report.md").write_text((ROOT / "scripts" / "templates" / "replication_report_zh.md").read_text(encoding="utf-8"), encoding="utf-8")
print(f"Wrote {REPORTS / 'replication_report.md'} from {len(checks)} raw-data checks.")
