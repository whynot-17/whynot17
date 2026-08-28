# Step 10C — CCHS public association demonstration

## Frozen analysis

The feasibility audit was completed before this model was run. The sole frozen demonstration is:

> **Current smoking status (`SMKDVSTY`) → self-reported high blood pressure (`CCC_80`) among adults**

The model uses categorical adult age group, sex, adult BMI classification, education, and household income. Point estimation uses `WTS_M`; precision uses the 1,000 supplied CCHS bootstrap replicate weights (`BSW1–BSW1000`). The bootstrap standard error is the empirical standard deviation with `ddof=1` around the replicate mean. No CCHS-specific Fay multiplier is applied.

## Result

- Analytic N: **55,186**
- Unweighted high-blood-pressure cases: **14,955**
- Unweighted current smokers: **5,905**
- Current smoker versus not-currently-smoking OR: **1.034**
- 95% CI: **0.909–1.175**
- Bootstrap percentile 95% CI sensitivity: **0.894–1.158**
- Normal-reference P value: **0.611656**
- Bootstrap fits converged: **1000/1000**

This is a source-native population-replacement demonstration, not a replication of the 29 NHANES urinary biomarker tests. It is also not a causal estimate: CCHS is cross-sectional and the exposure/outcome are contemporaneous survey measures.

## Firewall

The exposure–outcome pair was frozen in `CCHS_PUBLIC_AUDIT_QC_SUMMARY.json` before this script was run. No association estimate, P value, confidence interval, or FDR was used to choose the pair. No alternative disease was searched after seeing this result.

The raw CCHS archive and bootstrap matrix remain outside Git under `D:\whynot17\public_sources\cchs_2022`; only aggregate result/provenance files are versioned.
