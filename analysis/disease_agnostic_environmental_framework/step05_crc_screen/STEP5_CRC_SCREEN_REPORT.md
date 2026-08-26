# Step 5 CRC screen report

## Frozen scope

- Frozen unique NHANES tests entered: **29**.
- Models with finite P values: **28/29**.
- Nominal P<0.05 tests: **5**.
- BH-FDR<0.05 tests with fixed denominator 29: **2**.
- Model warning/non-ok statuses: **6**.

Primary model: `CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR`; urinary biomarkers additionally include `log2(creatinine)`. Each analyte uses its own laboratory/subsample weight and cycle-pooled strata/PSU.

## Key prespecified biomarkers

- URXCOP: OR=1.24551; 95% CI 1.07753–1.43968; P=0.00331136; BH-FDR=0.0480147; status=ok
- LBXPFHS: OR=0.624414; 95% CI 0.470702–0.828322; P=0.00146008; BH-FDR=0.0423423; status=converged_with_warning

## Nominal and FDR-supported signals

Nominal P<0.05: LBXPFHS, URXCOP, URXMHH, URXMOH, URXUCO.
BH-FDR<0.05: LBXPFHS, URXCOP.

No test was removed or re-ranked before the 29-test BH-FDR calculation. This screen is an association analysis of prevalent CRC and does not establish temporality or causality.
