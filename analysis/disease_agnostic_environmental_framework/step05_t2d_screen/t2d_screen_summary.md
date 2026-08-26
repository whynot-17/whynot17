# T2D disease plug-in screen

Generated (UTC): 2026-08-26T11:28:57.650382+00:00

## Frozen scope

- Frozen Step 4 tests entered: **29**.
- Models with finite P values: **29/29**.
- Nominal P<0.05: **14**.
- Suggestive BH-FDR q<0.10: **14**.
- Primary BH-FDR q<0.05: **14**.
- BH-FDR denominator is fixed at 29; no test was removed or re-ranked before correction.

## Outcome definition

Adults aged >=20 years were classified using the Diabetes Questionnaire and official glycohemoglobin files. Diagnosed diabetes was defined by DIQ010=1, with a conservative exclusion for likely early-onset insulin-dependent cases (current insulin use and reported diagnosis age <20). Probable undiagnosed diabetes was restricted to DIQ010=2 with available LBXGH >=6.5%. Controls were restricted to DIQ010=2 with available LBXGH <6.5%; borderline, unknown, missing, and other ambiguous responses remained indeterminate.

## Primary model

`T2D ~ log2(exposure) + age + sex + race/ethnicity + BMI + smoking + PIR`; urinary biomarkers additionally include `log2(urinary creatinine)`. Each test uses its own assay-specific laboratory file, cycle coverage and subsample weight. Outcome and covariates are constructed without using any exposure laboratory file.

## Primary screen result

| Biomarker | N | T2D cases | OR/doubling | 95% CI | P | q (29) | Status |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- |
| URXUUR | 11854 | 1715 | 1.16486 | 1.09944–1.23418 | 7.12723e-07 | 2.0669e-05 | ok |
| URXUMO | 14249 | 2103 | 1.18955 | 1.10663–1.27868 | 4.72494e-06 | 6.85117e-05 | ok |
| URXUTU | 14265 | 2110 | 1.1353 | 1.06901–1.20571 | 5.14598e-05 | 0.000497445 | ok |
| URXUPB | 14343 | 2120 | 0.849645 | 0.777042–0.929032 | 0.000423815 | 0.00307266 | ok |
| URXMIB | 12994 | 1940 | 1.09569 | 1.03588–1.15894 | 0.00160154 | 0.00928892 | ok |
| URXUBA | 14225 | 2095 | 0.922872 | 0.875578–0.97272 | 0.00301592 | 0.014577 | ok |
| URXP02 | 11519 | 1618 | 1.08497 | 1.0253–1.14811 | 0.00506748 | 0.0209938 | ok |
| URXUSN | 5842 | 977 | 1.12491 | 1.03466–1.22304 | 0.0065565 | 0.0237673 | ok |
| URXMHH | 12994 | 1940 | 1.05442 | 1.01278–1.09777 | 0.0103159 | 0.03324 | ok |
| URXCOP | 10236 | 1617 | 1.05054 | 1.0094–1.09337 | 0.0160421 | 0.0458677 | ok |
| URXMOH | 12994 | 1940 | 1.05218 | 1.00912–1.09707 | 0.0173981 | 0.0458677 | ok |
| URXUSR | 4485 | 712 | 0.886522 | 0.80212–0.979805 | 0.0193472 | 0.0465276 | ok |
| URXECP | 11581 | 1786 | 1.05444 | 1.00779–1.10326 | 0.0220479 | 0.0465276 | ok |
| LBXPFHS | 11500 | 1744 | 0.933812 | 0.880627–0.990209 | 0.0224616 | 0.0465276 | ok |
| URXP04 | 12908 | 1893 | 1.06055 | 0.993395–1.13225 | 0.0777726 | 0.146481 | ok |
| URXMBP | 14202 | 2073 | 1.06103 | 0.992679–1.13409 | 0.0808169 | 0.146481 | ok |
| LBXPFDE | 11500 | 1744 | 0.953303 | 0.893694–1.01689 | 0.145202 | 0.247697 | ok |
| URXUSB | 14255 | 2105 | 1.04516 | 0.981983–1.11239 | 0.16369 | 0.263723 | ok |
| URXP10 | 12881 | 1890 | 1.03796 | 0.973216–1.10701 | 0.254703 | 0.388757 | ok |
| URXMEP | 14198 | 2072 | 0.986768 | 0.953426–1.02128 | 0.445115 | 0.608222 | ok |
| URXUCO | 14343 | 2120 | 1.02893 | 0.955772–1.10768 | 0.446132 | 0.608222 | ok |
| URXUCS | 14343 | 2120 | 1.04608 | 0.927297–1.18007 | 0.46141 | 0.608222 | ok |
| LBXPFNA | 11500 | 1744 | 0.977527 | 0.911542–1.04829 | 0.520956 | 0.656857 | ok |
| URXMHP | 14202 | 2073 | 1.00963 | 0.966339–1.05487 | 0.666244 | 0.805044 | ok |
| URXUCD | 11882 | 1828 | 1.01323 | 0.938067–1.09442 | 0.736242 | 0.825558 | ok |
| URXMZP | 14202 | 2073 | 0.991881 | 0.944946–1.04115 | 0.740156 | 0.825558 | ok |
| URXP25 | 4406 | 752 | 0.984009 | 0.87288–1.10929 | 0.787675 | 0.846022 | ok |
| URXUTL | 14309 | 2109 | 0.991281 | 0.883108–1.1127 | 0.881179 | 0.892114 | ok |
| URXBPH | 7204 | 996 | 0.994929 | 0.92362–1.07174 | 0.892114 | 0.892114 | ok |

## Interpretation firewall

This is the first outcome-aware T2D demonstration of the frozen environmental test set. No GeneCards, disease-specific CTD, transcriptomic, literature, robustness, or mechanistic analysis was performed in this stage. The output is a screening result and does not establish causality.

Nominal signals: LBXPFHS, URXCOP, URXECP, URXMHH, URXMIB, URXMOH, URXP02, URXUBA, URXUMO, URXUPB, URXUSN, URXUSR, URXUTU, URXUUR.
Suggestive q<0.10 signals: LBXPFHS, URXCOP, URXECP, URXMHH, URXMIB, URXMOH, URXP02, URXUBA, URXUMO, URXUPB, URXUSN, URXUSR, URXUTU, URXUUR.
q<0.05 signals: LBXPFHS, URXCOP, URXECP, URXMHH, URXMIB, URXMOH, URXP02, URXUBA, URXUMO, URXUPB, URXUSN, URXUSR, URXUTU, URXUUR.
