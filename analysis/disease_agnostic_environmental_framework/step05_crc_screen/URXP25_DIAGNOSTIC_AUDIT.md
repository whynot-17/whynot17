# URXP25 diagnostic audit

Generated (UTC): 2026-08-26T10:37:08.275741+00:00

## Finding

URXP25 did not fail because it had too few CRC cases. The Step 5 outcome/covariate frame was rebuilt from the MBzP/phthalate laboratory file (`PHTHTE_*`), while URXP25 is measured in PAH files (`PAH_H/I/J`). The participant IDs in the three URXP25 source files had zero intersection with that phthalate-shaped frame.

- Registry/exposure source rows: **9,013**; non-missing positive exposure values: **8,459**.
- Pre-repair Step 5 merge rows: **0**; analytic N: **0**; status: **not_estimable**.
- Corrected cycle-matched outcome/covariate merge rows: **9,013**; complete-case N: **4,275**; CRC cases: **30**.

## Isolated corrected reanalysis

Using the same Step 5 estimator and covariate model on the corrected PAH-compatible frame: **OR=1.14991**, 95% CI **0.626286–2.11134**, P **0.645585**, status **ok**.

This isolated result is not inserted into `full_29_test_crc_screen.csv`, is not assigned a BH-FDR, and does not alter any of the other 28 tests.

## Important scope note

The zero-intersection diagnosis exposes a broader Step 5 architecture issue: the existing harmonized frame is phthalate-subsample shaped, so non-phthalate tests may not have been evaluated in their own laboratory subsamples. This audit intentionally does not rerun those tests; before the 29-test screen is treated as final, each non-phthalate test should be checked with its own cycle-matched outcome/covariate frame.
