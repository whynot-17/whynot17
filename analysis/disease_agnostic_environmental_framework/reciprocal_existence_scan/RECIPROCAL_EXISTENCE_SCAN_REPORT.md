# Reciprocal effect-existence screen

This screen starts from the 80 interaction-supported directional pairs in the frozen primary audit. It reads the frozen primary coefficients and SEs only, adds the pre-specified endpoint existence rule `z = beta / SE >= 1.96`, and applies no new FDR, score, robustness analysis, model refit, figure, quartile, spline, literature, or mechanism analysis.

## Rule

- Male endpoint A: `male_beta > 0`, `female_beta <= 0`, and `male_z >= 1.96`.
- Female endpoint B: `female_beta > 0`, `male_beta <= 0`, and `female_z >= 1.96`.
- Both pooled interaction coefficients must retain the expected directions: A interaction < 0 and B interaction > 0.
- Existing fixed-406 BH q-values are reported descriptively; both endpoint q < 0.05 is not required and no new multiplicity adjustment is calculated.

## Results

- Input: 80 interaction-supported pairs from the prior frozen-primary directional scan.
- Pairs passing both endpoint z thresholds: 2 across 2 exposures.
- No robustness results were read, and no model was refit.

- URXP02: male thyroid_disease (βM 0.122402, SE 0.043606, z 2.807; interaction β -0.104123, P 0.0228304, fixed-406 q 0.2439246) -> female hypertension (βF 0.104517, SE 0.026705, z 3.914; interaction β 0.072253, P 0.0195787, fixed-406 q 0.2148370).
- URXUSN: male any_cancer_history (βM 0.130305, SE 0.060664, z 2.148; interaction β -0.161470, P 0.0299013, fixed-406 q 0.2759071) -> female asthma (βF 0.090575, SE 0.039344, z 2.302; interaction β 0.117764, P 0.0341868, fixed-406 q 0.2891636).

These pairs satisfy the requested direction and effect-existence screen, but neither endpoint is automatically FDR-confirmed by this screen. They are the only current candidates for a separately pre-specified focused robustness review.
