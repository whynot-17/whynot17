# URXUBA focused check

This is a read-only focused check of the frozen URXUBA primary and completed robustness CSV files. No association model was refit; no quartiles, splines, figures, external database, literature, or mechanism analysis were run.

## Direct answer

- Obesity: male beta = 0.101688 (SE 0.026859; z 3.786); female beta = 0.160704 (SE 0.020885; z 7.695). Interaction beta = 0.049764, P = 0.1050002, fixed-406 q = 0.4855722; N = 15057, cases = 5436.
- Myocardial infarction: male beta = -0.072102 (SE 0.047494; z -1.518); female beta = 0.005638 (SE 0.058474; z 0.096). Interaction beta = 0.069561, P = 0.3052515, fixed-406 q = 0.7163706; N = 15270, cases = 659.
- The strongest positive male association in the seven-outcome context set is Obesity (beta 0.101688); the strongest positive female association is also Obesity (beta 0.160704). The largest context interaction is CHF (beta 0.277058, fixed-406 q 0.0147906), not MI.
- The formal interaction coefficients for both obesity and MI are positive (female minus male), so both target edges are female-enhanced. Because the top male and female destinations are the same outcome and MI is not the strongest female association, the reciprocal “male obesity vs female MI” story is not supported.

## Existing robustness diagnostics

- Obesity: creatinine-adjusted beta = 0.049764; 10 LOCO refits; LOCO beta range 0.031739 to 0.070414; sign reversals = 0; winsorized beta = 0.049872; upper-1%-deleted beta = 0.053524; above-LOD beta = 0.041938; cycle heterogeneity P = 0.1612350.
- Myocardial infarction: creatinine-adjusted beta = 0.069561; 10 LOCO refits; LOCO beta range 0.019088 to 0.106969; sign reversals = 0; winsorized beta = 0.070724; upper-1%-deleted beta = 0.054011; above-LOD beta = 0.082657; cycle heterogeneity P = 0.0090972.
- Both target interactions retain a positive direction in every existing diagnostic, with no LOCO sign reversal. This is directional robustness, not evidence that every sensitivity estimate is statistically significant.

## Final classification

**B. Female-dominant multi-outcome susceptibility.** URXUBA shows female-enhanced coefficients for both obesity and MI, but not a clean reciprocal disease split; its strongest context signal is female-enhanced CHF rather than female MI.

Recommendation: continue URXUBA as a female-enhanced multi-outcome candidate, while avoiding the unsupported “male obesity versus female MI” description.
