# MCOP–CRC Phase 3A：WHI power planning simulation

本文件是 WHI biospecimen access 确认前的规划模拟，不是 WHI 实际回归结果。模拟使用 NHANES Phase 2 女性 complete-case 的 MCOP/covariate 分布作为代理，模拟一个病例对应 1 或 2 个 matched controls，并用 conditional logistic likelihood 生成和拟合病例状态。

- Simulation seed: `20260822`
- Replicates per cell: `1000`
- Alpha: 0.05, two-sided Wald test
- NHANES proxy rows: `5065`
- Proxy log2(MCOP) SD: `2.07267`
- Target effect is OR per MCOP doubling

## Interpretation

Power is scenario planning only. It will change with the actual WHI MCOP distribution, residual urine availability, matching factors, missingness, assay batch structure and the number of analyzable CRC cases. The table should be used to decide whether the WHI sample can support directionally informative replication, not to claim that WHI is already analyzable.

## Results

| CRC cases | Controls/case | Total N | Target OR | Convergence | Power | Median estimated OR | 2.5–97.5% estimated OR |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 181 | 1 | 362 | 1.15 | 1.000 | 0.622 | 1.159 | 1.024–1.331 |
| 181 | 1 | 362 | 1.20 | 1.000 | 0.854 | 1.220 | 1.066–1.419 |
| 181 | 1 | 362 | 1.25 | 1.000 | 0.932 | 1.263 | 1.106–1.471 |
| 181 | 1 | 362 | 1.30 | 1.000 | 0.988 | 1.323 | 1.161–1.572 |
| 181 | 2 | 543 | 1.15 | 1.000 | 0.790 | 1.155 | 1.048–1.297 |
| 181 | 2 | 543 | 1.20 | 1.000 | 0.938 | 1.208 | 1.085–1.356 |
| 181 | 2 | 543 | 1.25 | 1.000 | 0.993 | 1.260 | 1.134–1.419 |
| 181 | 2 | 543 | 1.30 | 1.000 | 1.000 | 1.310 | 1.178–1.490 |
| 150 | 1 | 300 | 1.15 | 1.000 | 0.548 | 1.165 | 1.009–1.372 |
| 150 | 1 | 300 | 1.20 | 1.000 | 0.771 | 1.224 | 1.054–1.449 |
| 150 | 1 | 300 | 1.25 | 1.000 | 0.895 | 1.275 | 1.100–1.528 |
| 150 | 1 | 300 | 1.30 | 1.000 | 0.962 | 1.335 | 1.143–1.604 |
| 150 | 2 | 450 | 1.15 | 1.000 | 0.712 | 1.154 | 1.032–1.297 |
| 150 | 2 | 450 | 1.20 | 1.000 | 0.884 | 1.208 | 1.072–1.384 |
| 150 | 2 | 450 | 1.25 | 1.000 | 0.977 | 1.259 | 1.129–1.445 |
| 150 | 2 | 450 | 1.30 | 1.000 | 0.994 | 1.317 | 1.171–1.511 |
| 120 | 1 | 240 | 1.15 | 1.000 | 0.456 | 1.169 | 0.993–1.423 |
| 120 | 1 | 240 | 1.20 | 1.000 | 0.664 | 1.215 | 1.045–1.496 |
| 120 | 1 | 240 | 1.25 | 1.000 | 0.829 | 1.289 | 1.081–1.600 |
| 120 | 1 | 240 | 1.30 | 1.000 | 0.904 | 1.347 | 1.123–1.708 |
| 120 | 2 | 360 | 1.15 | 1.000 | 0.605 | 1.157 | 1.017–1.336 |
| 120 | 2 | 360 | 1.20 | 1.000 | 0.812 | 1.211 | 1.064–1.400 |
| 120 | 2 | 360 | 1.25 | 1.000 | 0.948 | 1.268 | 1.111–1.465 |
| 120 | 2 | 360 | 1.30 | 1.000 | 0.978 | 1.317 | 1.147–1.534 |
| 100 | 1 | 200 | 1.15 | 1.000 | 0.405 | 1.172 | 0.976–1.462 |
| 100 | 1 | 200 | 1.20 | 1.000 | 0.577 | 1.231 | 1.008–1.557 |
| 100 | 1 | 200 | 1.25 | 1.000 | 0.705 | 1.286 | 1.045–1.653 |
| 100 | 1 | 200 | 1.30 | 1.000 | 0.859 | 1.348 | 1.104–1.758 |
| 100 | 2 | 300 | 1.15 | 1.000 | 0.542 | 1.162 | 1.011–1.368 |
| 100 | 2 | 300 | 1.20 | 1.000 | 0.737 | 1.212 | 1.047–1.404 |
| 100 | 2 | 300 | 1.25 | 1.000 | 0.897 | 1.273 | 1.094–1.506 |
| 100 | 2 | 300 | 1.30 | 1.000 | 0.962 | 1.326 | 1.136–1.566 |
| 80 | 1 | 160 | 1.15 | 1.000 | 0.313 | 1.180 | 0.937–1.516 |
| 80 | 1 | 160 | 1.20 | 1.000 | 0.510 | 1.249 | 1.008–1.638 |
| 80 | 1 | 160 | 1.25 | 1.000 | 0.649 | 1.312 | 1.044–1.775 |
| 80 | 1 | 160 | 1.30 | 1.000 | 0.737 | 1.351 | 1.075–1.869 |
| 80 | 2 | 240 | 1.15 | 1.000 | 0.446 | 1.164 | 0.974–1.399 |
| 80 | 2 | 240 | 1.20 | 1.000 | 0.639 | 1.221 | 1.026–1.461 |
| 80 | 2 | 240 | 1.25 | 1.000 | 0.805 | 1.273 | 1.074–1.538 |
| 80 | 2 | 240 | 1.30 | 1.000 | 0.910 | 1.337 | 1.123–1.650 |

## Limitations frozen before WHI access

1. NHANES is an exposure-distribution proxy, not a substitute for WHI biomarker data.
2. The simulation sets nuisance outcome coefficients to zero; covariates are retained in the fitted model and their empirical correlation with MCOP is preserved through the NHANES proxy rows.
3. Matching is represented by conditional-logistic matched sets, but exact WHI matching-factor distributions are not yet known.
4. This is not a sample-size justification until WHI confirms analyzable urine samples and actual MCOP assay performance.

## Next gate

Confirm the number of CRC cases with sufficient prediagnostic urine volume, matched controls, sample timing and assay availability in WHI Query Builder/CCC. Then rerun this simulation with the observed WHI exposure variance and missingness assumptions before finalizing the sample request.