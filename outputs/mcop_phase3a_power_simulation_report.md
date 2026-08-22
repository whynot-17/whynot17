# MCOP–CRC Phase 3A：WHI power planning simulation

本文件是 WHI biospecimen access 确认前的规划模拟，不是 WHI 实际回归结果。模拟使用 NHANES Phase 2 complete-case 的 MCOP/covariate 分布作为代理，模拟一个病例对应 1 或 2 个 matched controls，并用 conditional logistic likelihood 生成和拟合病例状态。

- Simulation seed: `20260822`
- Replicates per cell: `1000`
- Alpha: 0.05, two-sided Wald test
- NHANES proxy rows: `9936`
- Proxy log2(MCOP) SD: `2.05863`
- Target effect is OR per MCOP doubling

## Interpretation

Power is scenario planning only. It will change with the actual WHI MCOP distribution, residual urine availability, matching factors, missingness, assay batch structure and the number of analyzable CRC cases. The table should be used to decide whether the WHI sample can support directionally informative replication, not to claim that WHI is already analyzable.

## Results

| CRC cases | Controls/case | Total N | Target OR | Convergence | Power | Median estimated OR | 2.5–97.5% estimated OR |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 181 | 1 | 362 | 1.15 | 1.000 | 0.628 | 1.162 | 1.025–1.340 |
| 181 | 1 | 362 | 1.20 | 1.000 | 0.841 | 1.215 | 1.070–1.410 |
| 181 | 1 | 362 | 1.25 | 1.000 | 0.942 | 1.273 | 1.110–1.492 |
| 181 | 1 | 362 | 1.30 | 1.000 | 0.980 | 1.327 | 1.146–1.557 |
| 181 | 2 | 543 | 1.15 | 1.000 | 0.785 | 1.154 | 1.043–1.287 |
| 181 | 2 | 543 | 1.20 | 1.000 | 0.934 | 1.207 | 1.085–1.363 |
| 181 | 2 | 543 | 1.25 | 1.000 | 0.994 | 1.261 | 1.135–1.405 |
| 181 | 2 | 543 | 1.30 | 1.000 | 0.999 | 1.312 | 1.173–1.475 |
| 150 | 1 | 300 | 1.15 | 1.000 | 0.558 | 1.170 | 1.009–1.366 |
| 150 | 1 | 300 | 1.20 | 1.000 | 0.757 | 1.226 | 1.055–1.449 |
| 150 | 1 | 300 | 1.25 | 1.000 | 0.881 | 1.271 | 1.090–1.530 |
| 150 | 1 | 300 | 1.30 | 1.000 | 0.958 | 1.336 | 1.138–1.615 |
| 150 | 2 | 450 | 1.15 | 1.000 | 0.703 | 1.159 | 1.034–1.299 |
| 150 | 2 | 450 | 1.20 | 1.000 | 0.889 | 1.209 | 1.072–1.368 |
| 150 | 2 | 450 | 1.25 | 1.000 | 0.982 | 1.265 | 1.132–1.446 |
| 150 | 2 | 450 | 1.30 | 1.000 | 0.994 | 1.316 | 1.154–1.499 |
| 120 | 1 | 240 | 1.15 | 1.000 | 0.443 | 1.166 | 0.982–1.402 |
| 120 | 1 | 240 | 1.20 | 1.000 | 0.652 | 1.225 | 1.032–1.481 |
| 120 | 1 | 240 | 1.25 | 1.000 | 0.809 | 1.274 | 1.077–1.594 |
| 120 | 1 | 240 | 1.30 | 1.000 | 0.910 | 1.346 | 1.131–1.695 |
| 120 | 2 | 360 | 1.15 | 1.000 | 0.590 | 1.160 | 1.011–1.339 |
| 120 | 2 | 360 | 1.20 | 1.000 | 0.810 | 1.210 | 1.062–1.394 |
| 120 | 2 | 360 | 1.25 | 1.000 | 0.942 | 1.273 | 1.104–1.473 |
| 120 | 2 | 360 | 1.30 | 1.000 | 0.982 | 1.325 | 1.155–1.533 |
| 100 | 1 | 200 | 1.15 | 1.000 | 0.402 | 1.178 | 0.970–1.476 |
| 100 | 1 | 200 | 1.20 | 1.000 | 0.575 | 1.237 | 1.022–1.577 |
| 100 | 1 | 200 | 1.25 | 1.000 | 0.740 | 1.298 | 1.066–1.718 |
| 100 | 1 | 200 | 1.30 | 1.000 | 0.818 | 1.346 | 1.097–1.765 |
| 100 | 2 | 300 | 1.15 | 1.000 | 0.494 | 1.155 | 0.993–1.345 |
| 100 | 2 | 300 | 1.20 | 1.000 | 0.741 | 1.217 | 1.041–1.427 |
| 100 | 2 | 300 | 1.25 | 1.000 | 0.881 | 1.271 | 1.088–1.507 |
| 100 | 2 | 300 | 1.30 | 1.000 | 0.961 | 1.323 | 1.139–1.586 |
| 80 | 1 | 160 | 1.15 | 1.000 | 0.329 | 1.176 | 0.945–1.564 |
| 80 | 1 | 160 | 1.20 | 1.000 | 0.496 | 1.245 | 0.994–1.682 |
| 80 | 1 | 160 | 1.25 | 1.000 | 0.627 | 1.301 | 1.037–1.794 |
| 80 | 1 | 160 | 1.30 | 1.000 | 0.743 | 1.365 | 1.071–1.834 |
| 80 | 2 | 240 | 1.15 | 1.000 | 0.451 | 1.166 | 0.977–1.403 |
| 80 | 2 | 240 | 1.20 | 1.000 | 0.651 | 1.226 | 1.030–1.488 |
| 80 | 2 | 240 | 1.25 | 1.000 | 0.796 | 1.266 | 1.089–1.565 |
| 80 | 2 | 240 | 1.30 | 1.000 | 0.898 | 1.337 | 1.128–1.626 |

## Limitations frozen before WHI access

1. NHANES is an exposure-distribution proxy, not a substitute for WHI biomarker data.
2. The simulation sets nuisance outcome coefficients to zero; covariates are retained in the fitted model and their empirical correlation with MCOP is preserved through the NHANES proxy rows.
3. Matching is represented by conditional-logistic matched sets, but exact WHI matching-factor distributions are not yet known.
4. This is not a sample-size justification until WHI confirms analyzable urine samples and actual MCOP assay performance.

## Next gate

Confirm the number of CRC cases with sufficient prediagnostic urine volume, matched controls, sample timing and assay availability in WHI Query Builder/CCC. Then rerun this simulation with the observed WHI exposure variance and missingness assumptions before finalizing the sample request.