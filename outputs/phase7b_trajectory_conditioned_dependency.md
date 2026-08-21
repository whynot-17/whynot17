# Phase 7B：OXA-R state-conditioned dependency mapping

已固定 6 条独立 OXA-R resistance trajectories；没有把它们合并成一个 consensus signature。

## 本轮已完成

- 为每条 trajectory 导出独立的 top-up/top-down directional signature；
- 导出 FAO、carnitine-entry、pyrimidine、UPR/ERAD相关、EMT、NRF2/redox、ferroptosis、ABC 等解释模块的 trajectory-level effect；
- 固定 DepMap CRISPR 的方向转换：原始 gene-effect 越负代表越依赖，报告中的 `vulnerability_rho = -rho` 越大代表 OXA-R-like state 越依赖该基因。

## 关键统计门槛

- 每条 trajectory 单独计算 Score；不把 HCT116、DLD1、HT29、LoVo 的轨迹强行合并。
- 每个 trajectory-gene 至少需要 15 个 CRC cell lines；严格 convergent vulnerability 要求至少 5/6 条轨迹方向一致、中位 vulnerability_rho ≥ 0.10、至少 2 条 trajectory 通过 permutation/bootstrap，且至少一条 trajectory 的 BH-q ≤ 0.10。
- 最终还需要 permutation/bootstrap 和药敏闭环；当前不能把转录投影结果称为功能依赖。

## 文件

- `phase7b_oxa_r_trajectory_signatures.csv`
- `phase7b_trajectory_signature_metadata.csv`
- `phase7b_module_trajectory_scores.csv`
- `phase7b_trajectory_dependency_manifest.json`
- `phase7b_crc_trajectory_scores.csv`
- `phase7b_trajectory_gene_dependency.csv`
- `phase7b_convergent_dependency_ranking.csv`
- `phase7b_mechanistic_dependency_audit.csv`

Top 50 genes per trajectory additionally receive 200 permutation/bootstrap resamples.

严格 convergent genes: 0；exploratory signals: 2。