# DINP axis：NHANES Phase 2A availability audit

本审计只回答能否进入 CRC 人群分析，不进行回归。分析框架为已有 Phase 2B 的 20 岁以上 CRC outcome frame 与 PHTHTE urine laboratory files 按 SEQN 连接。

NHANES comment code 按官方定义解释：0 = at or above detection limit，1 = below detection limit。

官方变量定义与 LOD 规则：[CDC NHANES PHTHTE_J codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PHTHTE_J.htm)。

## Summary

| Analyte | Variable | Cycles | Exposure + CRC outcome | CRC cases | Above LOD | Above LOD % |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| MiNP | URXMNP | 2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018 | 12127 | 76 | 3321 | 27.39% |
| MCOP | URXCOP | 2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018 | 12127 | 76 | 11933 | 98.40% |
| MONP | URXMONP | 2017-2018 | 1700 | 10 | 1435 | 84.41% |

## Interpretation boundary

Exposure + CRC outcome 是可用于当前 CRC outcome 分析的最低数据条件，不等同于足够的统计 power。MONP 若只覆盖单一 NHANES cycle，病例数必须先过稀疏性门槛，再决定是否运行回归。

跳过的文件和逐 cycle 明细见同目录 CSV/JSON。