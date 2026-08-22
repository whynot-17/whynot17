# CELLxGENE Census WSL2 环境

## 已安装位置

- WSL distribution: `Ubuntu-24.04`
- WSL distribution storage: `D:\WSL\Ubuntu`
- WSL version: 2
- Linux virtual environment: `/opt/cellxgene-census`
- Linux Python: 3.12.3

Windows 主 Python 环境没有安装旧版 TileDB-SOMA；Census 依赖只放在这个 Linux 虚拟环境中。

## 固定版本

- `cellxgene-census==1.17.0`
- `tiledbsoma==1.17.1`
- Reproducibility Census release: `2025-11-08`

对应的 pip requirements 见 `work/requirements_cellxgene_census_wsl.txt`。

## 验证结果

已通过：

```bash
/opt/cellxgene-census/bin/python -c \
  "import tiledbsoma, cellxgene_census; \
   print(tiledbsoma.__version__, cellxgene_census.__version__)"
```

并成功执行：

```python
import cellxgene_census

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    print(list(census["census_data"].keys()))
```

返回的 organism collections 为：`callithrix_jacchus`、`homo_sapiens`、`macaca_mulatta`、`mus_musculus`、`pan_troglodytes`。

## 使用方式

从 Windows PowerShell 调用：

```powershell
wsl -d Ubuntu-24.04 --user root -- /opt/cellxgene-census/bin/python /path/to/your/script.py
```

生产分析必须显式固定 `census_version`，并在细胞查询中加入：

```text
is_primary_data == True
```
