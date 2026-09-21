# DINP–PPARG 2D/3D FEL plotting workflow

这套代码用于把 GROMACS `gmx sham` 生成的二维自由能网格，重现为：

- Origin contour 2D 图；
- Python/matplotlib 3D free-energy surface；
- 可选：从 Excel 的三列轨迹直接生成 `gmx sham` 输入和 XPM。

## 输入格式

绘图脚本统一读取制表符分隔的 XYZ 文件，必须包含三列：

```text
Complex_RMSD_nm    Protein_Rg_nm    Free_energy_kJ_mol
```

每一行是一个二维网格单元。未采样区域建议写为 `NaN`，不要把它们当作真实高能值。

## 从 Excel 生成 FEL 网格

Excel 的指定 sheet 前三列应为：

1. time；
2. complex RMSD；
3. Protein Rg。

在 WSL/GROMACS 环境中运行：

```bash
bash make_fel_energyxlsx.sh /mnt/c/path/to/energy.xlsx
```

可用环境变量覆盖参数：

```bash
SHEET=Sheet2 TSHAM_K=310 OUTDIR=/path/to/output \
  bash make_fel_energyxlsx.sh /mnt/c/path/to/energy.xlsx
```

脚本输出 `gibbs.xpm`。然后将 XPM 转换为 Origin/绘图脚本使用的 XYZ：

```bash
python export_origin_data.py gibbs.xpm energy.xlsx output_dir
```

该步骤同时输出：

- `gibbs_origin_fel_xyz_masked.txt`：推荐用于绘图，未采样/最高背景 bin 为 `NaN`；
- `gibbs_origin_fel_xyz_all.txt`：完整网格备份；
- `energyxlsx_origin_time_rmsd_rg.txt`：原始轨迹三列导出。

## 生成 2D Origin 图

需要 Windows OriginPro 和 `originpro` Python 包：

```powershell
python plot_fel_in_origin.py `
  path\to\gibbs_origin_fel_xyz_masked.txt `
  F:\origin\CONTOUR.otpu
```

如果 Origin 安装路径不是 `F:\origin`，必须把模板路径作为第二个参数传入。

2D 脚本的固定处理规则：

- 使用有限网格值；
- 对采样区域做 normalized Gaussian smoothing；
- 色阶上界为实际显示最大值加 `0.01`；
- 色标标题使用 `SPECTRUM1.title$`；
- Origin `layer.cmap.colorAbove=1`，避免色标顶部出现白色 Above-Max 方块；
- 输出 PDF、PNG 和可编辑 OPJU。

## 生成 3D 图

```powershell
python plot_fel_3d.py path\to\gibbs_origin_fel_xyz_masked.txt
```

3D 脚本输出 PDF、SVG 和 PNG。它使用局部线性插值和轻度平滑，仅用于显示，不修改原始 FEL 网格；主版本使用实际平滑后的自由能尺度。此前的对比增强版本若保留，文件名带有 `_display_scaled`，不应作为定量主结果。

## 可复现性边界

- 不要提交原始 Excel、XPM、MD trajectory 或大型二进制结果到代码目录；
- 新数据只需保持上述列名和单位即可重跑；
- 2D/3D 的平滑属于可视化步骤，不应解释为重新估计自由能；
- 论文中应同时保留原始 XPM/XYZ 和运行参数（温度、网格数、`nlevels`）。
