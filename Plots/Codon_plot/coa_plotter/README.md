# COA Plot Plotter —— CodonW 对应分析（COA）补充图绘制工具

本工具用于响应审稿意见，绘制 **密码子使用模式对应分析（Correspondence Analysis, COA）** 的补充图（Supplementary Figure）与相关性统计表。

审稿意见要点：

1. 必须使用 **CodonW** 执行 COA（本文已完成）；
2. 需说明 **COA 产生一系列正交轴：第一轴解释最大的变异比例，后续各轴解释的变异依次递减**；
3. 补充图应展示 **基因在第一、二轴上的分布，并与核苷酸组成参数（GC3s、GC、GC12 等）做相关性分析**，在正文 4.2 节讨论。

程序界面与操作方式与上一级目录的 `pr2_plot_plotter.py` 保持一致（单文件 / 文件夹输入模式、按基因功能着色、离群基因标注、样式设置、图片导出）。

## 1. 环境要求与运行

Python 3.8+，依赖 `numpy`、`matplotlib`（**不需要** scipy，P 值由内置 t 分布精确计算）。

```powershell
pip install numpy matplotlib
python coa_plot_plotter.py
```

## 2. 需要输入的文件（重点）

程序读取 **CodonW 的 COA 输出文件**，共三类：

| 文件 | 必需？ | 作用 |
|------|--------|------|
| `genes.coa` | **必需** | 每个基因在 COA 前几轴上的坐标（散点数据来源） |
| `eigen.coa` | 可选（强烈建议） | 各轴特征值与解释的变异比例（用于轴标签 `Axis 1 (9.57%)` 和 Scree 图） |
| `coa_raw` 或 CodonW 批量输出 `.out` | 可选（建议） | 每个基因的核苷酸组成参数，用于相关性面板 |

> 选中 `genes.coa` 后，程序会**自动在同目录探测** `eigen.coa`、`coa_raw` 或非空 `.out` 并填入，一般无需手动逐个选择。

### 2.1 这些文件如何生成（CodonW 操作）

若尚未生成，可用 CodonW（上一级 `0925` 文件夹中已有 `CodonW.exe`）：

1. 打开 CodonW，载入 CDS 序列文件（如 `NC01.fasta`）；
2. 选择 **Codon usage → Correspondence Analysis (COA)**，轴数 ≥4（默认 4）；
3. 同时勾选 bulk output，得到 `.out` 文件（含 GC3s、GC、GC12、CAI、ENC 等列）；
4. 运行结束后生成：`genes.coa`、`eigen.coa`、`codon.coa`、`coa_raw`、`summary.coa`、`cusort.coa` 及 `.out`。

> 您的 `D:\零零碎碎\代码\NEW\0925` 文件夹中已存在上述全部文件（`summary.coa` 显示输入为 `NC01.fasta`，77 个基因，58 个轴），可直接使用。

### 2.2 文件格式

**`genes.coa`**（首行表头，之后每行 = 基因标签 + 各轴坐标）：

```
label Axis1 Axis2 Axis3 Axis4
1_psbD_[33312:34373](+\0925\NC01.bl 0.043438 -0.022571 0.170300 -0.049392
```

**`eigen.coa`**（每行 = 轴编号、特征值、解释比例、累计比例）：

```
01 6.8699E-002 0.0957 0.0957
02 5.9621E-002 0.0831 0.1788
```

**`coa_raw`**（每行 = 基因标签 + 59 个同义密码子使用次数，CodonW 标准顺序）：

```
1_psbD_[33312:34373](+\0925\NC01.blk 21	6	7	3	18	5	...
```

程序自动由 `coa_raw` 计算每个基因的 `GC3s`、`GC`、`GC12`、`L_aa`（0–1 分数，与 CodonW 约定一致）；若提供 `.out`，则直接读取表头中的 `GC3s / GC / GC12 / CAI / ENC / Fop / L_aa` 列。

### 2.3 文件夹模式（多个物种）

选择一个根目录，其中**每个子文件夹**放一套 CodonW 输出（`genes.coa` + `eigen.coa` + `coa_raw` 或 `.out`）：

```
species_folder\
  ├─ SpeciesA\  genes.coa  eigen.coa  coa_raw
  ├─ SpeciesB\  genes.coa  eigen.coa  coa_raw
  └─ ...
```

文件夹模式输出所有物种的 COA 散点网格图（共享图例）；统计表按物种分别导出并含汇总（pooled）行。

## 3. 输出内容

### 3.1 补充图（单物种模式，"Plot Graph"）

* **左上大图**：基因在 **Axis 1 × Axis 2** 上的分布；轴标签自动标注解释比例；原点画灰色虚线参考线；可按基因功能着色并标注离群基因。
* **右上 Scree 图**：各轴解释的变异比例（柱）与累计比例（折线），直观展示"第一轴解释最大、后续轴依次递减"。
* **下两行相关性面板**：`Axis 1 / Axis 2` 分别对 `GC3s / GC / GC12`（可选 CAI、ENC、Fop、L_aa）的散点图，含红色回归线及 Pearson r、P、Spearman ρ 注释。

可导出 PDF / SVG / PNG / EPS（300 dpi）。

### 3.2 统计表（"Export Stats (CSV)"）

CSV 列：`dataset, axis, parameter, n, pearson_r, pearson_P, spearman_rho, spearman_P`。
示例（NC01 实测）：`NC01,2,GC,77,-0.3273,0.0037,-0.3274,0.0037`。

## 4. 可直接用于回复审稿人 / 4.2 节的建议文字（英文模板）

> Correspondence analysis (COA) of relative synonymous codon usage was performed with CodonW. COA projects the codon usage of each gene onto a series of mutually orthogonal axes, which are ordered by the amount of variation they explain: the first axis accounts for the largest proportion of the total inertia, and each subsequent axis explains a progressively smaller fraction. In *species*, axis 1 explained X.XX% and axis 2 explained X.XX% of the total variation (Fig. S x). Most genes clustered near the origin of the plane formed by the first two axes, whereas a few genes occupied extreme positions (Fig. S xA). Because the axes are orthogonal, each axis describes an independent direction of codon-usage variation. The correlations between gene positions on the first two axes and nucleotide composition parameters are summarized in Table S x. The significant positive correlation between axis 1 and the GC content at synonymous third codon positions (GC3s; Pearson r = X.XXX, P < 0.001) indicates that the major axis of codon-usage variation is largely determined by mutational pressure acting on nucleotide composition, while the weaker correlation of axis 2 with GC3s suggests that a smaller, independent component of codon-usage variation may reflect additional influences such as gene expression level or photosystem-related functional constraints. These results are discussed in section 4.2.

使用时把图中的实际数值（轴解释比例、r、P、离群基因名）填入即可。示例数据（`example_data`，即 NC01）的实际结果：Axis 1 = 9.57%，Axis 2 = 8.31%（累计 17.88%）；Axis 1 vs GC3s：r = 0.126，P = 0.275；Axis 2 vs GC：r = −0.327，P = 0.004；Axis 2 vs GC12：r = −0.408，P < 0.001。**请以您实际数据运行后的输出为准。**

## 5. 目录结构

```
coa_plotter\
  ├─ coa_plot_plotter.py   主程序（GUI）
  ├─ README.md             本说明文档
  └─ example_data\
      ├─ genes.coa         示例：NC01 基因轴坐标
      ├─ eigen.coa         示例：NC01 轴特征值
      └─ coa_raw           示例：NC01 原始密码子使用
```

## 6. 常见问题

* **相关性面板空白**：未提供组成参数文件，或未勾选参数。请选择 `coa_raw` / `.out`，并在参数区勾选 `GC3s`、`GC`、`GC12`。
* **基因无法与组成参数匹配**：程序按基因编号（标签开头的数字）匹配；请确保 `genes.coa` 与 `coa_raw`/`.out` 来自**同一次** CodonW 运行（控制台会显示 `Composition match: n/m genes`）。
* **`NC01.out` 为空**：您的目录中该文件为空，程序会自动回退到 `coa_raw` 计算参数，无需处理。
* **数值 P 值计算**：使用 t 分布 + 不完全 Beta 函数精确计算，双尾；r 的自由度为 n−2。
