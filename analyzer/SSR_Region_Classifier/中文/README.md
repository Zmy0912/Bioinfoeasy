# SSR 区域分类分析工具

一个图形化界面（GUI）程序，用于根据 MISA 输出与 GenBank (GB) 注释文件，**将每条 SSR 按其基因组坐标归类到 CDS / Intron / IGS 三类区域**，并输出可直接用 Excel / Origin / R / Python 绘图的统计表。

## 功能特点

- **直观 GUI**：基于 Tkinter，所见即所得，无需命令行操作。
- **自动物种匹配**：MISA 文件 ID 列与 GB 文件主名一一对应（NC01 ↔ NC01.gb）。
- **区域分类规则**：
  - **CDS** —— SSR 完全落在 CDS 外显子区间内（或中心点在外显子内）
  - **Intron** —— SSR 落在由 `gene` 与 `CDS` 差集推断的内含子中
  - **IGS** —— 其余均为基因间间隔区
- **多格式输出**：明细表、长表、宽表、SSR 类型 × 区域矩阵，便于直接绘图。

## 使用方法

### 1. 启动程序

双击 `启动程序.bat`，或在命令行中执行：

```bash
python ssr_region_classifier.py
```

### 2. 在 GUI 中操作

1. 点击「选择…」指定 **MISA 文件**（如 `0142.MISA`）。
2. 选择 **GB 文件夹**（与 MISA 同编号的子文件夹，内含 `NC01.gb`、`NC02.gb` …）。
3. 选择 **输出目录**（默认会按 MISA 文件名自动推荐）。
4. 点击「**开始分析**」。

### 3. 查看结果

完成后点击「打开输出目录」即可查看生成的文件。

## 输出文件

| 文件 | 说明 |
| --- | --- |
| `SSR_with_region.txt` | MISA 原内容 + `Region`（CDS / Intron / IGS）与 `RegionDetail`（所在区间）列 |
| `SSR_region_summary.tsv` | 长表：`Species, Region, Count, Percent`，适合 ggplot 直接画堆叠柱状图 |
| `SSR_region_summary_wide.tsv` | 宽表：`Species, CDS, Intron, IGS, Total` |
| `SSR_type_by_region.tsv` | 行=区域，列=SSR 类型（p1/p2/.../c），便于画分组柱状图 |
| `log.txt` | 运行日志 |

## 输入格式要求

### MISA 文件
- 制表符分隔（`.misa` 或 `.txt`），表头包含 `ID、SSR nr.、SSR type、SSR、size、start、end`。
- 第 1 列 `ID`（如 `NC01`）与 GB 文件主名（`NC01.gb`）一一对应。

### GB 文件夹
- 每个 GB 文件需包含 FEATURES 表的 `gene` 与 `CDS` 特征。
- 程序会自动识别 `complement(join(...))` 等复杂位置。

## 依赖

- Python 3.8+
- 仅使用标准库（`tkinter`、`csv`、`re`、`threading` 等），**无需安装额外包**。

## 目录结构示例

```
0919\
├─ 0142.MISA                ← 输入
├─ 0142\                    ← GB 文件夹
│   ├─ NC01.gb
│   ├─ NC02.gb
│   └─ ...
└─ SSR_Region_Classifier\   ← 本程序
    ├─ ssr_region_classifier.py
    └─ 启动程序.bat
```