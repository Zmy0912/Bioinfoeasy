# COA Simple Plotter —— 仅正交轴主图的精简版 COA 绘图工具

本程序是 `..\coa_plotter\coa_plot_plotter.py` 的精简版：**只保留 COA 正交轴主图**（基因在 Axis 1 × Axis 2 上的分布散点图），去掉了 Scree 图与相关性面板；同时**完整保留多物种批量输入与组合绘图功能**（与 `pr2_plot_plotter.py` 的文件夹模式一致）。

## 1. 环境要求与运行

Python 3.8+，依赖 `numpy`、`matplotlib`。

```powershell
pip install numpy matplotlib
python coa_simple_plotter.py
```

## 2. 需要输入的文件

| 文件 | 必需？ | 作用 |
|------|--------|------|
| `genes.coa` | **必需** | CodonW COA 输出：每个基因在前几个轴上的坐标 |
| `eigen.coa` | 可选（建议） | 各轴解释的变异比例，用于坐标轴标签（如 `Axis 1 (9.57%)`） |

选中 `genes.coa` 后程序会**自动在同目录探测** `eigen.coa` 并填入。

文件格式示例：

**`genes.coa`**
```
label Axis1 Axis2 Axis3 Axis4
1_psbD_[33312:34373](+\0925\NC01.bl 0.043438 -0.022571 0.170300 -0.049392
2_rps2_[15496:16206](-\0925\NC01.bl 0.128910 -0.114170 0.078790 -0.121402
```

**`eigen.coa`**
```
01 6.8699E-002 0.0957 0.0957
02 5.9621E-002 0.0831 0.1788
```

这两个文件由 CodonW 的 **Correspondence Analysis (COA)** 功能生成（操作说明见 `..\coa_plotter\README.md` 第 2.1 节；`D:\零零碎碎\代码\NEW\0925` 中已有现成的 `genes.coa` 与 `eigen.coa`）。

## 3. 两种输入模式

### 3.1 单物种模式（Single species）

选择一个 `genes.coa`（可另选 `eigen.coa`），填入物种名，点击 **Plot Graph**，得到一张 7×7 英寸的正方形主图：

* X 轴 = Axis 1，Y 轴 = Axis 2，轴标签自动标注解释的变异比例；
* 原点处画灰色虚线参考线（正交轴交叉点）；
* 可按基因功能着色、标注离群基因（按距原点或数据中心的距离排序）。

### 3.2 文件夹模式（Folder，多物种组合图）

支持**两种摆放方式，可混用**：

**方式 A：多个 .coa 文件直接放在同一文件夹**（如 `D:\零零碎碎\代码\NEW\0925\gene`）：每个 `.coa` 基因坐标文件 = 一个物种子图，文件名（去扩展名）即子图标题；`eigen.coa`、`codon.coa` 等 CodonW 其它输出会被自动忽略。

**方式 B：每个物种一个子文件夹**（内含 `genes.coa`）：子文件夹名即子图标题。

* 每个子文件夹成为一个面板，组成一张网格组合图（列数可设，0 = 自动，默认每行 4 个）；
* 图例统一放在图右侧；各面板标题为子文件夹名；
* 勾选 **"Also export one image per file"** 后保存组合图时，会同时在 `individual_plots` 子文件夹中为每个物种单独导出一张图。

> 注意：若根目录本身直接含有 `genes.coa`（无子文件夹），则按单物种处理该目录。

## 4. 输出

* 图片：PDF / SVG / PNG / EPS，300 dpi；
* **轴统计表（"Export Axis Stats (CSV)"）**：
  * **Table 1**（每轴统计）：`species, n_genes, axis, eigenvalue, variation_explained_%, cumulative_%, total_inertia`，数值直接来自 CodonW 的 `eigen.coa`；
  * **Table 2**（每物种摘要）：`species, n_genes, n_axes, total_inertia, axis1–4_%, axis1+2_cumulative_%`，可直接用于正文 4.2 节。

> 轴解释比例**必须来自 CodonW 生成的 `eigen.coa`**。
* 与前一程序相同的样式设置（标题/轴标签/刻度字号与颜色、参考线、散点大小/透明度/颜色、图例、基因功能配色与显示开关）。

## 5. 目录结构

```
coa_simple_plotter\
  ├─ coa_simple_plotter.py     主程序（GUI，仅正交轴主图）
  ├─ README.md                 本说明文档
  └─ example_data\
      ├─ genes.coa             示例：单物种模式可直接选用
      ├─ eigen.coa
      └─ SpeciesA\             示例：文件夹模式（把多个物种子文件夹放在 example_data 下）
          ├─ genes.coa
          └─ eigen.coa
```

## 6. 与前两个程序的关系

| 程序 | 内容 |
|------|------|
| `pr2_plot_plotter.py` | PR2 图（G3/(G3+C3) × A3/(A3+T3)） |
| `..\coa_plotter\coa_plot_plotter.py` | 完整 COA 补充图：主图 + Scree 图 + 组成参数相关性面板 + 统计表导出（审稿意见要求的补充图） |
| 本程序 | **仅 COA 正交轴主图**，支持单物种与多物种批量组合绘图（适合正文使用的简洁版本） |
