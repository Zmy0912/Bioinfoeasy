# CodonW2RSCUPlot

> 将 **codonW**（或 MEGA）的 RSCU 结果汇总为
> [**Codon Usage Visualization / RSCU-Plot**](https://pcg-lab.shinyapps.io/RSCU-Plot/)
> 在线应用所需的输入格式。

一个轻量的、带图形界面的 Python 小工具；无第三方依赖。

---

## 1. 功能简介

RSCU-Plot 在线应用只负责**可视化**（热图、聚类、PCA 等），不负责 RSCU 计算；
它要求上传一个**宽表 CSV**，结构为：

| 列 1 | 列 2 | 列 3 ... |
|------|------|----------|
| Amino Acid（氨基酸） | Codon（密码子） | 各物种的 RSCU 数值 |

而我们日常用 **codonW** 得到的 RSCU 是**每个物种一个文件**（或每条序列一个文件），
无法直接上传。本工具做的就是这件事：

```
每个物种的 RSCU 文件 (.blk / .out / .txt / .csv)
            │
            │  CodonW2RSCUPlot（本工具）
            ▼
          宽表 CSV
            │
            ▼
  RSCU-Plot (pcg-lab.shinyapps.io/RSCU-Plot)
```

特性：
- 🖱️ 图形界面（tkinter），点点按钮即可完成；
- 📦 零第三方依赖（仅使用 Python 标准库）；
- 🔄 自动识别 **MEGA `.blk`** 与 **codonW `.out` / `.txt` / `.csv`** 风格；
- 🧬 自动将 DNA 字母 `T` 标准化为 RNA 字母 `U`，统一跨格式汇总；
- 🪪 文件名（去扩展名）自动作为物种列名；
- 🧹 行顺序按标准密码子表分组（与官方模板一致：Phe→Leu→Ile→…→Gly）；
- ⚠️ 对空文件 / 解析失败的文件给出明确警告并保留空列；
- 🛡️ 后台线程运行，处理期间界面不卡死。

---

## 2. 环境要求

- **Python ≥ 3.7**
- 自带 `tkinter`（Windows / macOS 官方安装包默认包含；Linux 发行版可能需要
  `sudo apt install python3-tk`）
- 无需任何 `pip install`

---

## 3. 使用方法

### 3.1 准备输入

把每个物种的 RSCU 结果文件放到同一个文件夹里，例如：

```
my_rscu_results/
├── Species_A.blk
├── Species_B.blk
├── Species_C.out
└── Species_D.txt
```

文件名（去扩展名）会成为输出 CSV 里的物种列名。

### 3.2 运行

```bash
python codonw2rscuplot.py
```

弹出窗口，按提示操作：

1. 选择包含 RSCU 结果的文件夹；
2. 设置输出 CSV 路径（默认会在输入文件夹下生成 `rscu_input_for_plot.csv`）；
3. （可选）取消勾选「包含终止密码子」可剔除 TER；
4. 点击 **开始转换**，下方日志区会显示解析进度与最终结果。

### 3.3 上传至 RSCU-Plot

打开 <https://pcg-lab.shinyapps.io/RSCU-Plot/>，在 *Data Upload* 区域：

- 选择文件类型 **CSV**，分隔符选 **Comma**；
- 勾选 *Header row*；
- 上传刚才生成的 CSV；
- 切换到 *Data Preview* 检查无误后即可在 *Visualization* 中绘图。

---

## 4. 输入格式说明

### 4.1 MEGA 风格（`.blk`）

每行 4 列密码子，每列第一行有氨基酸缩写，第二行无：

```
Phe UUU  791 1.27 Ser UCU  469 1.69 Tyr UAU  656 1.61 Cys UGU  197 1.54
    UUC  450 0.73     UCC  274 0.99     UAC  159 0.39     UGC   59 0.46
Leu UUA  754 1.96     UCA  335 1.21 TER UAA   36 1.40 TER UGA   20 0.78
    UUG  476 1.24     UCG  161 0.58     UAG   21 0.82 Trp UGG  393 1.00
...
22178 codons in NC01_joined_from (used Universal Genetic code)
```

### 4.2 codonW 风格（`.out` / `.txt` / `.csv`）

每行一个密码子条目：`AminoAcid  Codon  Count  RSCU`：

```
Phe  UUU  791  1.27
Phe  UUC  450  0.73
Leu  UUA  754  1.96
...
```

### 4.3 通用规则

| 字段 | 允许写法 |
|------|---------|
| 氨基酸缩写 | `Phe` / `Leu` / `TER`（大小写不敏感，首字母大写即可） |
| 密码子 | `UUU`、`TTT` 均可（`T` 会被自动规范化为 `U`） |
| 计数 | 任意整数 |
| RSCU | 浮点数（如 `1.27`） |

未列出在 STANDARD_ORDER 中的密码子（如自定义物种里出现的非标准编码）会被自动
追加到 CSV 末尾，不会丢失。

---

## 5. 输出格式

UTF-8 with BOM 编码的 CSV（Excel 双击打开不会乱码），示例：

```csv
"Amino Acid","Codon","Species_A","Species_B","Species_C"
"Phe","UUU",1.27,1.31,1.25
"Phe","UUC",0.73,0.69,0.75
"Leu","UUA",1.96,2.00,1.92
"Leu","UUG",1.24,1.20,1.28
...
```

列 1、列 2 固定为 `Amino Acid`、`Codon`，其余列依次为各物种 RSCU 值。

---

## 6. 命令行（无界面）使用

GUI 仅是外壳，核心函数 `merge_rscu_files` 可被其他脚本调用：

```python
from codonw2rscuplot import merge_rscu_files

n_species, n_codons, empty = merge_rscu_files(
    input_dir="my_rscu_results",
    output_csv="rscu_input_for_plot.csv",
    include_ter=True,
)
print(n_species, n_codons, empty)
```

---

## 7. 常见问题

**Q：弹出 "所选文件夹中没有找到任何 RSCU 结果文件"。**
A：工具会自动识别 `.blk` / `.out` / `.txt` / `.csv` 四种扩展名。请确认：
- 文件确实放在该文件夹下；
- 文件扩展名属于上述四种之一；
- 文件非空（≥ 1 字节且包含可识别的密码子条目）。

**Q：某个物种对应的列全是空白。**
A：通常是因为该文件为空、文件编码异常、或内容不含任何密码子条目。日志区会列出
   该物种名，请重新生成对应物种的 RSCU 结果。

**Q：终止密码子要不要保留？**
A：RSCU-Plot 应用说明明确要求包含 Stop 行，建议保留；只有在做"仅编码区"分析
   时才勾选剔除。

**Q：密码子顺序能否自定义？**
A：当前固定为标准顺序。若有强需求，可编辑源码中的 `STANDARD_ORDER` 常量。

**Q：能否汇总 codonW 的 `.out` 中**每个基因**的 RSCU？**
A：不可以。本工具假设每个输入文件**已经是一个物种的 RSCU 汇总**。请先用
   codonW 的 *Correspondence analysis* / *Codon usage* 模块生成每个物种的
   汇总表，再交给本工具做横向拼接。

---

## 8. 致谢

- 可视化功能由
  [**Codon Usage Visualization (RSCU-Plot)**](https://pcg-lab.shinyapps.io/RSCU-Plot/)
  提供（开发：Aboozar Soorni, pcg-lab.iut.ac.ir）。
  使用该工具请按其要求引用相关论文：
  - <https://doi.org/10.1038/s41598-025-29339-x>
  - <https://doi.org/10.1007/s00425-025-04875-y>
  - <https://doi.org/10.1186/s12864-025-11729-0>
  - <https://doi.org/10.1007/s11033-025-10789-5>
- RSCU 计算可使用
  [**codonW**](http://codonw.sourceforge.net/)
  或 [**MEGA**](https://www.megasoftware.net/) 完成。

---

## 9. 目录结构

```
codonw2rscuplot/
├── codonw2rscuplot.py   # 主程序
└── README.md            # 本文档
```