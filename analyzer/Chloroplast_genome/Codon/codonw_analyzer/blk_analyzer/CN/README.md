# CodonW密码子偏好性统计分析工具

一个用于处理CodonW输出结果的图形化工具，可生成热图可直接使用的表格。

## 功能特点

- 📊 批量解析CodonW生成的`.blk`文件
- 🔬 完整统计所有61个密码子的偏好性结果
- 🏷️ 自动标注密码子对应的氨基酸
- 📈 生成可直接用于热图绘制的CSV表格
- 🎨 图形化界面，操作简单直观
- 🔍 数据质量检测与清理建议
- ⚡ 支持按氨基酸分组排序

## 安装依赖

```bash
pip install pandas
```

## 使用方法

### 1. 运行程序

```bash
python codonw_analyzer.py
```

### 2. 界面操作

- **输入文件夹**：选择包含CodonW结果文件的文件夹（.blk文件）
- **输出文件夹**：选择结果保存位置（默认与输入文件夹相同）
- **统计指标**：选择RSCU（相对同义密码子使用度）或Count（密码子计数）
- **按氨基酸排序**：勾选后，密码子列将按氨基酸分组排列
- 点击"开始分析"按钮执行分析

## 输出文件

### 热图表格

- **未排序**：`codonw_heatmap_rscu.csv` 或 `codonw_heatmap_count.csv`
- **已排序**：`codonw_heatmap_sorted_rscu.csv` 或 `codonw_heatmap_sorted_count.csv`

### 数据质量报告

- `codonw_data_quality_report.txt`：包含数据质量分析和清理建议

## 热图表格格式

| 列 | 说明 |
|---|---|
| Species | 物种名称 |
| UUU(Phe) | 密码子(氨基酸) |
| UUC(Phe) | 密码子(氨基酸) |
| ... | 其他密码子 |
| 数据 | RSCU值或密码子计数 |

## 氨基酸排序说明

勾选"按氨基酸排序"后，密码子列将按照以下标准生物化学顺序排列：

```
Ala → Arg → Asn → Asp → Cys → Gln → Glu → Gly → His → Ile → 
Leu → Lys → Met → Phe → Pro → Ser → Thr → Trp → Tyr → Val
```

同一氨基酸内的多个密码子按字母顺序排列。

**优势**：
- 便于在热图中观察同一氨基酸的不同密码子偏好性
- 氨基酸分组排列使结果更直观
- 便于识别密码子使用模式

## 数据质量报告内容

1. **缺失密码子检测**：标记每个物种缺失的密码子
2. **零值密码子**：列出计数为0的密码子（可能影响统计分析）
3. **RSCU异常值**：识别RSCU值 > 2.0 的密码子（可能存在异常）
4. **清理建议**：提供数据处理和清理的建议

## 示例

假设CodonW结果文件夹包含以下文件：
```
codonW结果/
├── cdsNC01.blk
├── cdsNC02.blk
└── ...
```

运行工具后，将生成：
```
codonW结果/
├── cdsNC01.blk
├── cdsNC02.blk
├── codonw_heatmap_rscu.csv          # 热图表格（RSCU值）
├── codonw_heatmap_sorted_rscu.csv   # 热图表格（按氨基酸排序）
└── codonw_data_quality_report.txt   # 数据质量报告
```

## 热图绘制示例

使用Python绘制热图：

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 读取数据
df = pd.read_csv('codonw_heatmap_sorted_rscu.csv')

# 设置物种名称为索引
df.set_index('Species', inplace=True)

# 绘制热图
plt.figure(figsize=(16, 10))
sns.heatmap(df, cmap='viridis', center=1.0, 
            cbar_kws={'label': 'RSCU值'})
plt.title('密码子偏好性热图')
plt.xlabel('密码子(氨基酸)')
plt.ylabel('物种')
plt.tight_layout()
plt.savefig('codonw_heatmap.png', dpi=300)
plt.show()
```

## 注意事项

1. 确保输入文件夹中只包含`.blk`格式的CodonW结果文件
2. 缺失密码子的数据可能需要手动补充或排除
3. 计数为0的密码子在热图分析中可能导致偏差，根据分析需求决定是否保留
4. RSCU异常值建议检查原始数据或分析流程

## 技术支持

如有问题或建议，请查看数据质量报告或检查输入数据格式。
