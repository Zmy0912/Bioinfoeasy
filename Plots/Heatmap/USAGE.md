# Heatmap Visualization Tool - Detailed User Guide

This guide provides detailed instructions and best practices for using the heatmap visualization tool.

## Table of Contents

1. [Installation Guide](#installation-guide)
2. [Data Preparation](#data-preparation)
3. [Basic Usage](#basic-usage)
4. [Advanced Configuration](#advanced-configuration)
5. [Color Scheme Details](#color-scheme-details)
6. [Cluster Analysis](#cluster-analysis)
7. [Output Optimization](#output-optimization)
8. [Common Scenarios](#common-scenarios)
9. [Troubleshooting](#troubleshooting)

---

## Installation Guide

### System Requirements

- Python 3.7+
- Operating System: Windows / macOS / Linux

### Installation Steps

#### 1. Install Python Dependencies

Install using pip:

```bash
pip install pandas numpy matplotlib seaborn scipy
```

Or use the provided requirements.txt:

```bash
pip install -r requirements.txt
```

#### 2. Verify Installation

Run the test program:

```python
from heatmap_visualizer import HeatmapVisualizer
import matplotlib
import seaborn
import pandas
import numpy
import scipy

print(f"pandas: {pandas.__version__}")
print(f"numpy: {numpy.__version__}")
print(f"matplotlib: {matplotlib.__version__}")
print(f"seaborn: {seaborn.__version__}")
print(f"scipy: {scipy.__version__}")
```

#### 3. Chinese Font Configuration (Windows Users)

Ensure one of the following Chinese fonts is installed on your system:
- SimHei (Blackbody)
- Microsoft YaHei (Microsoft YaHei)
- Arial Unicode MS

The program will automatically detect and use these fonts.

---

## Data Preparation

### CSV File Format

#### Format Requirements

1. **First column**: Row index (must be unique identifiers, e.g., sample names, gene names)
2. **First row**: Column names (must be unique identifiers, e.g., gene names, time points)
3. **Data area**: All cells must be numeric
4. **Row/column consistency**: Each row must have the same number of columns

#### Correct Example

```csv
,Gene1,Gene2,Gene3,Gene4,Gene5
Sample1,1.23,2.45,3.67,4.89,5.12
Sample2,2.34,3.56,4.78,5.90,6.23
Sample3,3.45,4.67,5.89,6.12,7.34
Sample4,4.56,5.78,6.90,7.23,8.45
Sample5,5.67,6.89,7.12,8.34,9.56
```

#### Incorrect Examples

❌ **Missing first row/column names**:
```csv
Sample1,1.23,2.45,3.67,4.89,5.12
Sample2,2.34,3.56,4.78,5.90,6.23
```

❌ **Inconsistent column counts**:
```csv
,Gene1,Gene2,Gene3
Sample1,1.23,2.45,3.67
Sample2,2.34,3.56     ← Missing data
```

❌ **Non-numeric data**:
```csv
,Gene1,Gene2,Gene3
Sample1,1.23,2.45,N/A  ← N/A is not numeric
Sample2,2.34,3.56,4.78
```

### Excel File Format

Excel files have the same format requirements as CSV. Advantages of using Excel:
- More intuitive data editing
- Automatic format preservation
- Supports Chinese file names

### Data Preprocessing Recommendations

#### 1. Missing Value Handling

Handle missing values before importing:

```python
import pandas as pd

# Read data
df = pd.read_csv("your_data.csv", index_col=0)

# Method 1: Fill with mean
df = df.fillna(df.mean())

# Method 2: Fill with median
df = df.fillna(df.median())

# Method 3: Fill with 0
df = df.fillna(0)

# Method 4: Drop rows/columns with missing values
df = df.dropna()
# or
df = df.dropna(axis=1)

# Save processed data
df.to_csv("processed_data.csv")
```

#### 2. Data Normalization (Optional)

For data with different scales, normalization is recommended:

```python
from sklearn.preprocessing import StandardScaler

# Z-score normalization
scaler = StandardScaler()
df_normalized = pd.DataFrame(
    scaler.fit_transform(df),
    index=df.index,
    columns=df.columns
)
df_normalized.to_csv("normalized_data.csv")
```

#### 3. Log Transformation (Optional)

For data with large ranges (e.g., gene expression data):

```python
import numpy as np

# Log transformation (avoid log(0))
df_log = np.log2(df + 1)  # or np.log1p(df)
df_log.to_csv("log_transformed_data.csv")
```

---

## Basic Usage

### Simplest Usage

```python
from heatmap_visualizer import HeatmapVisualizer

# Step 1: Create visualizer
visualizer = HeatmapVisualizer()

# Step 2: Generate sample data (or load data)
visualizer.generate_sample_data(n_rows=10, n_cols=8)

# Step 3: Plot heatmap
visualizer.plot_heatmap(save_path="simple_heatmap.png")

# Step 4: Display image (optional)
visualizer.show()
```

### Using Your Own Data File

```python
from heatmap_visualizer import HeatmapVisualizer

# Create and load data
visualizer = HeatmapVisualizer(data_file="my_data.csv")
visualizer.load_data()

# Plot and save
visualizer.plot_heatmap(save_path="my_heatmap.png")
visualizer.show()
```

### Step-by-Step (More Flexible)

```python
from heatmap_visualizer import HeatmapVisualizer

# Step 1: Create visualizer
visualizer = HeatmapVisualizer(
    title="My Heatmap",
    xlabel="Genes",
    ylabel="Samples"
)

# Step 2: Load data
data = visualizer.load_data("my_data.csv")

# Step 3: View data information
print(f"Data shape: {data.shape}")
print(f"Data range: {data.min().min():.2f} - {data.max().max():.2f}")

# Step 4: Plot heatmap (can plot multiple configurations)
visualizer.plot_heatmap(save_path="heatmap_v1.png")

# Modify configuration and replot
visualizer.cmap = "coolwarm"
visualizer.plot_heatmap(save_path="heatmap_v2.png")
```

---

## Advanced Configuration

### Font Size Adjustment

#### Example: Large Fonts (For Presentations)

```python
visualizer = HeatmapVisualizer(
    title="Large Title Heatmap",
    xlabel="Column Label",
    ylabel="Row Label",
    title_size=24,      # Extra large title
    label_size=18,      # Large axis labels
    tick_size=14        # Large tick labels
)
```

#### Example: Small Fonts (For High-Density Data)

```python
visualizer = HeatmapVisualizer(
    title="Small Title Heatmap",
    xlabel="Columns",
    ylabel="Rows",
    title_size=12,      # Small title
    label_size=10,      # Small axis labels
    tick_size=6         # Small tick labels
)
```

### Color Bar Adjustment

#### Example: Wide Color Bar

```python
visualizer = HeatmapVisualizer(
    title="Wide Color Bar Heatmap",
    cbar_width_ratio=0.3,     # 30% width
    cbar_height_ratio=1.0     # 100% height
)
```

#### Example: Narrow Color Bar

```python
visualizer = HeatmapVisualizer(
    title="Narrow Color Bar Heatmap",
    cbar_width_ratio=0.05,    # 5% width
    cbar_height_ratio=0.5     # 50% height
)
```

#### Example: Medium Centered Bar

```python
visualizer = HeatmapVisualizer(
    title="Medium Color Bar Heatmap",
    cbar_width_ratio=0.15,    # 15% width (default)
    cbar_height_ratio=0.8     # 80% height
)
```

### Color Range Control

#### Example: Fixed Range (0-1)

```python
visualizer = HeatmapVisualizer(
    title="Fixed Range Heatmap",
    vmin=0,
    vmax=1,
    cmap="coolwarm"
)
```

#### Example: Zero-Centered Range (-10 to 10)

```python
visualizer = HeatmapVisualizer(
    title="Centered Heatmap",
    vmin=-10,
    vmax=10,
    cmap="RdBu"
)
```

#### Example: Automatic Range (Recommended)

```python
visualizer = HeatmapVisualizer(
    title="Auto Range Heatmap",
    vmin=None,  # Auto
    vmax=None   # Auto
)
```

### Output Resolution Adjustment

#### Example: Publication Quality (300 DPI)

```python
visualizer = HeatmapVisualizer(
    title="Publication Heatmap",
    dpi=300
)
```

#### Example: Poster Quality (150 DPI)

```python
visualizer = HeatmapVisualizer(
    title="Poster Heatmap",
    dpi=150
)
```

#### Example: Quick Preview (100 DPI)

```python
visualizer = HeatmapVisualizer(
    title="Preview Heatmap",
    dpi=100
)
```

---

## Color Scheme Details

### Sequential Colormaps

Suitable for unidirectional data such as intensity, concentration, expression level, etc.

| Scheme | Colors | Use Case | Example |
|--------|--------|----------|---------|
| `viridis` | Purple→Blue→Green→Yellow | Default, colorblind-friendly | Gene expression |
| `plasma` | Blue→Purple→Red→Yellow | High contrast | Temperature data |
| `inferno` | Black→Red→Yellow | High contrast, dark background | Intensity data |
| `magma` | Black→Purple→Yellow | High contrast | Density data |
| `cividis` | Blue→Green→Yellow | Colorblind-friendly | Public data |
| `rocket` | Black→Red→Orange | Unidirectional increase | Count data |

```python
# Sequential colormap example
visualizer = HeatmapVisualizer(
    title="Sequential - viridis",
    cmap="viridis"  # or plasma, inferno, magma, etc.
)
```

### Diverging Colormaps

Suitable for data with positive/negative values or meaningful center points, such as correlations, differences, etc.

| Scheme | Colors | Use Case | Example |
|--------|--------|----------|---------|
| `coolwarm` | Blue→White→Red | Correlation, difference | Correlation coefficient |
| `RdBu` | Red→White→Blue | Positive-negative comparison | Fold-change |
| `RdYlBu` | Red→Yellow→Blue | Multi-level classification | Risk assessment |
| `PiYG` | Purple→White→Green | Balanced comparison | Score analysis |
| `PRGn` | Purple→White→Green | Balanced comparison | Rating data |
| `Spectral` | Red→Yellow→Green→Blue | Multi-spectrum | Classification data |

```python
# Diverging colormap example (important: usually need to set vmin and vmax)
visualizer = HeatmapVisualizer(
    title="Diverging - coolwarm",
    cmap="coolwarm",
    vmin=-1,    # Set center value
    vmax=1      # Set center value
)
```

### Monochrome Colormaps

Suitable for mainly positive values where intensity needs to be emphasized.

| Scheme | Colors | Use Case | Example |
|--------|--------|----------|---------|
| `Greens` | Light→Dark Green | Positive, growth data | Growth curves |
| `Blues` | Light→Dark Blue | Positive, temperature data | Temperature distribution |
| `Reds` | Light→Dark Red | Positive, risk data | Risk assessment |
| `Greys` | Light→Dark Grey | Grayscale data | Image processing |

```python
# Monochrome colormap example
visualizer = HeatmapVisualizer(
    title="Monochrome - Greens",
    cmap="Greens"
)
```

### Color Scheme Selection Guide

```
Data Type → Recommended Colormap
───────────────────────────
General data → viridis
Correlation → coolwarm (set vmin=-1, vmax=1)
Differential analysis → RdBu (set symmetric range)
Expression → plasma or inferno
Positive data → Greens or Blues
Classification → Spectral
Colorblind-friendly → cividis or viridis
```

---

## Cluster Analysis

### Clustering Explanation

Cluster analysis automatically arranges similar rows or columns together, making it easier to discover patterns in the data.

- **Row clustering**: Reorders rows based on their similarity
- **Column clustering**: Reorders columns based on their similarity
- **Bidirectional clustering**: Clusters both rows and columns simultaneously

### Clustering Algorithm

This tool uses hierarchical clustering:
- Clustering method: average (average linkage)
- Distance metric: Euclidean distance (default)

### Clustering Configuration

#### Example 1: Bidirectional Clustering (Default)

```python
visualizer = HeatmapVisualizer(
    title="Bidirectional Clustering Heatmap",
    cluster_rows=True,   # Row clustering
    cluster_cols=True    # Column clustering
)
```

#### Example 2: Row Clustering Only

```python
visualizer = HeatmapVisualizer(
    title="Row Clustering Heatmap",
    cluster_rows=True,   # Row clustering
    cluster_cols=False   # No column clustering
)
```

#### Example 3: Column Clustering Only

```python
visualizer = HeatmapVisualizer(
    title="Column Clustering Heatmap",
    cluster_rows=False,   # No row clustering
    cluster_cols=True    # Column clustering
)
```

#### Example 4: No Clustering

```python
visualizer = HeatmapVisualizer(
    title="No Clustering Heatmap",
    cluster_rows=False,  # No row clustering
    cluster_cols=False   # No column clustering
)
```

### Dendrogram Explanation

When clustering is enabled, the heatmap will display dendrograms:

- **Column dendrogram** (top): Shows similarity relationships between columns
- **Row dendrogram** (left): Shows similarity relationships between rows

The height of the dendrogram indicates the distance between clusters - greater height means larger differences.

### Clustering Usage Recommendations

```
Data Feature → Clustering Strategy
─────────────────────────────
Discover patterns → Bidirectional clustering (recommended)
Columns have order → Row clustering only
Rows have order → Column clustering only
Keep original order → No clustering
Large sample size → Consider disabling clustering to improve speed
```

---

## Output Optimization

### Image Size Control

Image size is automatically calculated based on data dimensions:

```python
# Width = max(8, number of columns × 0.6 + 2)
# Height = max(6, number of rows × 0.6 + 1.5)
```

### Quality Optimization

#### Publication Quality

```python
visualizer = HeatmapVisualizer(
    title="Publication Figure",
    dpi=300,
    title_size=16,
    label_size=12,
    tick_size=8
)
```

#### Presentation Quality

```python
visualizer = HeatmapVisualizer(
    title="Presentation Figure",
    dpi=150,
    title_size=20,
    label_size=14,
    tick_size=10
)
```

### Batch Generate Multiple Configurations

```python
from heatmap_visualizer import HeatmapVisualizer

# Load data once
base_visualizer = HeatmapVisualizer(data_file="my_data.csv")
base_visualizer.load_data()
data = base_visualizer.data

# Define multiple color schemes
cmaps = ["viridis", "coolwarm", "plasma", "RdBu"]

# Batch generate
for cmap in cmaps:
    visualizer = HeatmapVisualizer(
        title=f"Color Scheme: {cmap}",
        xlabel="Genes",
        ylabel="Samples",
        cmap=cmap
    )
    visualizer.data = data  # Reuse data
    visualizer.plot_heatmap(save_path=f"heatmap_{cmap}.png")
    print(f"Generated: heatmap_{cmap}.png")
```

---

## Common Scenarios

### Scenario 1: Gene Expression Analysis

```python
from heatmap_visualizer import HeatmapVisualizer

visualizer = HeatmapVisualizer(
    data_file="gene_expression.csv",
    title="Gene Expression Heatmap Analysis",
    xlabel="Genes",
    ylabel="Samples",
    cmap="plasma",
    cluster_rows=True,
    cluster_cols=True,
    title_size=16,
    label_size=12,
    tick_size=8,
    dpi=300
)

visualizer.load_data()
visualizer.plot_heatmap(save_path="gene_expression.png")
visualizer.show()
```

### Scenario 2: Correlation Analysis

```python
from heatmap_visualizer import HeatmapVisualizer

visualizer = HeatmapVisualizer(
    data_file="correlation_matrix.csv",
    title="Variable Correlation Analysis",
    xlabel="Variables",
    ylabel="Variables",
    cmap="coolwarm",
    vmin=-1,
    vmax=1,
    cluster_rows=True,
    cluster_cols=True,
    title_size=18,
    label_size=14,
    tick_size=10,
    dpi=300
)

visualizer.load_data()
visualizer.plot_heatmap(save_path="correlation.png")
visualizer.show()
```

### Scenario 3: Time Series Data

```python
from heatmap_visualizer import HeatmapVisualizer

# Time series usually keeps column order, no column clustering
visualizer = HeatmapVisualizer(
    data_file="time_series.csv",
    title="Time Series Analysis",
    xlabel="Time Points",
    ylabel="Samples",
    cmap="viridis",
    cluster_rows=True,    # Row clustering
    cluster_cols=False,   # No column clustering (keep time order)
    title_size=16,
    label_size=12,
    tick_size=8,
    dpi=300
)

visualizer.load_data()
visualizer.plot_heatmap(save_path="time_series.png")
visualizer.show()
```

### Scenario 4: Differential Analysis (Fold Change)

```python
from heatmap_visualizer import HeatmapVisualizer

visualizer = HeatmapVisualizer(
    data_file="fold_change.csv",
    title="Differential Expression Analysis (Fold Change)",
    xlabel="Genes",
    ylabel="Samples",
    cmap="RdBu",
    vmin=-2,   # Down-regulated 2-fold
    vmax=2,    # Up-regulated 2-fold
    cluster_rows=True,
    cluster_cols=True,
    title_size=16,
    label_size=12,
    tick_size=8,
    dpi=300
)

visualizer.load_data()
visualizer.plot_heatmap(save_path="fold_change.png")
visualizer.show()
```

### Scenario 5: Quick Prototype (Using Sample Data)

```python
from heatmap_visualizer import HeatmapVisualizer

# Quickly generate sample data for testing
visualizer = HeatmapVisualizer(
    title="Prototype Test",
    xlabel="Test Columns",
    ylabel="Test Rows",
    cmap="viridis"
)

visualizer.generate_sample_data(n_rows=15, n_cols=12)
visualizer.plot_heatmap(save_path="prototype.png")
visualizer.show()
```

### Scenario 6: Multi-Configuration Comparison

```python
from heatmap_visualizer import HeatmapVisualizer

# Load data
base = HeatmapVisualizer(data_file="data.csv")
base.load_data()
data = base.data

# Configuration 1: Standard
v1 = HeatmapVisualizer(
    title="Standard Configuration",
    cmap="viridis",
    cluster_rows=True,
    cluster_cols=True
)
v1.data = data
v1.plot_heatmap(save_path="config1_standard.png")

# Configuration 2: High contrast
v2 = HeatmapVisualizer(
    title="High Contrast",
    cmap="plasma",
    cluster_rows=True,
    cluster_cols=True
)
v2.data = data
v2.plot_heatmap(save_path="config2_contrast.png")

# Configuration 3: Correlation colormap
v3 = HeatmapVisualizer(
    title="Correlation Colormap",
    cmap="coolwarm",
    cluster_rows=True,
    cluster_cols=True
)
v3.data = data
v3.plot_heatmap(save_path="config3_correlation.png")

print("All configurations generated!")
```

---

## Troubleshooting

### Problem 1: Import Error

**Error Message**:
```
ModuleNotFoundError: No module named 'pandas'
```

**Solution**:
```bash
pip install pandas numpy matplotlib seaborn scipy
```

### Problem 2: CSV Encoding Error

**Error Message**:
```
UnicodeDecodeError: 'utf-8' codec can't decode byte...
```

**Solutions**:

Option 1: Convert using Excel
1. Open CSV file with Excel
2. Save As `.xlsx` format
3. Use the Excel file

Option 2: Modify CSV encoding
1. Open CSV with Notepad
2. When saving, select UTF-8 encoding

Option 3: The program automatically tries multiple encodings (utf-8, gbk, gb18030, etc.)

### Problem 3: Chinese Characters Display as Squares

**Solution**:

Add the following code at the beginning of your program:

```python
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False
```

### Problem 4: Inconsistent Column Count Error

**Error Message**:
```
ParserError: Error tokenizing data. C error: Expected 15 fields in line 4, saw 16
```

**Solutions**:

Option 1: Check with Excel
1. Open CSV with Excel
2. Check if each row has consistent column counts
3. Fix and save as `.xlsx`

Option 2: Fix data with Python

```python
import pandas as pd

# Read data and skip malformed lines
df = pd.read_csv("data.csv", index_col=0, on_bad_lines='skip')
df.to_csv("data_fixed.csv")
```

### Problem 5: Out of Memory

**Error Message**:
```
MemoryError
```

**Solutions**:

Option 1: Reduce data size
```python
# Only read partial data
df = pd.read_csv("data.csv", index_col=0, nrows=1000)
```

Option 2: Disable clustering
```python
visualizer = HeatmapVisualizer(
    cluster_rows=False,
    cluster_cols=False
)
```

### Problem 6: Image Too Large

**Solutions**:

Lower DPI:
```python
visualizer = HeatmapVisualizer(dpi=150)  # Reduce from 300 to 150
```

Reduce font sizes:
```python
visualizer = HeatmapVisualizer(
    title_size=12,
    label_size=10,
    tick_size=6
)
```

### Problem 7: Color Bar Position Incorrect

**Solutions**:

Adjust color bar ratios:
```python
visualizer = HeatmapVisualizer(
    cbar_width_ratio=0.15,   # 15% width (default)
    cbar_height_ratio=None   # Auto height (default)
)
```

---

## Advanced Tips

### Tip 1: Custom Colormaps

Use matplotlib's colormap:

```python
import matplotlib.pyplot as plt
from heatmap_visualizer import HeatmapVisualizer

# Create custom colormap
colors = ["#2E86AB", "#A23B72", "#F18F01", "#C73E1D"]
custom_cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
    "custom", colors
)

visualizer = HeatmapVisualizer(
    cmap=custom_cmap
)
```

### Tip 2: Mask Specific Values

```python
import numpy as np
from heatmap_visualizer import HeatmapVisualizer

visualizer = HeatmapVisualizer(data_file="data.csv")
visualizer.load_data()

# Mask NaN values or specific ranges
visualizer.data = visualizer.data.mask(
    (visualizer.data < -10) | (visualizer.data > 10)
)

visualizer.plot_heatmap(save_path="masked.png")
```

### Tip 3: Add Annotations

```python
import seaborn as sns

visualizer = HeatmapVisualizer(data_file="data.csv")
visualizer.load_data()

# Plot with seaborn (with annotations)
sns.heatmap(
    visualizer.data,
    annot=True,        # Show values
    fmt=".2f",        # 2 decimal places
    cmap="viridis"
)
```

---

## Appendix

### A. Complete Parameter List

```python
HeatmapVisualizer(
    data_file=None,           # Data file path
    title="Heatmap Analysis", # Chart title
    xlabel="Columns",         # X-axis label
    ylabel="Rows",            # Y-axis label
    cmap="viridis",           # Color scheme
    dpi=300,                  # Output resolution
    cluster_rows=True,        # Row clustering
    cluster_cols=True,        # Column clustering
    vmin=None,                # Color minimum
    vmax=None,                # Color maximum
    title_size=16,            # Title font size
    label_size=12,            # Axis label font size
    tick_size=8,              # Tick font size
    cbar_width_ratio=0.15,    # Color bar width ratio
    cbar_height_ratio=None    # Color bar height ratio
)
```

### B. Supported File Formats

| Format | Extension | Encoding Support |
|--------|-----------|------------------|
| CSV | .csv | UTF-8, GBK, GB18030, GB2312, UTF-16 |
| Excel | .xlsx | All Excel-supported encodings |
| Excel | .xls | All Excel-supported encodings |

### C. Complete Colormap List

**Sequential**:
viridis, plasma, inferno, magma, cividis, rocket, mako

**Diverging**:
coolwarm, RdBu, RdYlBu, PiYG, PRGn, Spectral, seismic, bwr

**Monochrome**:
Greens, Blues, Reds, Greys, Purples, Oranges

### D. Font Size Reference Values

| Purpose | title_size | label_size | tick_size |
|---------|-----------|-----------|-----------|
| Publication | 16-18 | 10-12 | 7-9 |
| Presentation | 20-24 | 14-18 | 10-12 |
| Poster | 24-30 | 18-24 | 12-16 |
| Preview | 12-14 | 8-10 | 6-8 |

### E. Color Bar Adjustment Reference Values

| Purpose | cbar_width_ratio | cbar_height_ratio |
|---------|-----------------|-------------------|
| Narrow bar | 0.05-0.10 | 0.5-0.7 |
| Medium bar | 0.12-0.18 | 0.7-0.9 |
| Wide bar | 0.20-0.30 | 0.9-1.0 |

---

## Changelog

### v1.0.0 (2024)
- Initial release
- Support multiple data formats
- Cluster analysis features
- Rich color schemes
- Highly customizable
