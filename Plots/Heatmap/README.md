# Heatmap Data Visualization Tool

A professional heatmap data visualization tool supporting CSV/Excel data loading, cluster analysis, and multiple color schemes. Suitable for gene expression analysis, correlation analysis, data pattern recognition, and various other scenarios.

## Features

- **Multiple Data Formats**: CSV, Excel (.xlsx, .xls)
- **Automatic Encoding Detection**: Supports UTF-8, GBK, GB18030 and more
- **Cluster Analysis**: Supports row clustering, column clustering, or bidirectional clustering
- **Rich Color Schemes**: 20+ professional color schemes
- **Highly Customizable**: Adjust font sizes, color bar dimensions, output resolution, etc.
- **High Quality Output**: Supports high DPI image export (suitable for academic publication)

## Installation

```bash
pip install pandas numpy matplotlib seaborn scipy
```

Or install using requirements.txt:

```bash
pip install -r requirements.txt
```

## Quick Start

### Option 1: Using Sample Data

```python
from heatmap_visualizer import HeatmapVisualizer

# Create visualizer and generate sample data
visualizer = HeatmapVisualizer(
    title="Gene Expression Heatmap Analysis",
    xlabel="Genes",
    ylabel="Samples",
    cmap="viridis"
)

# Generate sample data
visualizer.generate_sample_data(n_rows=20, n_cols=15)

# Plot and save heatmap
visualizer.plot_heatmap(save_path="heatmap.png")
visualizer.show()
```

### Option 2: Using Your Own Data

```python
from heatmap_visualizer import HeatmapVisualizer

# Create visualizer
visualizer = HeatmapVisualizer(
    data_file="your_data.csv",  # or "your_data.xlsx"
    title="Custom Title",
    xlabel="Column Label",
    ylabel="Row Label",
    cmap="coolwarm"
)

# Load data and plot
visualizer.load_data()
visualizer.plot_heatmap(save_path="my_heatmap.png")
visualizer.show()
```

## Data Format Requirements

### CSV File Format

CSV files should meet the following requirements:
- First column: row index (sample names)
- First row: column names (gene names, etc.)
- All cells: numeric data
- Consistent number of columns per row

Example CSV file content:

```csv
,Gene1,Gene2,Gene3,Gene4
Sample1,1.2,2.3,3.4,4.5
Sample2,2.1,3.2,4.3,5.4
Sample3,3.0,4.1,5.2,6.3
```

### Excel File Format

Excel files have the same requirements as CSV. First column is row index, first row is column names.

## Parameter Reference

### Core Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_file` | str | None | Data file path (CSV or Excel) |
| `title` | str | "Heatmap Analysis" | Chart title |
| `xlabel` | str | "Columns" | X-axis label |
| `ylabel` | str | "Rows" | Y-axis label |
| `cmap` | str | "viridis" | Color scheme |
| `dpi` | int | 300 | Output image resolution |

### Clustering Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cluster_rows` | bool | True | Enable row clustering |
| `cluster_cols` | bool | True | Enable column clustering |

### Color Range Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `vmin` | float | None | Color scale minimum (None=auto) |
| `vmax` | float | None | Color scale maximum (None=auto) |

### Font Size Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `title_size` | int | 16 | Title font size (recommended 14-24) |
| `label_size` | int | 12 | Axis label font size (recommended 10-16) |
| `tick_size` | int | 8 | Tick label font size (recommended 6-12) |

### Color Bar Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cbar_width_ratio` | float | 0.15 | Color bar width ratio (recommended 0.05-0.3) |
| `cbar_height_ratio` | float | None | Color bar height ratio (recommended 0.3-1.0, None=auto) |

## Color Schemes

### Sequential Colormaps (Recommended for General Data)

| Scheme | Colors | Use Case |
|--------|--------|----------|
| `viridis` | Purple-Blue-Green-Yellow | Default, colorblind-friendly |
| `plasma` | Blue-Red-Yellow | High contrast |
| `inferno` | Black-Red-Yellow | High contrast, dark background |
| `magma` | Black-Purple-Yellow | High contrast |

### Diverging Colormaps (Recommended for Correlation/Diff Analysis)

| Scheme | Colors | Use Case |
|--------|--------|----------|
| `coolwarm` | Blue-White-Red | Correlation analysis, difference analysis |
| `RdBu` | Red-White-Blue | Positive-negative comparison |
| `RdYlBu` | Red-Yellow-Blue | Multi-level classification |
| `PiYG` | Purple-White-Green | Balanced comparison |
| `PRGn` | Purple-White-Green | Balanced comparison |

### Single Color Colormaps

| Scheme | Colors | Use Case |
|--------|--------|----------|
| `Greens` | Light-Dark Green | Positive values |
| `Blues` | Light-Dark Blue | Positive values |
| `Reds` | Light-Dark Red | Positive values |

## Usage Examples

### Example 1: Gene Expression Heatmap (viridis)

```python
visualizer = HeatmapVisualizer(
    data_file="gene_expression.csv",
    title="Gene Expression Heatmap Analysis",
    xlabel="Genes",
    ylabel="Samples",
    cmap="viridis",
    cluster_rows=True,
    cluster_cols=True
)
visualizer.load_data()
visualizer.plot_heatmap(save_path="gene_heatmap.png")
visualizer.show()
```

### Example 2: Correlation Heatmap (coolwarm)

```python
visualizer = HeatmapVisualizer(
    data_file="correlation.csv",
    title="Correlation Analysis Heatmap",
    xlabel="Variables",
    ylabel="Variables",
    cmap="coolwarm",
    vmin=-1,
    vmax=1,
    cluster_rows=True,
    cluster_cols=True
)
visualizer.load_data()
visualizer.plot_heatmap(save_path="correlation_heatmap.png")
visualizer.show()
```

### Example 3: Large Font + Wide Color Bar

```python
visualizer = HeatmapVisualizer(
    data_file="data.csv",
    title="Large Title Heatmap",
    xlabel="Columns",
    ylabel="Rows",
    cmap="plasma",
    title_size=20,         # Large title
    label_size=14,         # Large axis labels
    tick_size=10,          # Large ticks
    cbar_width_ratio=0.25, # Wide color bar
    cbar_height_ratio=0.9, # Tall color bar
    dpi=300
)
visualizer.load_data()
visualizer.plot_heatmap(save_path="large_heatmap.png")
visualizer.show()
```

### Example 4: No Clustering

```python
visualizer = HeatmapVisualizer(
    data_file="data.csv",
    title="No Clustering Heatmap",
    xlabel="Columns",
    ylabel="Rows",
    cmap="viridis",
    cluster_rows=False,
    cluster_cols=False
)
visualizer.load_data()
visualizer.plot_heatmap(save_path="no_cluster_heatmap.png")
visualizer.show()
```

### Example 5: Using Sample Data for Quick Testing

```python
from heatmap_visualizer import HeatmapVisualizer

visualizer = HeatmapVisualizer(
    title="Sample Data Heatmap",
    xlabel="Genes",
    ylabel="Samples",
    cmap="coolwarm"
)

# Generate 20 rows × 15 columns of sample data
visualizer.generate_sample_data(n_rows=20, n_cols=15)
visualizer.plot_heatmap(save_path="sample_heatmap.png")
visualizer.show()
```

## Running Example Program

The project includes a complete example program. Run it directly to generate multiple example heatmaps:

```bash
python heatmap_visualizer.py
```

This program will generate the following examples:
- `heatmap_example_viridis.png` - Default viridis colormap
- `heatmap_example_coolwarm.png` - coolwarm colormap + large font + wide bar
- `heatmap_example_plasma.png` - plasma colormap + small font + narrow bar

## FAQ

### Q1: CSV File Reading Error

**Problem**: `ParserError` or encoding error

**Solutions**:
- Ensure file uses UTF-8, GBK, or GB18030 encoding
- Open with Excel and save as `.xlsx` format
- Check that all rows have consistent column counts

### Q2: Chinese Characters Display as Squares

**Problem**: Chinese characters not displaying properly

**Solutions**:
- Ensure Chinese fonts are installed (SimHei, Microsoft YaHei)
- Manually set fonts:
  ```python
  plt.rcParams['font.sans-serif'] = ['SimHei']
  plt.rcParams['axes.unicode_minus'] = False
  ```

### Q3: Image Size Inappropriate

**Solutions**:
- Image size adjusts automatically based on data dimensions
- Manually adjust font sizes to optimize layout
- Output DPI recommended to be 300 (for publication)

### Q4: Color Bar Too Large or Too Small

**Solutions**:
- Adjust `cbar_width_ratio` (width, 0.05-0.3)
- Adjust `cbar_height_ratio` (height, 0.3-1.0 or None)

## API Reference

### HeatmapVisualizer Class

#### Methods

- `__init__(...)` - Initialize visualizer
- `load_data(data_file=None)` - Load data file
- `generate_sample_data(n_rows=15, n_cols=12, seed=42)` - Generate sample data
- `plot_heatmap(save_path=None)` - Plot heatmap
- `save_heatmap(save_path)` - Save heatmap
- `show()` - Display heatmap

## Project Structure

```
Heatmap/
├── heatmap_visualizer.py   # Main program file
├── use_heatmap.py          # Usage example
├── README.md               # This file
├── USAGE.md                # Detailed usage guide
├── requirements.txt        # Dependency list
└── examples/               # Example folder
    └── data_sample.csv     # Sample data file
```

## License

MIT License

## Changelog

### v1.0.0
- Initial release
- Support CSV/Excel data loading
- Support row/column clustering
- Support 20+ color schemes
- Support custom font sizes and color bar dimensions
- Automatic encoding detection
