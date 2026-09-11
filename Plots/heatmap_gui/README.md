# Heatmap GUI

A graphical interface for creating clustered heatmaps from tabular data.

## Contents

| File | Description |
|------|-------------|
| `heatmap_gui_v3.py` | GUI entry point |
| `heatmap_visualizer.py` | Heatmap plotting backend |
| `start.bat` | Windows launcher |
| `requirements.txt` | Dependencies |
| `examples/data_sample.csv` | Sample input data |

## Requirements

- Python 3.7+

```bash
pip install -r requirements.txt
```

- `pandas>=1.5.0`
- `numpy>=1.23.0`
- `matplotlib>=3.6.0`
- `seaborn>=0.12.0`
- `scipy>=1.9.0`

## Usage

```bash
python heatmap_gui_v3.py
```

or double-click `start.bat` on Windows.

> Looking for the command-line version? See [`../Heatmap/`](../Heatmap/).

Back to [Plots](../README.md) · [project root](../../README.md).
