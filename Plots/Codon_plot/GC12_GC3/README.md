# GC12-GC3 Plotter

A GUI application for plotting GC12-GC3 neutrality graphs from codon usage
tables (e.g. CodonW `*.out` files), used to test the mutation-selection
equilibrium hypothesis of codon usage bias.

## Contents

| Version | Description |
|---------|-------------|
| [`v1.0.0/`](v1.0.0/) | Initial release: neutrality plot, linear regression, neutral reference line, style customization |
| [`v2.0.0/`](v2.0.0/) | Adds gene-function coloring with per-category colors and visibility |
| [`v3.0.0/`](v3.0.0/) | Current release — adds folder (batch) mode with combined grid figures (per-panel regression stats, shared side legend, per-file export) and outlier gene labels |

Each version contains `gc12_gc3_plotter.py` and a `README.md` with full usage
instructions and the interpretation of the plot.

## Companion tool

| Tool | Description |
|------|-------------|
| [`gc12_gc3_regression_table/`](gc12_gc3_regression_table/) | Uses exactly the same calculation as `v3.0.0`, but reports the statistics as a sortable table (one row per species) instead of drawing figures. Supports single-file and folder batch modes, a name prefix, and UTF-8 BOM CSV export. |

## Requirements

- Python 3.6+
- `matplotlib`, `numpy`
- `tkinter` (ships with CPython)

## Usage

```bash
# Plot version
cd v3.0.0
python gc12_gc3_plotter.py

# Table version (batch regression statistics)
cd gc12_gc3_regression_table
python gc12_gc3_regression_table.py
```

Back to [Codon Plot](../README.md) · [project root](../../../README.md).
