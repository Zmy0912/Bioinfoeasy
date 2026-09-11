# PR2-Plot Plotter

A GUI application for generating PR2 (Parity Rule 2) plots, used to examine
the AT/GC and purine/pyrimidine asymmetry of codon usage.

## Contents

| Version | Description |
|---------|-------------|
| [`v1.0.0/`](v1.0.0/) | Initial release |
| [`v2.0.0/`](v2.0.0/) | Intermediate release |
| [`v3.0.0/`](v3.0.0/) | GC12-GC3 plotter variant (kept for historical continuity) |
| [`v4.0.0/`](v4.0.0/) | Intermediate release |
| [`v5.0.0/`](v5.0.0/) | Current release |

Each version contains its plotter script and a `README.md` with the full usage
instructions.

> **Note:** `v3.0.0/` in this directory contains a `gc12_gc3_plotter.py`
> (a GC12-GC3 plotter) rather than a PR2 plotter, reflecting the historical
> evolution of this tool family. Use `v5.0.0` for PR2 plots.

## Requirements

- Python 3.6+
- `matplotlib`, `numpy`
- `tkinter` (ships with CPython)

## Usage

```bash
cd v5.0.0
python pr2_plot_plotter.py
```

Back to [Codon Plot](../README.md) · [project root](../../../README.md).
