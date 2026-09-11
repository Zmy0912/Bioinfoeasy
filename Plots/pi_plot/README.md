# Pi Plot

Tools for computing and visualizing nucleotide polymorphism (Pi) along a
sequence, with gene annotation.

## Contents

| Module | Description |
|--------|-------------|
| [`pi_gene_analyzer/`](pi_gene_analyzer/) | Gene nucleotide polymorphism (Pi) analyzer |
| [`Plot/`](Plot/) | Pi sliding-window gene annotation program, including a GUI (v2.0.0) |

### `Plot/` versions

| Version | Description |
|---------|-------------|
| [`Plot/v1.0.0/`](Plot/v1.0.0/) | Command-line pi plot |
| [`Plot/v1.1.0/`](Plot/v1.1.0/) | Command-line pi plot (updated) |
| [`Plot/v2.0.0/`](Plot/v2.0.0/) | Command-line (`pi_plot/`) and GUI (`pi_plot_gui/`, CN/EN) versions |

Each version includes a `README.md`, `USAGE_GUIDE.md` and `CHANGELOG.md` where
applicable.

## Requirements

- Python 3.6+
- `pandas`, `numpy`, `matplotlib`

## Usage

```bash
# Compute Pi per gene
cd pi_gene_analyzer
python pi_gene_analyzer.py

# Plot (GUI, English)
cd ../Plot/v2.0.0/pi_plot_gui/EN
python pi_plot_gui.py
```

Back to [Plots](../README.md) · [project root](../../README.md).
