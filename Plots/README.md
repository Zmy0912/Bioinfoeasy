# Plots

Visualization tools for sequence, codon usage, and phylogenetic data.

## Contents

| Module | Description |
|--------|-------------|
| [`Codon_plot/`](Codon_plot/) | Codon usage and GC analysis plots (codon usage, GC12-GC3, GC3-ENC, PR2-Plot) |
| [`Heatmap/`](Heatmap/) | Heatmap data visualization (command line) |
| [`heatmap_gui/`](heatmap_gui/) | Heatmap data visualization (graphical interface) |
| [`pi_plot/`](pi_plot/) | Pi (nucleotide polymorphism) gene analysis and sliding-window plots |
| [`time_scale_generator/`](time_scale_generator/) | Phylogenetic time scale generator (CN/EN) |

## Dependency notes

Each tool lists its own dependencies in a `requirements.txt`. Several tools
rely on `tkinter` (bundled with CPython) or `PyQt5` for their GUI.

Back to the [project root](../README.md).
