# Codon Plot

Plotting tools for codon usage bias analysis. Each tool is a standalone GUI
application that reads codon usage tables (e.g. CodonW `*.out` files) and
produces publication-quality figures.

## Contents

| Tool | Description | Versions |
|------|-------------|----------|
| [`draw_codon_usage/`](draw_codon_usage/) | Codon usage bias plots (CN/EN) | `v0.1.0`, `v0.1.1` |
| [`GC12_GC3/`](GC12_GC3/) | GC12 vs GC3 neutrality plots, plus a regression statistics table | `v1.0.0` – `v3.0.0` |
| [`GC3_ENC/`](GC3_ENC/) | GC3 vs ENC (effective number of codons) plots | `v1.0.0`, `v2.0.0`, `v3.0.0` |
| [`PR2_plot/`](PR2_plot/) | PR2 (parity rule 2) plots | `v1.0.0` – `v5.0.0` |
| [`coa_plotter/`](coa_plotter/) | CodonW correspondence analysis (COA) supplementary figures, with correlation tests against GC3s / GC / GC12 | — |
| [`coa_simple_plotter/`](coa_simple_plotter/) | Streamlined COA variant — main orthogonal-axis panel only, batch input retained | — |

> Older versions superseded by the tools above are kept under
> [`archive/codon_plot_legacy/`](../../archive/codon_plot_legacy/).

## Related analysis tools

The corresponding statistics modules live under
[`analyzer/Chloroplast_genome/Codon/`](../../analyzer/Chloroplast_genome/Codon/).

Back to [Plots](../README.md) · [project root](../../README.md).
