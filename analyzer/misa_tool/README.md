# MISA Tool

Toolkit for working with the output of
[MISA](https://webblast.ipk-gatersleben.de/misa/) (MIcroSAtellite) searches.

## Contents

| Module | Description |
|--------|-------------|
| [`misa_parser/`](misa_parser/) | Parses a `.misa` result file into tidy tables (per-species motif matrix, per-category counts, statistics). |
| [`misa_ssr_viewer/`](misa_ssr_viewer/) | Turns a `.misa` file into publication-style figures — motif heatmap, SSR-category pie chart and per-species stacked bar — and exports the underlying data tables. |

## Typical workflow

1. Run MISA on your genome/sequence set to obtain a `.misa` file.
2. Use `misa_parser` to turn it into tabular data (or feed the `.misa` file
   straight into `misa_ssr_viewer`).
3. Optionally classify each SSR by genomic region with
   [`SSR_Region_Classifier`](../SSR_Region_Classifier/), or count repeats with
   [`RepeatAnalyzer`](../RepeatAnalyzer/).

## Requirements

Each sub-module lists its own dependencies in a `requirements.txt`:

```bash
cd misa_parser      # pandas, openpyxl
pip install -r requirements.txt
python misa_parser.py
```

```bash
cd misa_ssr_viewer  # matplotlib, numpy (tkinter ships with CPython)
python misa_ssr_viewer.py
```

Back to [Analyzer](../README.md) · [project root](../../README.md).
