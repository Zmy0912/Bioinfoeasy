# ENC-GC3s Neutral Expectation Curve Analyzer

A GUI tool that quantifies the proportion of genes falling **below**, **on**, or **above** Wright's (1990) neutral expectation curve in an ENC-GC3s plot, for a whole folder of codonw `.out` files (one file per species).

This tool was written to address the following reviewer comment:

> 3- ENC-GC3s analysis does not specify how expected ENC values under neutrality were calculated; authors must provide Wright's formula for expected ENC calculation in Methods, quantify the proportion of genes falling below/on/above the expected curve, and use this to strengthen their conclusion about selection dominance.

## Requirements

- Python 3.7+
- matplotlib

```bash
pip install matplotlib
```

## Usage

```bash
python ENC_GC3s_Analyzer.py
```

1. Click **Browse...** and select the folder containing your `.out` files (codonw output, tab-delimited, must contain `Nc` and `GC3s` columns). Each `.out` file is treated as one species.
2. Optionally adjust the tolerance (see below).
3. Click **Run analysis**.
4. View results in the four tabs:
   - **Species summary** — per-species gene counts and percentages below/on/above the curve, with a TOTAL row.
   - **Gene detail** — per-gene GC3s, observed Nc, expected ENC, difference, and class.
   - **ENC-GC3s plot** — scatter plot with the expected curve (select a species row in the summary tab to plot that species alone).
   - **Methods** — ready-to-use Methods text and conclusion wording for the manuscript.
5. Export the summary table and/or gene detail as CSV (UTF-8 with BOM, opens directly in Excel).

## Method

The expected ENC under neutrality (no selection) is calculated with Wright's (1990) formula:

```
ENC_exp = 2 + GC3s + 29 / (GC3s^2 + (1 - GC3s)^2)
```

where `GC3s` is the GC content at synonymous third codon positions (codonw output).

Each gene is classified by comparing its observed Nc with its expected ENC:

| Class | Criterion | Interpretation |
|---|---|---|
| Below | `Nc_obs < ENC_exp - tol` | Codon usage bias stronger than the neutral expectation — indicates natural selection |
| On | `\|Nc_obs - ENC_exp\| <= tol` | Consistent with the neutral (mutation-pressure) expectation |
| Above | `Nc_obs > ENC_exp + tol` | Bias weaker than the neutral expectation — consistent with mutation pressure dominance |

**Tolerance modes** (configurable in the GUI):

- `relative` (default, 0.05): tolerance is 5% of the expected ENC value.
- `absolute`: tolerance is a fixed ENC difference (e.g., 1.0 ENC units).

Genes with invalid codonw values (e.g., `*****` in the Nc column) are skipped automatically.

## Interpreting the results for the selection-dominance conclusion

- If the **majority** of genes lie **on or above** the expected curve and only a minority fall below it, codon usage is primarily shaped by **mutation pressure**, with natural selection acting on a subset of genes.
- If the **majority** of genes fall **below** the curve, **natural selection** dominates the codon usage bias.

## Output columns

Summary CSV:

| Column | Description |
|---|---|
| Species | File name (without extension) of the `.out` file |
| Total genes | Number of genes with valid Nc and GC3s values |
| Below (n) / Below (%) | Genes below the expected curve: count / percentage |
| On (n) / On (%) | Genes on the curve: count / percentage |
| Above (n) / Above (%) | Genes above the curve: count / percentage |

Gene detail CSV: `Species, Gene, GC3s, Observed_Nc, Expected_ENC, Nc_minus_expected, Class`

## Reference

Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23-29. https://doi.org/10.1016/0378-1119(90)90488-K
