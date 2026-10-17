# GC12-GC3 Regression Table

A Python GUI application that computes the GC12-GC3 (neutrality plot) regression
statistics for one or many species and presents them as a **table**. It uses exactly
the same calculation method as the `GC12-GC3 Plotter` (v3.0.0), but instead of
drawing graphs it summarizes every species in one row: linear equation, slope,
intercept, r, R², and more.

## Features

- **Two Input Modes**:
  - **Single file**: analyze one `.out` file
  - **Folder (batch)**: analyze every `.out` file in a folder at once
- **Table Output**: one row per species (file), with the following columns:
  - Species / File (the file name without extension, optionally prefixed)
  - Linear equation (`y = slope·x + intercept`)
  - Slope, Intercept
  - r (Pearson correlation coefficient), R² (coefficient of determination)
  - Number of genes (n)
  - GC3 mean (%), GC12 mean (%)
  - Standard error of the slope
- **Sortable Table**: click any column header to sort ascending; click again to sort descending
- **Name Prefix**: prepend a common prefix (e.g. a genus name) to every result name
- **CSV Export**: save the whole table as a UTF-8 (BOM) CSV that opens cleanly in Excel
- **Robust Parsing**: invalid values (`*****`) are skipped automatically; both fractional (`0–1`) and percentage (`0–100`) values are accepted

## Requirements

- Python 3.6 or higher
- tkinter (usually included with Python)
- numpy

```bash
pip install numpy
```

## Usage

1. Run the program:

```bash
python gc12_gc3_regression_table.py
```

2. **Choose the Input Mode** (top of the window):
   - **Single file**: click "Browse File..." and select one `.out` file
   - **Folder**: click "Browse Folder..." and select a folder that directly contains
     `.out` files (sub-folders are not scanned). The default path points to the
     `out` folder next to this script, if it exists.

3. **(Optional) Name prefix**: enter a prefix that is prepended to every file
   name, e.g. `Asparagus_` turns `NC01` into `Asparagus_NC01`.

4. Click **Analyze**. The table shows one row per species / file.

5. Click any **column header** to sort the table by that column
   (click again to reverse the order).

6. Click **Export CSV** to save the table (choose the destination and file name).
   The CSV is written with a UTF-8 BOM so that Excel displays it correctly.

## Input File Format

The input `.out` file is a whitespace-separated text file (codon usage table):

- First line: header (automatically skipped)
- Required columns (0-indexed):

| Index | Column | Meaning |
|-------|--------|---------|
| 0     | title  | Gene name / identifier |
| 9     | GC3s   | GC content at the 3rd codon position |
| 10    | GC     | Overall GC content of the gene |

Example:

```
title                      T3s     C3s     A3s     G3s     CAI     CBI     Fop     Nc      GC3s    GC
rbcL_[56338:57780](+)|01.  0.5090  0.1860  0.4169  0.1295  0.264   0.059   0.459   47.80   0.250   0.313
psbB_[73610:75136](+)|01.  0.4839  0.1590  0.3780  0.1753  0.187  -0.057   0.383   48.64   0.270   0.320
```

## Theoretical Background

The GC12-GC3 (neutrality) plot tests the mutation-selection equilibrium
hypothesis of codon usage bias. For every gene two quantities are used:

```
X axis = GC3 (%)                     -> GC content at the 3rd codon position
Y axis = GC12 (%) = (GC1 + GC2) / 2  -> mean GC content of the 1st and 2nd positions
```

Because the overall GC content is the average over the three positions:

```
GC = (GC1 + GC2 + GC3) / 3
=> GC1 + GC2 = 3 × GC − GC3
=> GC12 = (GC1 + GC2) / 2 = (3 × GC − GC3) / 2
```

so GC12 is derived from the `GC` and `GC3s` columns of the `.out` file.

- The line `y = x` represents the neutral expectation (slope = 1).
- The linear regression `GC12 = slope × GC3 + intercept` describes the actual trend.

## Table Interpretation

- **Slope close to 1**: weak mutation pressure / near-neutrality
- **Shallow slope (far below 1)**: strong selection constraint relative to mutation pressure
- **R²**: goodness of fit of the linear regression; a low R² indicates that GC3
  explains little of the GC12 variation
- **Slope SE**: standard error of the slope, useful to judge whether the slope
  differs significantly from 1
- **Genes (n)**: number of valid genes that entered the regression of that species

## Technical Details

- **Framework**: Tkinter (Treeview table)
- **Regression**: ordinary least squares via `numpy.polyfit`
- **Parsing**: identical to the `GC12-GC3 Plotter` (v3.0.0), so both tools give
  exactly the same slope / intercept / R² for the same input file
- **CSV encoding**: UTF-8 with BOM (`utf-8-sig`)

## Troubleshooting

**Issue**: "Please select a valid input file!"
- Ensure you have selected an existing `.out` file (single-file mode).

**Issue**: "Please select a valid input folder!" / "No .out files found in the selected folder!"
- Switch to **Folder** mode and select a directory that directly contains `.out`
  files (sub-folders are not scanned).

**Issue**: "No valid data could be extracted from any file!"
- Verify the files have at least 11 columns with valid `GC3s` (column 9) and
  `GC` (column 10) data.

**Issue**: Excel shows garbled characters when opening the CSV
- The CSV is written as UTF-8 with BOM; use "Data → From Text/CSV" and pick
  UTF-8 if your Excel version ignores the BOM.

## License

This tool is provided as-is for research and educational purposes.

## Citation

If you use this tool in your research, please cite the theoretical foundation:

- Sueoka, N. (1988). Directional mutation pressure and neutral molecular evolution. *Proceedings of the National Academy of Sciences*, 85(8), 2653-2657.
- Sueoka, N. (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. *Journal of Molecular Evolution*, 40(3), 318-325.
- Lobry, J. R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution*, 13(5), 660-665.

## Version History

- **v1.0**: Initial release — single-file / folder batch analysis, per-species
  regression statistics in a sortable table, optional name prefix, CSV export
