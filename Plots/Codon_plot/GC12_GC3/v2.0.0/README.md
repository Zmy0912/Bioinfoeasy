# GC12-GC3 Plotter

A Python-based GUI application for plotting GC12-GC3 neutrality graphs from `.out` files, designed for codon usage bias analysis in molecular biology research. It is part of the same series as the `GC3-ENC Plotter` and `PR2-Plot Plotter`, sharing the same layout logic and interface style.

## Features

- **Data Visualization**: Plot GC12 (%) (Y axis) against GC3 (%) (X axis) as a scatter of observed genes
- **Neutral Reference**: Display the neutral expectation line `y = x`, the baseline for neutrality analysis
- **Linear Regression**: Compute and draw the linear regression line `GC12 = slope × GC3 + intercept`
- **Regression Annotation**: Show the regression equation (`y = slope·x + intercept`) and `R²` in the top-left corner
- **Symmetric Plot**: Both axes share the identical `0–100 (%)` range, and the canvas keeps an equal square aspect ratio
- **Species Annotation**: Add species names as the graph title
- **Multiple Output Formats**: Support for PDF, SVG, PNG, and EPS formats
- **Fixed Canvas Size**: Maintains square aspect ratio regardless of window size
- **Gene Function Coloring**: Optionally color every scatter point by the functional category of its gene (e.g. photosystem, ATP synthase, NADH dehydrogenase, ribosomal proteins, RNA polymerase, ...), with a per-category legend, customizable category colors and per-category visibility
- **Comprehensive Style Customization**: Full control over all graph elements including:
  - Title, axis labels, and tick labels (size and color)
  - Neutral line (line width, style, color, visibility)
  - Regression line (line width, style, color, visibility)
  - Scatter points (size, color, transparency, edge color, visibility)
  - Gene function categories (color and visibility of each category)
  - Regression statistics annotation (font size, color, visibility)
  - Legend (font size, color, visibility)

## Requirements

- Python 3.6 or higher
- tkinter (usually included with Python)
- matplotlib
- numpy

## Installation

1. Ensure Python is installed on your system
2. Install required packages:

```bash
pip install matplotlib numpy
```

3. Download the `gc12_gc3_plotter.py` file

## Usage

1. Run the program:

```bash
python gc12_gc3_plotter.py
```

2. **Load Data File**:
   - Click "Browse..." to select your `.out` file
   - The file should contain codon usage columns including `GC3s` and overall `GC`

3. **Set Species Name**:
   - Enter the species name in the "Species name" field
   - This will appear as the graph title

4. **Plot Graph**:
   - Click "Plot Graph" to generate the visualization
   - The graph displays:
     - Red dashed line: the neutral expectation `y = x`
     - Blue dashed line: the linear regression line
     - Blue points: observed data from your input file
     - Text (top-left): regression equation and `R²`

5. **Color Points by Gene Function** (optional):
   - Tick the checkbox **"Color points by gene function"** in the top panel
   - Each point is colored according to the functional category of its gene and the legend lists every category with its counts
   - Untick the checkbox to go back to the single-color scatter

6. **Customize Style**:
   - Click "Style Settings" to open the customization panel
   - Use the tabs to adjust different elements:
     - **Title & Labels**: Adjust title, axis labels, and tick labels
     - **Neutral Line**: Modify the neutral expectation line appearance (width, style, color, visibility)
     - **Regression Line**: Modify the regression line appearance
     - **Scatter Points**: Customize data point appearance
     - **Annotation**: Control the regression equation / R² text
     - **Legend**: Control legend visibility and style
     - **Gene Function**: Toggle function-based coloring, set the color of every functional category and show/hide individual categories
   - Click "Apply" to save changes (graph will redraw automatically)
   - Click "Reset to Default" to restore original settings

7. **Save Image**:
   - Click "Save Image" to export the graph
   - Choose from PDF, SVG, PNG, or EPS formats
   - Saved images maintain the fixed square aspect ratio

## Input File Format

The input `.out` file should be a whitespace-separated text file with the following structure:

- First line: Header (automatically skipped)
- Subsequent lines: Data rows
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
psbB_[73610:75136](+)|01. 0.4839  0.1590  0.3780  0.1753  0.187   -0.057  0.383   48.64   0.270   0.320
```

Invalid data points (marked as `*****`) are automatically skipped. Both fractional (`0–1`) and percentage (`0–100`) values are accepted; values are internally converted to percentages.

## Gene Function Categories

When "Color points by gene function" is enabled, each gene is assigned to a category
using the leading identifier of its name (the part before the first `_`, so
`rpl20_[69725:70078](-)|01` is recognised as `rpl20`). The standard chloroplast
genome classification is used:

| Category | Genes |
|----------|-------|
| Photosystem I | `psaA`–`psaJ`, `ycf3`, `ycf4` |
| Photosystem II | `psbA`–`psbZ` |
| Cytochrome b/f complex | `petA`, `petB`, `petD`, `petG`, `petL`, `petN` |
| ATP synthase | `atpA`, `atpB`, `atpE`, `atpF`, `atpH`, `atpI` |
| NADH dehydrogenase | `ndhA`–`ndhK` |
| Rubisco (rbcL) | `rbcL` |
| Ribosomal proteins (SSU) | `rps*` |
| Ribosomal proteins (LSU) | `rpl*` |
| RNA polymerase | `rpoA`, `rpoB`, `rpoC1`, `rpoC2` |
| Maturase (matK) | `matK` |
| Protease (clpP) | `clpP` |
| Envelope membrane (cemA) | `cemA` |
| C-type cytochrome synthesis (ccsA) | `ccsA` |
| Translation initiation factor (infA) | `infA` |
| Acetyl-CoA carboxylase (accD) | `accD` |
| Conserved ORFs (ycf) | `ycf1`, `ycf2`, other `ycf*` / `orf*` |
| Transfer RNA | `trn*` |
| Ribosomal RNA | `rrn*` |
| Other / Unclassified | any gene matching none of the above |

## Theoretical Background

The GC12-GC3 (neutrality) plot is used to test the **mutation-selection equilibrium hypothesis** of codon usage bias. Under neutrality, the GC content at the three codon positions should be equal, so a plot of `GC12` against `GC3` is expected to follow the diagonal `y = x`.

For each gene, two quantities are used:

```
X axis = GC3 (%)                       -> GC content at the 3rd codon position
Y axis = GC12 (%) = (GC1 + GC2) / 2    -> GC content at the 1st and 2nd codon positions
```

Since the overall GC content is the average over the three positions:

```
GC = (GC1 + GC2 + GC3) / 3
=> GC1 + GC2 = 3 × GC − GC3
=> GC12 = (GC1 + GC2) / 2 = (3 × GC − GC3) / 2
```

Thus GC12 is derived from the `GC` and `GC3s` columns of the `.out` file.

- The line `y = x` represents the neutral expectation.
- The linear regression `GC12 = slope × GC3 + intercept` describes the actual trend.

## Graph Interpretation

- **Regression slope close to 1 and points near `y = x`**: Weak mutation pressure / near-neutrality
- **Shallow slope (far below 1)**: Strong selection constraint relative to mutation pressure
- **`R²` value**: Goodness of fit of the linear regression; a low `R²` indicates that GC3 explains little of the GC12 variation
- **Points far from `y = x`**: The three codon positions have unequal GC content, suggesting non-neutral evolution

## Style Settings Reference

### Title & Labels Tab
- **Title**: Graph title font size (8-24) and color
- **X-axis Label**: X-axis title font size (8-20) and color
- **Y-axis Label**: Y-axis title font size (8-20) and color
- **X-axis Tick Labels**: X-axis numbers font size (6-18) and color
- **Y-axis Tick Labels**: Y-axis numbers font size (6-18) and color

### Neutral Line Tab
- **Visibility**: Toggle the neutral expectation line (`y = x`)
- **Line Width**: 0.5 to 5
- **Line Style**: Solid (-), Dashed (--), Dotted (:), Dash-dot (-.)
- **Color**: Customizable

### Regression Line Tab
- **Visibility**: Toggle the linear regression line display
- **Line Width**: 0.5 to 5
- **Line Style**: Solid (-), Dashed (--), Dotted (:), Dash-dot (-.)
- **Color**: Customizable

### Scatter Points Tab
- **Visibility**: Toggle points display
- **Point Size**: 10 to 200
- **Transparency**: 0.1 to 1.0
- **Edge Width**: 0 to 3
- **Fill Color**: Customizable
- **Edge Color**: Customizable

### Annotation Tab
- **Visibility**: Toggle the regression equation / `R²` text
- **Font Size**: 8 to 20
- **Text Color**: Customizable

### Legend Tab
- **Visibility**: Toggle legend display
- **Font Size**: 6 to 16
- **Text Color**: Customizable

### Gene Function Tab
- **Color scatter points by gene function**: Enable / disable function-based coloring
- **Show (per category)**: Hide a single category from the graph without disabling the whole feature
- **Color (per category)**: Choose the color used for each functional category

## Troubleshooting

**Issue**: "Please select a valid input file!"
- **Solution**: Ensure you have selected a `.out` file with the correct format

**Issue**: "Could not extract valid data from file!"
- **Solution**: Verify your file has at least 11 columns and contains valid `GC3s` (column 9) and `GC` (column 10) data

**Issue**: Graph appears distorted
- **Solution**: The graph is designed to be 600×600 pixels. Adjust window size if needed.

**Issue**: Style changes not visible
- **Solution**: Click "Apply" in Style Settings. Changes only appear after clicking Apply.

## Technical Details

- **Framework**: Tkinter with matplotlib backend
- **Canvas Size**: 600×600 pixels (fixed)
- **Output Resolution**: 300 DPI for raster formats
- **Axis Range**: Both axes fixed to `0–100 (%)` (equal numeric span)
- **Aspect Ratio**: Equal (square) via `set_aspect('equal')`
- **Regression**: Ordinary least squares via `numpy.polyfit`

## License

This tool is provided as-is for research and educational purposes.

## Citation

If you use this tool in your research, please cite the theoretical foundation:

- Sueoka, N. (1988). Directional mutation pressure and neutral molecular evolution. *Proceedings of the National Academy of Sciences*, 85(8), 2653-2657.
- Sueoka, N. (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. *Journal of Molecular Evolution*, 40(3), 318-325.
- Lobry, J. R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution*, 13(5), 660-665.

## Contact

For questions or suggestions, please contact the developer.

## Version History

- **v1.0**: Initial release with GC12-GC3 neutrality plotting, linear regression, neutral reference line and style customization
- **v1.1**: Added gene function classification coloring (per-category colors, per-category visibility, category legend)
