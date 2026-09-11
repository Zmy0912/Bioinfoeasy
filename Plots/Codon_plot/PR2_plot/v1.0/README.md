# PR2-Plot Plotter

A Python-based GUI application for plotting PR2-bias (Parity Rule 2) graphs from `.out` files, designed for codon usage bias analysis in molecular biology research. It is part of the same series as the `GC3-ENC Plotter`, sharing the same layout logic and interface style.

## Features

- **Data Visualization**: Plot PR2-bias scatter of `G3/(G3+C3)` (X axis) against `A3/(A3+T3)` (Y axis)
- **Neutral Reference**: Display the reference crosshair at the unbiased center `(0.5, 0.5)`
- **Symmetric Plot**: X and Y axes share one identical, center-symmetric range, and the canvas keeps an equal square aspect ratio
- **Species Annotation**: Add species names as the graph title
- **Multiple Output Formats**: Support for PDF, SVG, PNG, and EPS formats
- **Fixed Canvas Size**: Maintains square aspect ratio regardless of window size
- **Gene Labels**: Automatically label the genes that deviate most from the center
- **Comprehensive Style Customization**: Full control over all graph elements including:
  - Title, axis labels, and tick labels (size and color)
  - Reference lines (line width, style, color, visibility)
  - Scatter points (size, color, transparency, edge color, visibility)
  - Legend (font size, color, visibility)
  - Gene labels (number of genes, font size, color, visibility)

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

3. Download the `pr2_plot_plotter.py` file

## Usage

1. Run the program:

```bash
python pr2_plot_plotter.py
```

2. **Load Data File**:
   - Click "Browse..." to select your `.out` file
   - The file should contain codon usage columns including T3s, C3s, A3s and G3s

3. **Set Species Name**:
   - Enter the species name in the "Species name" field
   - This will appear as the graph title

4. **Plot Graph**:
   - Click "Plot Graph" to generate the visualization
   - The graph displays:
     - Red dashed crosshair: the neutral reference at `(0.5, 0.5)`
     - Blue points: observed data from your input file

5. **Customize Style**:
   - Click "Style Settings" to open the customization panel
   - Use the tabs to adjust different elements:
     - **Title & Labels**: Adjust title, axis labels, and tick labels
     - **Reference Lines**: Modify the neutral crosshair appearance (width, style, color, visibility)
     - **Scatter Points**: Customize data point appearance
     - **Legend**: Control legend visibility and style
     - **Gene Labels**: Control the labels of the most biased genes
   - Click "Apply" to save changes (graph will redraw automatically)
   - Click "Reset to Default" to restore original settings

6. **Save Image**:
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
| 1     | T3s    | T at 3rd codon position |
| 2     | C3s    | C at 3rd codon position |
| 3     | A3s    | A at 3rd codon position |
| 4     | G3s    | G at 3rd codon position |

Example:

```
title                      T3s     C3s     A3s     G3s     CAI     CBI     Fop     Nc      GC3s    ...
rbcL_[56338:57780](+)|01.  0.5090  0.1860  0.4169  0.1295  0.264   0.059   0.459   47.80   0.250   ...
psbB_[73610:75136](+)|01. 0.4839  0.1590  0.3780  0.1753  0.187   -0.057  0.383   48.64   0.270   ...
```

Invalid data points (marked as `*****`) are automatically skipped.

## Theoretical Background

The PR2-bias plot is based on **Parity Rule 2** (Chargaff's second parity rule), which states that, within a single DNA strand, the content of A approximately equals T and the content of G approximately equals C. Deviations from this rule reveal strand-specific mutation and selection pressures.

For each gene, the third codon position is used to compute:

```
X axis = G3 / (G3 + C3)      -> GC bias at the 3rd codon position
Y axis = A3 / (A3 + T3)      -> AT bias at the 3rd codon position
```

- The center `(0.5, 0.5)` represents no bias (perfect parity).
- A point far from the center indicates a strong bias between the two strands.

## Graph Interpretation

- **Points at the center (0.5, 0.5)**: No or weak bias
- **Points away from the center**: Strand-biased base composition
- **Quadrants**: The direction of deviation indicates which base is favored (e.g. higher G than C, or higher A than T)
- **Axis ranges**: Both axes share the same, center-symmetric range so the plot is symmetric left/right and top/bottom

## Style Settings Reference

### Title & Labels Tab
- **Title**: Graph title font size (8-24) and color
- **X-axis Label**: X-axis title font size (8-20) and color
- **Y-axis Label**: Y-axis title font size (8-20) and color
- **X-axis Tick Labels**: X-axis numbers font size (6-18) and color
- **Y-axis Tick Labels**: Y-axis numbers font size (6-18) and color

### Reference Lines Tab
- **Visibility**: Toggle the neutral crosshair display
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

### Legend Tab
- **Visibility**: Toggle legend display
- **Font Size**: 6 to 16
- **Text Color**: Customizable

### Gene Labels Tab
- **Visibility**: Toggle gene labels for the most deviated genes
- **Number of Genes**: 3 to 20 (ranked by distance from the center)
- **Font Size**: 6 to 14
- **Label Color**: Customizable

## Troubleshooting

**Issue**: "Please select a valid input file!"
- **Solution**: Ensure you have selected a `.out` file with the correct format

**Issue**: "Could not extract valid data from file!"
- **Solution**: Verify your file has at least 5 columns and contains valid T3s (column 1), C3s (column 2), A3s (column 3) and G3s (column 4) data

**Issue**: Graph appears distorted
- **Solution**: The graph is designed to be 600×600 pixels. Adjust window size if needed.

**Issue**: Style changes not visible
- **Solution**: Click "Apply" in Style Settings. Changes only appear after clicking Apply.

## Technical Details

- **Framework**: Tkinter with matplotlib backend
- **Canvas Size**: 600×600 pixels (fixed)
- **Output Resolution**: 300 DPI for raster formats
- **Axis Range**: One shared, center-symmetric range for both axes
- **Aspect Ratio**: Equal (square) via `set_aspect('equal')`

## License

This tool is provided as-is for research and educational purposes.

## Citation

If you use this tool in your research, please cite the theoretical foundation:

- Sueoka, N. (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. *Journal of Molecular Evolution*, 40(3), 318-325.
- Lobry, J. R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution*, 13(5), 660-665.

## Contact

For questions or suggestions, please contact the developer.

## Version History

- **v1.0**: Initial release with PR2-bias plotting, style customization and gene labels
