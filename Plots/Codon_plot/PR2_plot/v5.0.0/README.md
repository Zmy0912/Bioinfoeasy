# PR2-Plot Plotter

A Python-based GUI application for plotting PR2-bias (Parity Rule 2) graphs from `.out` files, designed for codon usage bias analysis in molecular biology research. It reads either a **single `.out` file** or **every `.out` file of a folder** (one panel per file in a combined figure), and it is part of the same series as the `GC3-ENC Plotter`, sharing the same layout logic and interface style.

## Features

- **Data Visualization**: Plot PR2-bias scatter of `G3/(G3+C3)` (X axis) against `A3/(A3+T3)` (Y axis)
- **Neutral Reference**: Display the reference crosshair at the unbiased center `(0.5, 0.5)`
- **Symmetric Plot**: X and Y axes share one identical, center-symmetric range, and the canvas keeps an equal square aspect ratio
- **Species Annotation**: Add species names as the graph title
- **Multiple Output Formats**: Support for PDF, SVG, PNG, and EPS formats
- **Batch Folder Mode**: Point the tool at a folder and every `.out` file inside it is read and plotted at once (invalid files are skipped with a message in the console)
- **Combined Multi-Panel Figure**: One panel per input file, arranged in a grid whose number of columns is either chosen automatically or fixed by the user; all panels share one axis range and one unified legend on the right, the preview is scrollable so that large figures keep their true size, and one individual image per input file can be exported as well
- **Fixed Canvas Size**: Maintains square aspect ratio regardless of window size (the plot area grows for combined figures)
- **Outlier Gene Labels**: Automatically label the genes that lie farthest from the bulk of the data — either from the parity center `(0.5, 0.5)` or from the center of the point cloud — printing the **bare gene symbol** only, with an adjustable number of labels, font size, color and collision-aware placement
- **Gene Function Coloring**: Optionally color every scatter point by the functional category of its gene (e.g. photosystem, ATP synthase, NADH dehydrogenase, ribosomal proteins, RNA polymerase, ...), with a per-category legend, customizable category colors and per-category visibility
- **Comprehensive Style Customization**: Full control over all graph elements including:
  - Title, axis labels, and tick labels (size and color)
  - Reference lines (line width, style, color, visibility)
  - Scatter points (size, color, transparency, edge color, visibility)
  - Gene function categories (color and visibility of each category)
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

3. Download the `pr2_plot_plotter.py` file

## Usage

1. Run the program:

```bash
python pr2_plot_plotter.py
```

2. **Choose the Input Mode** (`Input mode` radio buttons):
   - **Single file** (default): click "Browse File..." and select one `.out` file (it must contain the codon usage columns T3s, C3s, A3s and G3s)
   - **Folder (all .out files, combined figure)**: click "Browse Folder..." and select a folder; every `*.out` file inside it is read and drawn as one panel of a combined figure (files without valid data are skipped and listed in the console)
   - **Subplot columns (0 = auto)**: number of panel columns of the combined figure — `0` picks up to 4 columns automatically, any other value (1–12) forces that layout

3. **Set Species Name** (single-file mode only):
   - Enter the species name in the "Species name" field; it becomes the graph title
   - In folder mode the panels are titled with the file names instead

4. **Plot Graph**:
   - Click "Plot Graph" to generate the visualization
   - The graph displays:
     - Red dashed crosshair: the neutral reference at `(0.5, 0.5)`
     - Blue points: observed data from your input file
     - Dark-green gene symbols next to the most outlying points (see the next step)
   - In folder mode one combined figure is produced: all panels share the same axis range (so they are directly comparable) and one unified legend is placed on the right; the preview area becomes scrollable for large figures (see [Folder Mode & Combined Figure](#folder-mode--combined-figure))

5. **Label the Outlying Genes** (`Gene Labels (outliers)` row):
   - Tick / untick **"Label outlying genes"** (ticked by default)
   - **Rank by** selects how the labelled genes are chosen:
     - *Farthest from parity center (0.5, 0.5)* – the classic PR2 criterion (default)
     - *Farthest from data center* – the points that stick out of the cloud
   - **Number** = how many genes receive a label (1–30, default 5)
   - **Font size** (4–16) and **Color** of the labels
   - Every change repaints the graph immediately (no need to click "Plot Graph" again)

6. **Color Points by Gene Function** (optional):
   - Tick the checkbox **"Color points by gene function"** in the top panel
   - Each point is colored according to the functional category of its gene and the legend lists every category with its counts
   - Untick the checkbox to go back to the single-color scatter

7. **Customize Style**:
   - Click "Style Settings" to open the customization panel
   - Use the tabs to adjust different elements:
     - **Title & Labels**: Adjust title, axis labels, and tick labels
     - **Reference Lines**: Modify the neutral crosshair appearance (width, style, color, visibility)
     - **Scatter Points**: Customize data point appearance
     - **Legend**: Control legend visibility and style
     - **Gene Function**: Toggle function-based coloring, set the color of every functional category and show/hide individual categories
   - Click "Apply" to save changes (graph will redraw automatically)
   - Click "Reset to Default" to restore original settings (including the gene labels and the combined-figure options)

8. **Save Image**:
   - Click "Save Image" to export the graph (the combined figure or the single graph)
   - Choose from PDF, SVG, PNG, or EPS formats
   - Saved images maintain the fixed square aspect ratio (the combined figure keeps its grid layout)
   - In folder mode, tick **"Also export one image per file (folder mode)"** in the options row to additionally write one single-graph image per input file into an `individual_plots` folder next to the file you saved

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

## Folder Mode & Combined Figure

Selecting **Folder (all .out files, combined figure)** and choosing a folder that
contains `.out` files makes the program read every one of them (sorted by file name)
and render a single combined figure — one panel per file, which is the quickest way to
compare many genomes side by side.

| Setting | Effect |
|---------|--------|
| *Subplot columns (0 = auto)* | `0` → automatically up to 4 columns, `1–12` → that exact number of columns; the number of rows follows from the file count |
| *Also export one image per file* | When saving in folder mode, additionally write one single-graph image per input file into `individual_plots/` next to the saved figure |

**Layout of the combined figure**

- Panel size is about 4 × 4 inches, so the figure grows with the grid; the preview area
  is a scrollable viewport that shows the figure at its true size instead of shrinking it.
- Every panel uses the **same axis range** (the center-symmetric range that covers all
  files) and the same square aspect ratio, so panels are directly comparable.
- Each panel carries its own title (the file name), axis labels, tick labels, grid and
  gene labels, but **no in-panel legend**: one unified legend is placed on the right side
  of the figure. It lists the neutral reference line plus either "Observed data" or, with
  gene function coloring enabled, all visible functional categories.
- The status bar summarises the result, e.g. `Combined figure: 20 file(s), 1520 data points, 5×4 grid`.

**Individual panels**

The optional `individual_plots` export renders every file on its own, titled with the
file name and with its own in-panel legend, so the overview figure can be complemented
with publication-ready single graphs.

## Outlier Gene Labels

With **Label outlying genes** (ticked by default) the points lying farthest from the
bulk of the data are marked with the name of their gene, so the most biased or most
unusual genes can be found without opening the data table.

**Which genes get a label**

| Rank by | Score |
|---------|-------|
| *Farthest from parity center (0.5, 0.5)* (default) | Distance from the unbiased parity point, i.e. `hypot(G3/(G3+C3) − 0.5, A3/(A3+T3) − 0.5)` — the classic PR2 selection criterion |
| *Farthest from data center* | Distance from the center of the point cloud (mean of all points) |

Both scores are divided by the visible axis span (the shared, center-symmetric
range), so they are fractions of the plotted area and the two methods are directly
comparable. The **Number** field defines how many genes receive a label (the most
extreme ones first). Labels are placed automatically: several candidate offsets
around the point are tried and the first one that does not collide with an already
placed label is used, connected to its point by a thin arrow. A small number (2–5)
is recommended when the gene-function legend is enabled as well.

**Which part of the identifier is used as the gene name**

Only the gene symbol is printed — the location, the orientation and the file counter
that the pipeline appends to the header are dropped. Everything in front of the last
`>` of a FASTA style header is removed as well, and the symbol ends at the first
character that cannot be part of a gene name:

| Raw identifier in the `.out` file | Label shown |
|-----------------------------------|-------------|
| `rpl20_[69725:70078](-)\|01` | `rpl20` |
| `>psbA_[73610:75136](+)\|01.` | `psbA` |
| `atpF_join{[13161:13305](-` | `atpF` |
| `psbK_[7404:7589](+)\|01.fa` | `psbK` |
| `trnL-UAA_[107000:107063]\|01.` | `trnL-UAA` |
| `>rbcL` | `rbcL` |

The terminators are `_` `(` `[` `{` `:` `|` `,` `;` `.` and the space character,
while letters, digits and `-` are kept (so tRNA names such as `trnH-GUG` stay intact).

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

### Gene Function Tab
- **Color scatter points by gene function**: Enable / disable function-based coloring
- **Show (per category)**: Hide a single category from the graph without disabling the whole feature
- **Color (per category)**: Choose the color used for each functional category

## Troubleshooting

**Issue**: "Please select a valid input file!"
- **Solution**: Ensure you have selected a `.out` file with the correct format

**Issue**: "Could not extract valid data from file!"
- **Solution**: Verify your file has at least 5 columns and contains valid T3s (column 1), C3s (column 2), A3s (column 3) and G3s (column 4) data

**Issue**: "Please select a valid input folder!"
- **Solution**: Check that the *Input folder* field points to an existing directory (use "Browse Folder..."), and that folder mode is the selected input mode

**Issue**: "No .out files found in the selected folder!"
- **Solution**: The folder must contain files ending in `.out`; sub-folders are not searched

**Issue**: Only some of the files appear in the combined figure
- **Solution**: Files without a single valid data row are skipped on purpose; their names are printed in the console window

**Issue**: The combined figure is hard to read
- **Solution**: Use the scrollbars of the graph area to move around; set *Subplot columns* to a small value for bigger panels or a larger value to fit more panels per row

**Issue**: Graph appears distorted
- **Solution**: The graph is designed to be 600×600 pixels. Adjust window size if needed.

**Issue**: The gene labels are missing
- **Solution**: Tick **"Label outlying genes"** in the `Gene Labels (outliers)` row and make sure **Number** is at least 1.

**Issue**: Style changes not visible
- **Solution**: Click "Apply" in Style Settings. Changes only appear after clicking Apply.

## Technical Details

- **Framework**: Tkinter with matplotlib backend
- **Input Modes**: a single `.out` file, or every `.out` file of a folder (sorted by name)
- **Canvas Size**: 600×600 pixels for a single graph; the container grows to the figure's true pixel size for combined figures (scrollable viewport)
- **Output Resolution**: 300 DPI for raster formats
- **Axis Range**: One shared, center-symmetric range for both axes (one range common to all panels of a combined figure)
- **Aspect Ratio**: Equal (square) via `set_aspect('equal')`
- **Combined Figure Layout**: one panel per file, about 4 in per panel, automatic (≤ 4) or user-defined column count, unified side legend instead of per-panel legends
- **Individual Export**: optional `individual_plots/` folder containing one image per input file
- **Gene Label Ranking**: distances are normalised by the visible axis span, so both ranking methods return comparable scores (fractions of the plotted area)

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
- **v1.1**: Added gene function classification coloring (per-category colors, per-category visibility, category legend)
- **v1.2**: Reworked the gene labels — the old fixed "distance from the center" annotation (and its Style Settings tab) was removed and replaced by a dedicated `Gene Labels (outliers)` row in the main panel, with two ranking methods (parity center / data center), a free number of labels (1–30), adjustable font size and color, collision-aware label placement, and printing of the bare gene symbol only (location, orientation and file suffix stripped, `>` FASTA headers supported)
- **v1.3**: Added batch folder mode with a combined multi-panel figure — the input panel now offers Single file / Folder radio buttons, a folder picker, an optional fixed subplot-column count, a scrollable preview that keeps large figures at their true size, one shared axis range and a unified side legend for all panels, and an option to export one individual image per input file
