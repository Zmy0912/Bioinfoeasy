# GC3-ENC Plotter

A Python-based GUI application for plotting GC3-ENC relationship graphs from `.out` files, specifically designed for codon usage bias analysis in molecular biology research.

## Features

- **Data Visualization**: Plot GC3 (GC content at third codon position) against ENC (Effective Number of Codons)
- **Theoretical Curve**: Display the theoretical ENC curve based on Wright (1990) representing no selection pressure
- **Species Annotation**: Add species names as title annotations
- **Two Input Modes**:
  - **Single file**: analyze one `.out` file at a time (original behavior)
  - **Folder (batch)**: select a folder, and **every `.out` file inside it** is parsed and plotted automatically
- **Combined Figure**: All files of a folder are drawn into **one combined grid figure**, with the panel titles defaulting to the **file names** and the layout (number of columns) configurable or automatic
- **Combined Figure Legend**: Each panel keeps its own theoretical curve and scatter points, while the **legend is drawn once on the right side of the whole figure** (instead of being repeated inside every panel)
- **Batch Export**: When saving the combined figure, optionally also export one image per input file into an `individual_plots` sub-folder
- **Multiple Output Formats**: Support for PDF, SVG, PNG, and EPS formats
- **Fixed Canvas Size**: Maintains square aspect ratio regardless of window size
- **Gene Function Coloring**: Optionally color every scatter point by the functional category of its gene (e.g. photosystem, ATP synthase, NADH dehydrogenase, ribosomal proteins, RNA polymerase, ...), with a per-category legend, customizable category colors and per-category visibility
- **Outlier Gene Labels**: Optionally write the **gene symbol** next to the points that lie farthest from the bulk of the data, so that the genes with an unusual codon usage can be identified at a glance (number of labels, ranking method, font size and label color are adjustable)
- **Comprehensive Style Customization**: Full control over all graph elements including:
  - Title, axis labels, and tick labels (size and color)
  - Theoretical curve (line width, style, color, visibility)
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

3. Download the `gc3_enc_plotter.py` file

## Usage

1. Run the program:

```bash
python gc3_enc_plotter.py
```

2. **Choose the Input Mode** (top of the window):
   - **Single file**: analyze one `.out` file
   - **Folder (all .out files, combined figure)**: analyze every `.out` file in a folder and combine them into one figure

3. **Load Data**:
   - *Single-file mode*: click "Browse File..." to select your `.out` file
   - *Folder mode*: click "Browse Folder..." to select the folder containing `.out` files
     (the default path points to the `out` folder next to this script, if it exists)
   - The file should contain columns with NC (ENC) and GC3s (GC3) data
   - NC column should be at column index 8
   - GC3s column should be at column index 9

4. **Configure the Combined Figure (folder mode only)**:
   - **Subplot columns (0 = auto)**: number of panels per row. `0` chooses automatically (up to 4 per row)
   - Optionally tick **"Also export one image per file (folder mode)"** to also write an image per input file when saving

5. **Set Species Name** (single-file mode):
   - Enter the species name in the "Species name" field
   - This will appear as the graph title

6. **Plot Graph**:
   - Click "Plot Graph" to generate the visualization
   - *Single-file mode*: one square graph is drawn; the graph displays:
     - Red dashed line: Theoretical ENC curve (no selection pressure)
     - Blue points: Observed data from your input file
   - *Folder mode*: one combined grid figure is drawn; the graph displays:
     - One panel per `.out` file, titled with the **file name**
     - Each panel keeps its own theoretical curve and scatter points
     - The **shared legend** is placed on the **right side of the whole figure**, so it is drawn only once instead of in every panel

7. **Label the Outlying Genes** (optional, `Gene Labels (outliers)` row):
   - Tick **"Label outlying genes"**
   - **Rank by** selects how the genes that deserve a label are chosen:
     - *Farthest from data center* – distance to the center of the point cloud (default)
     - *Farthest from theoretical curve* – |ENC observed − ENC expected| of the Wright curve
   - **Number** = how many genes per panel get a label (1–30, default 5)
   - **Font size** (4–16) and **Color** of the labels
   - Every change redraws the figure immediately, in single-file **and** folder mode

8. **Color Points by Gene Function** (optional):
   - Tick the checkbox **"Color points by gene function"** in the top panel
   - Each point is colored according to the functional category of its gene and the legend lists every category with its counts
   - In folder mode the categories of all files are shared and listed in the **side legend**
   - Untick the checkbox to go back to the single-color scatter

9. **Customize Style**:
   - Click "Style Settings" to open the customization panel
   - Use the tabs to adjust different elements:
     - **Title & Labels**: Adjust title, axis labels, and tick labels
     - **Theoretical Curve**: Modify the red theoretical curve appearance
     - **Scatter Points**: Customize data point appearance
     - **Legend**: Control legend visibility and style
     - **Gene Function**: Toggle function-based coloring, set the color of every functional category and show/hide individual categories
   - Click "Apply" to save changes (graph will redraw automatically)
   - Click "Reset to Default" to restore original settings

10. **Save Image**:
   - Click "Save Image" to export the graph
   - Choose from PDF, SVG, PNG, or EPS formats
   - *Folder mode*: the combined figure is saved; if **"Also export one image per file"** is ticked, one image per input file is additionally written to an `individual_plots` sub-folder next to the saved file
   - Saved images maintain the fixed square aspect ratio

## Input File Format

The input `.out` file should be a whitespace-separated text file with the following structure:

- First line: Header (automatically skipped)
- Subsequent lines: Data rows
- Required columns (0-indexed):
  - Column 0: Gene name / identifier (used for gene function classification)
  - Column 8: NC (ENC values)
  - Column 9: GC3s (GC3 values)

Example:
```
Gene    CDS_len  ENC  ...  NC    GC3s
Gene1   1200     48   ...  52.3  0.45
Gene2   850      55   ...  58.7  0.38
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

The theoretical ENC curve is calculated using the Wright (1990) formula:

```
ENC = 2 + GC3 + 29 / (GC3² + (1-GC3)²)
```

This curve represents the expected ENC values under no selection pressure, serving as a reference for evaluating codon usage bias. Points below the curve indicate selection pressure favoring preferred codons.

## Graph Interpretation

- **Points on or above the curve**: No or weak selection pressure
- **Points below the curve**: Selection pressure favoring specific codons
- **Y-axis range**: 20 to 64.05 (ENC values)
- **X-axis range**: 0 to 1 (GC3 proportion)

## Folder Mode & the Combined Figure

Folder mode turns the tool into a batch analyzer:

1. Select **Folder (all .out files, combined figure)** as the input mode.
2. Click **Browse Folder...** and pick a directory — every `*.out` file in it is parsed.
3. Optionally set **Subplot columns** (`0` = auto, up to 4 panels per row).
4. Click **Plot Graph**.

The result is a single grid figure where:

- Each panel is one `.out` file and its **title is the file name** by default.
- Every panel keeps its own theoretical ENC curve and scatter points (same axis ranges as the single-file view: GC3 `0–1`, ENC `20–64.05`).
- The **legend is drawn once on the right side of the whole figure** (theoretical curve, plus either "Observed data" or every gene functional category) instead of being repeated inside all panels.
- The figure size grows with the number of panels so that details stay readable (`~4.0 in` per panel, plus a `3.0 in` strip on the right).
- The in-window preview is placed in a **scrollable area** and shows the figure at its **true size**, so large combined figures can be inspected without squeezing.

Files that contain no valid data are skipped and reported in the console.

When saving in folder mode you can tick **"Also export one image per file (folder mode)"** to additionally write a full-size single-file figure for every input file into an `individual_plots` sub-folder.

## Outlier Gene Labels

With **Label outlying genes** the points that sit far away from the rest of the data
are marked with the name of their gene, so the most interesting (or most suspicious)
genes can be spotted without opening the data table.

**Which genes get a label**

| Rank by | Score |
|---------|-------|
| *Farthest from data center* (default) | Euclidean distance from the center of the point cloud, computed in units of the visible axes (GC3 and ENC are both normalised) |
| *Farthest from theoretical curve* | `\|ENC(observed) − ENC(Wright curve at the same GC3)\|`, i.e. how far the gene deviates from the no-selection reference |

The **Number** field defines how many genes per panel receive a label (the most
extreme ones first). Labels are placed automatically: several candidate offsets are
tested around the point and the first one that does not collide with an already
placed label is used, connected to its point by a thin arrow.

The same settings apply to the single-file view, to every panel of the combined
folder figure, and to the images written by **Save Image** / the per-file export.
For large combined figures a small number (2–5) is recommended.

**Which part of the identifier is used as the gene name**

Only the gene symbol is printed — the location, the orientation and the file counter
that the pipeline appends to the header are dropped. Everything in front of the last
`>` of a FASTA style header is removed as well, and the symbol ends at the first
character that cannot be part of a gene name:

| Raw identifier in the `.out` file | Label shown |
|-----------------------------------|-------------|
| `rpl20_[69725:70078](-)|01` | `rpl20` |
| `>psbA_[73610:75136](+)|01.` | `psbA` |
| `atpF_join{[13161:13305](-` | `atpF` |
| `psbK_[7404:7589](+)|01.fa` | `psbK` |
| `trnL-UAA_[107000:107063]|01.` | `trnL-UAA` |
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

### Theoretical Curve Tab
- **Visibility**: Toggle curve display
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
- **Solution**: Verify your file has at least 10 columns and contains valid NC (column 8) and GC3s (column 9) data

**Issue**: Graph appears distorted
- **Solution**: The single-file graph is designed to be 600×600 pixels. Adjust window size if needed. Combined figures are shown 1:1 inside a scrollable area, so use the scrollbars to inspect them.

**Issue**: The gene labels are missing
- **Solution**: Tick **"Label outlying genes"** in the `Gene Labels (outliers)` row — the feature is off by default — and make sure **Number** is at least 1. Plot the graph afterwards (or change any label setting, which redraws automatically).

**Issue**: Style changes not visible
- **Solution**: Click "Apply" in Style Settings. Changes only appear after clicking Apply.

## Technical Details

- **Framework**: Tkinter with matplotlib backend
- **Canvas Size**: 600×600 pixels for single-file mode; for the combined figure the graph area is a **scrollable viewport** that renders the figure at its **true full size** (so nothing is squeezed into the window)
- **Layout**: Single file → one `6×6 in` square figure; folder mode → grid figure of `~4.0 in` per panel plus a `3.0 in` side strip for the shared legend
- **Output Resolution**: 300 DPI for raster formats
- **Y-axis**: 20 to 64.05 (includes 5% top padding)
- **X-axis**: 0 to 1 (full GC3 range)

## License

This tool is provided as-is for research and educational purposes.

## Citation

If you use this tool in your research, please cite the theoretical foundation:

Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23-29.

## Contact

For questions or suggestions, please contact the developer.

## Version History

- **v1.0**: Initial release with basic plotting functionality
- **v1.1**: Added comprehensive style customization
- **v1.2**: Added X/Y axis tick label size control
- **v1.3**: Added gene function classification coloring (per-category colors, per-category visibility, category legend)
- **v1.4**: Added folder (batch) mode — every `.out` file in a folder is plotted into one combined grid figure with file-name titles, a configurable column count, a single side legend, a scrollable full-size preview, and optional per-file image export
- **v1.5**: Added outlier gene labels — the genes farthest from the data center (or from the theoretical curve) are printed next to their point with the bare gene symbol (location/orientation/suffix stripped, `>` FASTA headers supported), with adjustable number, font size, color and collision-aware placement
