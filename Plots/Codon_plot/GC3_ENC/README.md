# GC3-ENC Plotter

A Python-based GUI application for plotting GC3-ENC relationship graphs from `.out` files, specifically designed for codon usage bias analysis in molecular biology research.

## Features

- **Data Visualization**: Plot GC3 (GC content at third codon position) against ENC (Effective Number of Codons)
- **Theoretical Curve**: Display the theoretical ENC curve based on Wright (1990) representing no selection pressure
- **Species Annotation**: Add species names as title annotations
- **Multiple Output Formats**: Support for PDF, SVG, PNG, and EPS formats
- **Fixed Canvas Size**: Maintains square aspect ratio regardless of window size
- **Comprehensive Style Customization**: Full control over all graph elements including:
  - Title, axis labels, and tick labels (size and color)
  - Theoretical curve (line width, style, color, visibility)
  - Scatter points (size, color, transparency, edge color, visibility)
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

2. **Load Data File**:
   - Click "Browse..." to select your `.out` file
   - The file should contain columns with NC (ENC) and GC3s (GC3) data
   - NC column should be at column index 8
   - GC3s column should be at column index 9

3. **Set Species Name**:
   - Enter the species name in the "Species name" field
   - This will appear as the graph title

4. **Plot Graph**:
   - Click "Plot Graph" to generate the visualization
   - The graph displays:
     - Red dashed line: Theoretical ENC curve (no selection pressure)
     - Blue points: Observed data from your input file

5. **Customize Style**:
   - Click "Style Settings" to open the customization panel
   - Use the tabs to adjust different elements:
     - **Title & Labels**: Adjust title, axis labels, and tick labels
     - **Theoretical Curve**: Modify the red theoretical curve appearance
     - **Scatter Points**: Customize data point appearance
     - **Legend**: Control legend visibility and style
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
  - Column 8: NC (ENC values)
  - Column 9: GC3s (GC3 values)

Example:
```
Gene    CDS_len  ENC  ...  NC    GC3s
Gene1   1200     48   ...  52.3  0.45
Gene2   850      55   ...  58.7  0.38
```

Invalid data points (marked as `*****`) are automatically skipped.

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

## Troubleshooting

**Issue**: "Please select a valid input file!"
- **Solution**: Ensure you have selected a `.out` file with the correct format

**Issue**: "Could not extract valid data from file!"
- **Solution**: Verify your file has at least 10 columns and contains valid NC (column 8) and GC3s (column 9) data

**Issue**: Graph appears distorted
- **Solution**: The graph is designed to be 600×600 pixels. Adjust window size if needed.

**Issue**: Style changes not visible
- **Solution**: Click "Apply" in Style Settings. Changes only appear after clicking Apply.

## Technical Details

- **Framework**: Tkinter with matplotlib backend
- **Canvas Size**: 600×600 pixels (fixed)
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
