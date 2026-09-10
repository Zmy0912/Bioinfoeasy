# GC12-GC3 Distribution Plotter

A Python GUI application for generating GC12-GC3 distribution charts from CSV data files. The program provides interactive visualization with customizable styling options and support for high-quality vector graphics output.

## Features

- **Interactive GUI**: User-friendly tkinter-based interface
- **Data Visualization**: Plot GC12-GC3 distribution with scatter points
- **Statistical Analysis**: 
  - Neutral expectation line (theoretical 1:1 relationship)
  - Linear regression line with R² value and regression equation
- **Customizable Styling**:
  - Scatter point color, size, and transparency
  - Line color, style, and width for both neutral and regression lines
  - Annotation font size and color
  - Title font size and color
- **Flexible Title Options**: Use species name or custom title
- **High-Quality Output**: Save as PDF, SVG, EPS (vector formats) or PNG (raster)
- **Fixed Plot Size**: Consistent 8x8 inch square output regardless of window size
- **Scrollable Interface**: Navigate through controls and plot independently

## Requirements

- Python 3.6+
- tkinter (usually included with Python)
- pandas
- matplotlib
- scipy
- numpy

## Installation

1. Clone or download this repository

2. Install required dependencies:
```bash
pip install -r requirements.txt
```

Or install manually:
```bash
pip install pandas matplotlib scipy numpy
```

## Usage

### Starting the Program

Run the application from the command line:

```bash
python gc_analysis_plotter.py
```

### Loading Data

1. Click the **"Select CSV File"** button
2. Navigate to your CSV file and select it
3. The file name will be displayed in the "Current File" section
4. Required CSV format:
   - Must contain columns: `GC3 (%)` and `GC12 (%)`
   - Data should be in percentage values (0-100)

### Creating a Plot

1. Enter the **Species Name** in the input field (default: "NC01")
2. Optionally check **"Custom Title"** to use a custom title instead of the species name
3. Click **"Plot Chart"** to generate the visualization

### Saving the Plot

1. Click **"Save Plot"** to export your chart
2. Choose the desired format and location:
   - **PDF** - Best for publications and printing (vector format)
   - **SVG** - Scalable vector graphics (vector format)
   - **EPS** - Encapsulated PostScript (vector format)
   - **PNG** - Raster image format at 300 DPI

### Customizing Styles

Click **"Style Settings"** to open the style customization panel. All style changes are applied when you click **"Plot Chart"**. The style settings window can remain open while you make multiple adjustments.

#### Scatter Points Settings
- **Color**: Hex color code (e.g., "#1f77b4")
- **Size**: Point size (10-200)
- **Alpha**: Transparency level (0.1-1.0)

#### Neutral Expectation Line
- **Color**: Line color
- **Style**: Line style (`-`, `--`, `:`, `-.`)
- **Width**: Line thickness (1-5)

#### Regression Line
- **Color**: Line color
- **Style**: Line style (`-`, `--`, `:`, `-.`)
- **Width**: Line thickness (1-5)

#### Annotation Settings
- **Font Size**: Text size for R² and equation (10-25)
- **Color**: Text color

#### Title Settings
- **Font Size**: Title text size (10-30)
- **Color**: Title text color

#### Axis Settings
- **Tick Label Font Size**: Font size for axis numeric labels (8-20)

### Clearing the Plot

Click **"Clear"** to reset the application and remove the current plot and data.

## CSV File Format Example

```csv
GC3 (%),GC12 (%)
45.2,43.8
52.1,51.3
38.7,39.5
41.3,40.8
48.9,47.6
...
```

## Chart Elements

The generated plot includes:

1. **Scatter Points**: Individual data points with customizable appearance
2. **Neutral Expectation Line**: Theoretical line showing 1:1 relationship (y = x)
3. **Linear Regression Line**: Best-fit line calculated from the data
4. **Annotations**: R² value and regression equation displayed in plain text
5. **Title**: Species name or custom title at the top
6. **Axis Labels**: GC3 (%) for x-axis and GC12 (%) for y-axis
7. **Axis Tick Labels**: Numeric labels along both axes (customizable font size)
8. **Grid Lines**: Light gray grid for reference
9. **Legend**: Identifying the two lines

## Troubleshooting

### Common Issues

**"Missing required columns in CSV file"**
- Ensure your CSV contains exactly `GC3 (%)` and `GC12 (%)` columns
- Check for extra spaces or case sensitivity

**"Failed to load file"**
- Verify the file is a valid CSV format
- Check that the file is not corrupted or in use by another program

**Plot not appearing**
- Make sure you have selected a CSV file first
- Check that the CSV data contains valid numeric values
- Try clearing and reloading the data

**Style changes not visible**
- Click **"Plot Chart"** again after changing style settings
- Style settings take effect only when regenerating the plot

## Technical Details

### Plot Configuration
- Figure size: 8×8 inches (square aspect ratio)
- X-axis: GC3 (%), range 0-100
- Y-axis: GC12 (%), range 0-100
- Origin at (0, 0)

### Statistical Methods
- Linear regression calculated using `scipy.stats.linregress`
- R² = correlation coefficient squared

## License

This project is open source and available for educational and research purposes.

## Contact

For questions or issues, please contact the development team.
