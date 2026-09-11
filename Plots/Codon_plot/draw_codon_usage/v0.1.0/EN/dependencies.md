# Codon Usage Analysis Tool - Dependencies

## Project Overview

This is a Python tool for visualizing codon usage frequency, capable of generating stacked bar charts to display the RSCU value distribution of codons corresponding to different amino acids. It supports multiple color schemes, custom text formatting, chart size adjustments, and other advanced features, with independent display zooming and panning capabilities.

## Key Features

- **Data Import**: Supports Excel files in .xlsx and .xls formats
- **Visualization**: Generates interactive stacked bar charts displaying codon RSCU values
- **Multiple Color Schemes**: Offers 11 color schemes (tab20c, tab20b, tab20, Set1, Set2, Set3, Pastel1, Pastel2, viridis, plasma, inferno)
- **Text Customization**: Adjustable font sizes for codon labels, axis labels, tick labels, and legend text
- **Chart Customization**: Adjustable canvas size, block height/width, border thickness, and block distance from X-axis
- **Interactive Operations**: Supports zoom, pan, and other operations
- **Multi-format Export**: Supports PDF, SVG, and PNG format exports

## Python Version Requirements

- Python 3.7 or higher
- Recommended version: Python 3.9 - 3.12

## Required Dependencies

### 1. pandas
- **Version**: >= 1.3.0
- **Purpose**: Reading and processing Excel data
- **Installation Command**:
  ```bash
  pip install pandas
  ```

### 2. matplotlib
- **Version**: >= 3.3.0
- **Purpose**: Drawing charts, creating graphical interfaces
- **Installation Command**:
  ```bash
  pip install matplotlib
  ```

### 3. numpy
- **Version**: >= 1.19.0
- **Purpose**: Numerical computation, array operations
- **Installation Command**:
  ```bash
  pip install numpy
  ```

### 4. openpyxl
- **Version**: >= 3.0.0
- **Purpose**: Reading Excel (.xlsx) files
- **Installation Command**:
  ```bash
  pip install openpyxl
  ```

### 5. xlrd
- **Version**: >= 1.2.0
- **Purpose**: Reading legacy Excel (.xls) files (optional)
- **Installation Command**:
  ```bash
  pip install xlrd
  ```
- **Note**: xlrd 2.0+ no longer supports .xlsx format, only for compatibility with old files

## Built-in Modules (No Installation Required)

The following are Python standard libraries, no additional installation needed:

- **tkinter** - GUI framework (included with Python)
- **tkinter.ttk** - Enhanced Tkinter components (included with Python)
- **tkinter.filedialog** - File selection dialog (included with Python)
- **tkinter.messagebox** - Message prompt box (included with Python)
- **os** - Operating system interface (included with Python)
- **sys** - System-related parameters and functions (included with Python, for command-line argument processing)
- **matplotlib.backends.backend_tkagg** - Matplotlib's Tkinter backend (installed with matplotlib)
  - **FigureCanvasTkAgg** - Embeds Matplotlib figures in Tkinter
  - **NavigationToolbar2Tk** - Provides zoom, pan, and other toolbar functions

## System Requirements

### Windows
- Windows 7 or higher
- System should have Chinese fonts installed (e.g., SimHei, Microsoft YaHei) to support Chinese display

### macOS
- macOS 10.13 or higher
- Need to install Chinese font support

### Linux
- Any mainstream Linux distribution
- May need to install tkinter: `sudo apt-get install python3-tk`

## Complete Installation Commands

### Install all dependencies at once using pip

```bash
pip install pandas matplotlib numpy openpyxl xlrd
```

### Install using requirements.txt

If a requirements.txt file is provided, use the following command:

```bash
pip install -r requirements.txt
```

## requirements.txt Content

```
pandas>=1.3.0
matplotlib>=3.3.0
numpy>=1.19.0
openpyxl>=3.0.0
xlrd>=1.2.0,<2.0.0
```

**Note**: xlrd 2.0.0+ no longer supports .xlsx format. If you need to read legacy .xls files, use version 1.2.0. For .xlsx files, the program will prioritize using the openpyxl engine.

## Verify Installation

Run the following command to verify all dependencies are installed correctly:

```bash
python -c "import pandas; import matplotlib; import numpy; import openpyxl; import xlrd; import tkinter; print('All dependencies installed successfully!')"
```

Or run the program's built-in test:

```bash
python draw_codon_usage.py --debug
```

**Note**: The `--debug` parameter will read the `新建文件夹/A1.xlsx` file and display debug information. Ensure the file exists before using this mode.

## Troubleshooting

### 1. tkinter Related Errors

**Error Message**: `ModuleNotFoundError: No module named 'tkinter'`

**Solutions**:
- **Windows**: Reinstall Python and check "tcl/tk and IDLE" option
- **macOS**: `brew install python-tk`
- **Linux**: `sudo apt-get install python3-tk`

### 2. openpyxl Related Errors

**Error Message**: `ModuleNotFoundError: No module named 'openpyxl'`

**Solution**:
```bash
pip install openpyxl
```

### 3. Chinese Font Display Issues

**Error Message**: Chinese characters in charts display as squares or garbled text

**Solutions**:
- **Windows**: Ensure system has SimHei or Microsoft YaHei fonts installed
- **macOS**: Install Chinese font package
- **Linux**: Install Chinese fonts: `sudo apt-get install fonts-wqy-zenhei`

### 4. Excel File Read Failures

**Error Message**: `FileNotFoundError` or `PermissionError`

**Solutions**:
- Confirm file path is correct
- Ensure file is not being used by other programs
- Try saving file as .xlsx format

## Usage Methods

### Graphical Interface Mode (Recommended)

Run the program directly to launch the graphical interface:

```bash
python draw_codon_usage.py
```

In the graphical interface:
1. Click "Select Excel File" button to import data
2. Adjust color scheme, canvas size, text format, and other parameters
3. Click "Redraw" to update the chart
4. Use navigation toolbar for zoom, pan, and other operations
5. Click save buttons to export charts

### Command Line Mode

```bash
python draw_codon_usage.py <Excel_file_path> [output_file_path]
```

Example:
```bash
python draw_codon_usage.py 新建文件夹/A1.xlsx codon_usage_chart.pdf
```

### Debug Mode

```bash
python draw_codon_usage.py --debug
```

Debug mode will display detailed information about the Excel file, including column names, data types, first 10 rows, etc.

**Note**: Debug mode defaults to reading the `新建文件夹/A1.xlsx` file.

## Data Format Requirements

Excel files should contain the following columns (supporting English or Chinese column names):

| English Name | Chinese Name | Description |
|--------------|--------------|-------------|
| Codon | 密码子 | Codon sequence (e.g., UUU, AUG, etc.) |
| Frequency | 频数/数量 | Codon occurrence frequency |
| RSCU | RSCU | Relative Synonymous Codon Usage |
| AminoAcid | 氨基酸 | Corresponding amino acid (optional, will be automatically derived if missing) |

The program will automatically recognize column names. If column name format differs, ensure data columns are arranged in the above order.

**Supported Codon Formats**: The program uses RNA bases (U), also compatible with DNA bases (T), will automatically convert.

**Amino Acid Sorting**: The program arranges 20 amino acids in standard order: Phe, Leu, Ile, Met, Val, Ser, Pro, Thr, Ala, Tyr, His, Gln, Asn, Lys, Asp, Glu, Cys, Trp, Arg, Gly

## Input Table Examples

### Example 1: Complete Format (with Amino Acid Column)

| Codon | Frequency | RSCU | AminoAcid |
|-------|-----------|------|-----------|
| UUU   | 150       | 1.20 | Phe       |
| UUC   | 100       | 0.80 | Phe       |
| UUA   | 80        | 0.80 | Leu       |
| UUG   | 120       | 1.20 | Leu       |
| CUU   | 100       | 1.00 | Leu       |
| CUC   | 100       | 1.00 | Leu       |
| CUA   | 50        | 0.50 | Leu       |
| CUG   | 150       | 1.50 | Leu       |
| AUU   | 120       | 1.20 | Ile       |
| AUC   | 80        | 0.80 | Ile       |
| AUA   | 50        | 0.50 | Ile       |
| AUG   | 100       | 1.00 | Met       |
| GUU   | 100       | 1.00 | Val       |
| GUC   | 100       | 1.00 | Val       |
| GUA   | 80        | 0.80 | Val       |
| GUG   | 120       | 1.20 | Val       |
| ...   | ...       | ...  | ...       |

### Example 2: Simplified Format (without Amino Acid Column)

| Codon | Frequency | RSCU |
|-------|-----------|------|
| UUU   | 150       | 1.20 |
| UUC   | 100       | 0.80 |
| UUA   | 80        | 0.80 |
| UUG   | 120       | 1.20 |
| AUU   | 120       | 1.20 |
| AUC   | 80        | 0.80 |
| AUG   | 100       | 1.00 |
| ...   | ...       | ...  |

**Note**: If there is no amino acid column, the program will automatically derive the corresponding amino acid based on the codon.

### Example 3: DNA Format (using T base)

| Codon | Frequency | RSCU |
|-------|-----------|------|
| TTT   | 150       | 1.20 |
| TTC   | 100       | 0.80 |
| TTA   | 80        | 0.80 |
| TTG   | 120       | 1.20 |
| ATT   | 120       | 1.20 |
| ATC   | 80        | 0.80 |
| ATG   | 100       | 1.00 |
| ...   | ...       | ...  |

**Note**: The program will automatically convert DNA bases (T) to RNA bases (U) for processing.

### Data Description

- **Frequency**: Actual occurrence count of the codon in the sequence
- **RSCU (Relative Synonymous Codon Usage)**: Indicator measuring codon usage bias
  - RSCU = 1.0: No bias (equal usage)
  - RSCU > 1.0: Preference for this codon
  - RSCU < 1.0: Dispreference for this codon
- Each amino acid typically corresponds to 2-6 codons (except for Met and Trp)

## Interface Function Description

### Toolbar Functions

**First Row - File Operations**
- Select Excel File: Import data file
- Save as PDF: Export as vector PDF format
- Save as SVG: Export as editable SVG format
- Save as PNG: Export as high-resolution PNG image

**Second Row - Chart Settings**
- Color Scheme: Choose from 11 color schemes
- Show Codon Names: Control whether codon labels are displayed in the bar chart
- Show Stacked Blocks: Control whether bottom codon stacked block legend is displayed
- Canvas Size: Adjust chart width and height (default 18x13)
- Block Height: Adjust bottom stacked block height (default 0.12, range 0.05-0.3)
- Block Width: Adjust bottom stacked block width (default 0.6, range 0.1-1.0)
- Redraw: Apply all parameter settings and redraw chart

**Third Row - Text and Border Settings**
- Bar Chart Text Size: Adjust codon name size in bar chart (default 10, range 6-24)
- Axis Label Size: Adjust X-axis and Y-axis label sizes (default 14, range 8-30)
- Tick Label Size: Adjust X-axis and Y-axis tick label sizes (default 12, range 6-20)
- Border Thickness: Adjust chart border and tick line thickness (default 1.0, range 0.5-5.0)
- Block Position: Adjust distance between bottom stacked block and X-axis (default 0.22, range 0.05-0.5)
- Block Text Size: Adjust text size in bottom stacked block (default 8, range 4-16)

**Drawing Area**
- Displays generated codon usage frequency chart
- Supports mouse wheel zoom and drag pan
- Navigate through large charts using scrollbars

### Navigation Toolbar Functions

The drawing area above provides a navigation toolbar with the following functions:
- Home: Return to original view
- Back/Forward: Navigate through view history
- Pan: Pan tool (click and drag chart)
- Zoom: Zoom tool (select area to zoom)
- Save: Save current view
- Subplots: Adjust subplot parameters
- Configure: Configure plot parameters

## Output Format Support

The program supports the following output formats:

- **PDF** - Vector format, editable (Adobe Illustrator, Inkscape, etc.), 300 DPI
- **SVG** - Vector format, fully editable, 300 DPI
- **PNG** - High-resolution bitmap, 300 DPI, suitable for presentations and web use

## Performance Recommendations

- For large datasets (>1000 rows), use newer matplotlib version (>= 3.5)
- If memory issues occur, try batch processing data or using smaller chart sizes
- Use scrollbars and navigation toolbar to efficiently navigate large charts
- Recommended to use PDF or SVG formats for best editing and print quality

## Development Environment Recommendations

### Recommended IDEs
- Visual Studio Code
- PyCharm
- Jupyter Notebook (for debugging)

### Virtual Environment (Recommended)

Use virtual environment to avoid dependency conflicts:

```bash
# Create virtual environment
python -m venv venv

# Activate virtual environment
# Windows:
venv\Scripts\activate
# macOS/Linux:
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

## Contact Support

If you have any questions or suggestions, please contact the developer or refer to the project documentation.

## FAQ

### Q: Chart not displaying after generation?
A: Ensure Excel file contains correct column names (Codon, RSCU, etc.), you can run `--debug` mode to check file contents.

### Q: How to modify chart color scheme?
A: Select a different color scheme from the "Color Scheme" dropdown in the second row of the graphical interface, then click "Redraw".

### Q: How to adjust exported chart quality?
A: The program defaults to 300 DPI export. If you need higher or lower resolution, you can modify the `dpi` parameter in the code (located in save_as_pdf, save_as_png methods).

### Q: How to adjust codon label display?
A: Check/uncheck "Show Codon Names" checkbox in the second row of the graphical interface, or adjust "Bar Chart Text Size" in the third row.

### Q: What codon encoding formats are supported?
A: The program supports RNA encoding (U) and DNA encoding (T). For example, both UUU and TTT will be recognized as phenylalanine codons.

### Q: Why doesn't the chart resize with the window?
A: This is by design. The chart is displayed at a fixed size with scroll and zoom support for better viewing of chart details. You can use the navigation toolbar's zoom function or scrollbars to navigate large charts.

### Q: How to install tkinter on Linux?
A: Install using package manager, for example:
- Ubuntu/Debian: `sudo apt-get install python3-tk`
- Fedora: `sudo dnf install python3-tkinter`
- Arch: `sudo pacman -S tk`

## Changelog

### v1.1 (March 2026)
- Added independent display feature, chart no longer scales with window
- Added scrollbar support, supports up/down/left/right scrolling
- Added navigation toolbar, supports zoom and pan operations
- Added codon stacked block width adjustment feature
- Added show/hide codon stacked block feature
- Added text size adjustment features (codon labels, axis labels, tick labels, legend)
- Added Y-axis tick label size adjustment
- Added chart border thickness adjustment
- Added stacked block distance from X-axis adjustment
- Removed all font bold features

### v1.0 (Initial Version)
- Basic codon usage frequency visualization function
- Support for multiple color schemes
- Support for Excel data import
- Support for PDF/SVG/PNG export

---

**Last Updated**: March 2026
**Document Version**: 1.1
**Program Version**: 1.1
