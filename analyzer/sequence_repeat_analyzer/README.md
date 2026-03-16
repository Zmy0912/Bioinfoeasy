# Sequence Repeat Analyzer

A graphical application for identifying and analyzing scattered repeats in sequence files.

## Features

### Core Functionality
- **File Recognition**: Automatically identifies and parses sequence repeat text files
- **Batch Analysis**: Supports analyzing multiple text files in a folder
- **Multi-dimensional Statistics**: Provides various statistical dimensions and visualization charts

### Statistical Analysis
- **Total Repeats**: Counts total repeats for each file
- **By Type**: Statistics by repeat type (P/F/R)
- **By Length**: Statistics grouped by repeat sequence length
- **By Length Range**: Statistics by length range (30-39, 40-49, 50-59, >=60 bp)

### Graphical Interface
- **Control Panel**: Folder selection, file list, real-time statistics
- **Chart Display**: Scrollable visualization area with interactive charts
- **Navigation Toolbar**: Zoom, pan, save charts and other interactive features

### Visualization Charts
- Total repeats comparison bar chart
- By type comparison stacked bar chart
- By length range comparison stacked bar chart
- File type distribution pie chart

### Data Display
- **Details**: Statistics by length and type
- **Comparison**: Cross-file comparison analysis
- **Type Statistics**: Detailed statistics of different repeat types within each file

## Requirements

- Python 3.6 or higher
- Windows/Linux/macOS operating system

## Installation

```bash
pip install -r requirements.txt
```

### Dependencies
- matplotlib >= 3.5.0
- (Other dependencies are part of Python standard library)

## Usage

### 1. Run the Program

```bash
python sequence_repeat_analyzer_en.py
```

### 2. Operation Steps

1. **Select Folder**
   - Click "Browse..." button to select a folder containing sequence repeat files
   - Or directly type the folder path in the text box

2. **Analyze Files**
   - Click "Analyze Folder" button to start analysis
   - The program will automatically identify all `.txt` format sequence files in the folder

3. **View Statistics**
   - The left "Statistics" panel displays statistics for the currently selected file
   - Click different files in "File List" to switch views

4. **View Charts**
   - The right "Charts" area displays visualization charts
   - Use scrollbars to browse charts (up/down/left/right)
   - Use the toolbar for zoom, pan, save, and other operations

5. **View Detailed Data**
   - Click bottom tabs to switch views:
     - **Details**: Detailed statistics by length and type
     - **Comparison**: Comparison statistics across multiple files
     - **Type Statistics**: Detailed distribution of different repeat types within each file

## File Format

The program supports parsing text files in the following format:

```
# Comment line (starts with #)
26553  85030 P 26553 130233  0 0.00e+00
   52  29749 P    52  29749  0 3.41e-22
   49  90331 F    49  90352 -1 3.21e-18
   49  90331 P    49 151415 -1 3.21e-18
   49 151415 F    49 151436 -1 3.21e-18
```

### Field Description (by column)
1. **Repeat Length**: Number of base pairs in the repeat sequence
2. **Start Position 1**: Position of the first repeat
3. **Repeat Type**:
   - `P` (Perfect): Perfect match (no errors)
   - `F` (Forward): Forward match
   - `R` (Reverse): Reverse match
4. **Start Position 2**: Position of the second repeat
5. **End Position**: End position of the second repeat
6. **Errors**: Number of matching errors
7. **E-value**: Expected value (statistical significance)

## Chart Features

### Navigation Toolbar
- **Home**: Reset view to initial state
- **Back/Forward**: Navigate view history
- **Pan**: Pan the chart
- **Zoom**: Zoom in/out of the chart
- **Save**: Save chart as image file

### Scroll Operations
- Use vertical scrollbar to scroll up/down
- Use horizontal scrollbar to scroll left/right
- Support mouse wheel vertical scrolling

## Technical Features

### English Support
- Full English interface
- Standard fonts for better compatibility
- Clear English labels and titles

### Performance
- Efficient file parsing algorithms
- Support for large-scale data processing
- Smooth graphical interface response

### Code Structure
- Modular design for easy maintenance and extension
- Clear code comments
- Follows Python coding standards

## FAQ

### Q: What should I do if charts are too small to see clearly?
A: Use the zoom feature in the bottom toolbar or mouse wheel to zoom in/out.

### Q: What file formats are supported?
A: Currently supports `.txt` format text files. See "File Format" section for details.

### Q: How many files can be analyzed at once?
A: The program supports batch analysis of all text files in a folder. The limit depends on your computer's performance.

## Project Structure

```
.
├── sequence_repeat_analyzer_en.py  # Main program (English version)
├── requirements.txt               # Python dependencies
└── README_EN.md                   # Project documentation (this file)
```

## Development Information

### Development Language
- Python 3.x

### Main Libraries
- tkinter: Graphical user interface
- matplotlib: Data visualization
- collections: Data structures

### Code Standards
- Follows PEP 8 coding standards
- Uses English comments and docstrings

## Changelog

### Version 1.0.0
- Initial release
- Basic sequence repeat recognition and statistics
- Graphical interface and visualization charts
- English interface
- Chart scrolling and navigation features

## License

This project is for educational and research purposes only.

## Contact

For questions or suggestions, please contact the developer.

---

**Note**: Please ensure input files follow the specified format for accurate analysis results.
