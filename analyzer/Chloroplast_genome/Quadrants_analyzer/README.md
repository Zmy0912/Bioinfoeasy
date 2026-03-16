# Chloroplast Genome Four-Region Analyzer

A graphical tool for identifying and visualizing the four regions (LSC, IRA, SSC, IRB) of chloroplast genomes.

## Features

- Supports GenBank format file input
- Automatically identifies the four regions of chloroplast genomes with position and length
- Graphical interface for intuitive display of region information
- Provides visual diagram of genome structure
- Clean English interface

## Installation

Install the required dependencies:

```bash
pip install biopython matplotlib numpy
```

### Dependencies

- **biopython**: For parsing GenBank files
- **matplotlib**: For plotting visualization charts
- **numpy**: For numerical calculations
- **tkinter**: Built-in Python GUI library (usually pre-installed)

## Usage

1. Run the program:
```bash
python chloroplast_analyzer.py
```

2. Click the "Select GenBank File" button to load your data file

3. Click the "Analyze and Display" button to perform region analysis

4. View the detailed region information in the table on the left and the genome diagram on the right

## Region Description

Chloroplast genomes typically contain four main regions:

- **LSC (Large Single Copy)**: Large single-copy region, typically occupying 50-60% of the genome
- **IRA (Inverted Repeat A)**: Inverted repeat region A, identical sequence to IRB but in opposite orientation
- **SSC (Small Single Copy)**: Small single-copy region, typically occupying 10-15% of the genome
- **IRB (Inverted Repeat B)**: Inverted repeat region B, identical sequence to IRA but in opposite orientation

## Program Features

### Intelligent Region Identification

The program uses two strategies to identify regions:
1. **Marker gene-based identification**: Searches for typical chloroplast marker genes
2. **Ratio-based estimation**: Uses typical chloroplast genome structure ratios

### Visual Display

- Uses different colors to distinguish the four regions
- Displays start position, end position, and length for each region
- Intuitive linear genome diagram
- Supports axis labels and grid lines

## Notes

- Supported file formats: .gb, .gbk (GenBank format)
- The program automatically detects genome length
- If marker genes cannot be found, default ratios will be used to estimate regions
- Analysis results are for reference only; please verify with other methods

## Common Issues

1. **Program won't run?**
   - Ensure all dependency packages are installed
   - Check Python version (3.8+ recommended)

2. **Regions not identified?**
   - Confirm input file is in GenBank format
   - Check if file contains necessary gene annotation information
   - The program will use default ratios as a fallback option

3. **Display not clear?**
   - Adjust window size to change chart display
   - Program uses matplotlib for automatic window sizing

## License

This program is for academic research use only.
