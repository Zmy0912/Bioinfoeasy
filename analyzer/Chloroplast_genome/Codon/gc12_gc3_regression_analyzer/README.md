# GC12-GC3 Regression Analyzer

A Python GUI application for analyzing GC12 and GC3 (GC content at first+second and third codon positions) regression across multiple species.

## Features

- **Batch Analysis**: Analyze CSV files from multiple species in a single folder
- **Linear Regression**: Calculate regression coefficients using scipy.stats.linregress
- **Comprehensive Statistics**: Slope, intercept, R², P-value, standard error, Pearson correlation, and mean values
- **Results Export**: Save results to CSV format with regression equation
- **User-Friendly Interface**: Built with tkinter for easy interaction

## Requirements

- Python 3.7+
- NumPy
- SciPy

Install dependencies:

```bash
pip install numpy scipy
```

## Usage

1. Run the program:

```bash
python gc12_gc3_regression_analyzer.py
```

2. Click "Browse..." to select a folder containing `_GC_analysis.csv` files
3. Click "Analyze" to process all files
4. View results in the table
5. Click "Save Results" to export to CSV

## Input File Format

The program expects CSV files with tab-separated or comma-separated columns. Required columns:
- `GC3 (%)` - GC content at third codon position
- `GC12 (%)` - GC content at first and second codon positions

Example header:
```
Source File,Gene ID,Gene Name,Sequence Length,GC1 (%),GC2 (%),GC3 (%),GC12 (%),Overall GC (%)
```

## Output

### Regression Statistics
For each species file:
- **Species**: Species/gene file name
- **Gene Count**: Number of genes analyzed
- **Slope**: Regression coefficient (GC12 = slope × GC3 + intercept)
- **Intercept**: Regression intercept
- **R²**: Coefficient of determination
- **P-value**: Statistical significance of regression
- **Std Error**: Standard error of the slope
- **Pearson r**: Pearson correlation coefficient
- **Mean GC12**: Average GC12 value
- **Mean GC3**: Average GC3 value
- **Regression Equation**: Full equation string

### CSV Export
Comma-separated values with all regression statistics per species.

## Regression Model

The program performs linear regression with GC3 as the independent variable (X) and GC12 as the dependent variable (Y):

```
GC12 = slope × GC3 + intercept
```

This matches the methodology used in CodonO and similar codon usage analysis tools.

## Related Tools

- **gc_analysis_plotter.py**: Single-file GC12-GC3 plotting tool
- **multi_species_analyzer.py**: ENC-GC3 analysis for multiple species
- **gc3_enc_plotter.py**: ENC-GC3 scatter plot visualization

## License

MIT License
