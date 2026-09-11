# Multi-Species ENC-GC3 Analyzer

A Python GUI application for analyzing ENC (Effective Number of Codons) and GC3 (GC content at third codon position) data from multiple species.

## Features

- **Batch Analysis**: Analyze .out files from multiple species in a single folder
- **Comprehensive Statistics**: Calculate mean, standard deviation, min, and max for both ENC and GC3 values
- **Summary Statistics**: Overall statistics across all species including Spearman correlation
- **Data Export**: Export results to CSV format for further analysis
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
python multi_species_analyzer.py
```

2. Click "Browse..." to select a folder containing `.out` files
3. Click "Analyze" to process all files
4. View results in the table and summary section
5. Optionally click "Export CSV" to save results

## Input File Format

The program expects `.out` files with tab-separated columns. The file should contain:
- Header row (automatically skipped)
- Gene data rows with ENC (column 9) and GC3 (column 10)

Example format:
```
title          T3s    C3s    A3s    G3s    CAI   CBI   Fop   Nc    GC3s   GC    L_sym  ...
gene_name      ...    ...    ...    ...    ...   ...   ...   44.90 0.264  ...    ...    ...
```

## Output

### Table Results
For each species file:
- Species name
- Gene count
- Mean ENC (with std, min, max)
- Mean GC3 (with std, min, max)

### Summary Section
- Total species and genes analyzed
- Overall ENC and GC3 statistics
- Spearman correlation between mean ENC and mean GC3

### CSV Export
Comma-separated values with all statistics per species.

## License

MIT License
