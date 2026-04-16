# RSCU Kruskal-Wallis Analysis Tool

## Overview

This software performs Kruskal-Wallis non-parametric tests on RSCU (Relative Synonymous Codon Usage) data to analyze differences in codon usage preferences between species groups.

## Features

- **Data Loading**: Support for xlsx format RSCU data files
- **Species Grouping**: Auto-grouping or manual group definition
- **Flexible Analysis**: All codons or filter by amino acids
- **Kruskal-Wallis Test**: Compare codon usage across species groups
- **Result Interpretation**: Detailed statistical results and biological insights
- **Visualization**: Export boxplots for significant differences
- **Export Options**: Full results and significant results (Excel/CSV)

## Installation

```bash
pip install -r requirements.txt
```

Dependencies:
- pandas >= 1.3.0
- numpy >= 1.20.0
- scipy >= 1.7.0
- openpyxl >= 3.0.0
- matplotlib >= 3.5.0

## Usage

### 1. Run the Program

```bash
python rscu_kruskal_wallis.py
```

### 2. Load Data

Click "Browse" to select your RSCU data file (xlsx format), then click "Load Data".

### 3. Set Up Groups

#### Auto Grouping
The software automatically groups species by name prefix. Suitable for data with clear grouping characteristics.

#### Manual Grouping
For custom grouping (e.g., geographic distribution, phylogenetic relationships), prepare a group file:

**Group File Format** (CSV or Excel, two columns):

| Species | Group |
|---------|-------|
| Species_A | Group1 |
| Species_B | Group1 |
| Species_C | Group2 |
| ... | ... |

### 4. Analysis Options

- **Analysis Scope**: All codons or filter by amino acids
- **Significance Level**: Default is 0.05, adjustable

### 5. Run Test

Click "Run Kruskal-Wallis Test" to start analysis.

### 6. View Results

Results display:
- H statistic and p-value for each codon
- Significance markers: `***` (p<0.001), `**` (p<0.01), `*` (p<0.05), `ns` (not significant)
- Comparison of group means and detailed interpretation

### 7. Export Results

- **Export Full Results**: Save all codon test results and group statistics
- **Export Significant**: Save only significantly different codons
- **Plot Boxplots**: Generate boxplots for significantly different codons

## Interpreting Results

### About Kruskal-Wallis Test

The Kruskal-Wallis test is a non-parametric method for comparing independent samples from three or more groups. It is an extension of the Mann-Whitney U test and does not assume normal distribution.

### Significance Levels

| Significance | p-value | Biological Meaning |
|--------------|---------|-------------------|
| Highly significant | p < 0.001 | Very obvious difference in codon usage |
| Very significant | 0.001 ≤ p < 0.01 | Significant difference between groups |
| Significant | 0.01 ≤ p < 0.05 | Some difference between groups |
| Not significant | p ≥ 0.05 | No statistically significant difference |

### Application Scenarios

1. **Species Adaptation Studies**: Compare codon usage preferences across ecologically different species
2. **Phylogenetic Analysis**: Verify codon usage differences across evolutionary branches
3. **Gene Expression Studies**: Analyze codon preferences in high vs. low expression genes
4. **Comparative Genomics**: Study translational selection pressure across species

## Data Format Requirements

Input should be an Excel file with columns:

- **First column**: Species names (recommended column name: "Species")
- **Subsequent columns**: RSCU values for each codon, format: `Codon(AminoAcid)`, e.g., `GCA(Ala)`

Example data format:

| Species | GCA(Ala) | GCC(Ala) | GCG(Ala) | GCU(Ala) | ... |
|---------|----------|----------|----------|----------|-----|
| Species_A | 1.19 | 0.62 | 0.40 | 1.79 | ... |
| Species_B | 1.18 | 0.63 | 0.41 | 1.78 | ... |

## Notes

1. Minimum 3 samples per group recommended for reliable results
2. Grouping should have biological significance; avoid random grouping
3. Significant results should be interpreted in biological context
4. Bonferroni correction can be used as reference for multiple comparisons

## File Structure

```
260416软件/
├── rscu_kruskal_wallis.py    # Main program
├── requirements.txt           # Dependencies
├── example_groups.csv        # Example group file
└── README.md                 # This file
```

## References

1. Kruskal, W.H. & Wallis, W.A. (1952). Use of Ranks in One-Criterion Variance Analysis. Journal of the American Statistical Association.
2. Wright, F. (1990). The 'effective number of codons' used in a gene. Gene.

---

For questions or suggestions, please contact the developer.
