# GC Content Analyzer

A graphical tool for analyzing GC content at different codon positions (GC1, GC2, GC3) in gene sequences from FASTA files.

## Features

- 📊 Analyzes GC content at first (GC1), second (GC2), and third (GC3) codon positions
- 🔬 Calculates overall GC content and average GC12 (mean of GC1 and GC2)
- 🏷️ Processes multiple genes from a single FASTA file
- 📈 Generates CSV tables ready for statistical analysis and visualization
- 🎨 User-friendly graphical interface
- 🔍 Includes analysis report with statistics and warnings
- 📁 Output filename automatically linked to input filename

## Installation

Install required dependencies:

```bash
pip install biopython pandas
```

## Usage

### 1. Run the Program

```bash
python gc_analyzer.py
```

### 2. Interface Operation

- **Input FASTA File**: Select a FASTA file containing gene sequences
- **Output CSV File**: Choose where to save the analysis results (auto-generated based on input file name)
- Click "Analyze" button to execute the analysis

### 3. Output File Format

The output CSV file contains the following columns:

| Column | Description |
|--------|-------------|
| Source File | Name of the source FASTA file |
| Gene ID | Full gene identifier from FASTA header |
| Gene Name | Extracted gene name (cleaned identifier) |
| Sequence Length | Length of the gene sequence |
| GC1 (%) | GC content at first codon position |
| GC2 (%) | GC content at second codon position |
| GC3 (%) | GC content at third codon position |
| GC12 (%) | Average of GC1 and GC2 |
| Overall GC (%) | Overall GC content of the entire sequence |

### 4. GC Content Calculation

**GC1**: GC content at the first nucleotide of each codon
- Nucleotides at positions: 0, 3, 6, 9, ...

**GC2**: GC content at the second nucleotide of each codon
- Nucleotides at positions: 1, 4, 7, 10, ...

**GC3**: GC content at the third nucleotide of each codon
- Nucleotides at positions: 2, 5, 8, 11, ...

**GC12**: Average of GC1 and GC2 content

**Overall GC**: GC content of the entire sequence

## Example

### Input FASTA File (NC01.fasta)

```
>rps3_[83935:84597](-)|NC01.fasta|NC01.fasta
ATGGGACAAAAAATAAATCCACTTGGTTTCAGACTTGGTACAACCCAAAACCATCATTCCTTTTGGTTCGC...

>rps11_[80174:80590](-)|NC01.fasta|NC01.fasta
ATGGCAAAACCTATACCGAGAATTGGTTCGCGTAGGAATGGACGTATTGGTTTACGTAAGAATGGACGTAG...
```

### Output CSV File (NC01_GC_analysis.csv)

```
Source File,Gene ID,Gene Name,Sequence Length,GC1 (%),GC2 (%),GC3 (%),GC12 (%),Overall GC (%)
NC01.fasta,rps3_[83935:84597](-)|NC01.fasta|NC01.fasta,rps3,636,45.28,42.45,62.26,43.87,50.00
NC01.fasta,rps11_[80174:80590](-)|NC01.fasta|NC01.fasta,rps11,417,48.44,48.68,65.23,48.56,54.20
...
```

### Analysis Report

The program also displays an analysis report in the interface, including:
- Total number of genes analyzed
- Average GC1, GC2, and GC3 content
- GC3/GC1 ratio (indicator of codon usage bias)
- Warnings for sequences with lengths not divisible by 3
- Data preview of the first 10 rows

## Output File Naming

The output CSV filename is automatically generated based on the input filename:
- Input: `NC01.fasta`
- Output: `NC01_GC_analysis.csv`

This ensures easy association between input and output files.

## Visualization Example

Using Python to visualize GC content distribution:

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read analysis results
df = pd.read_csv('NC01_GC_analysis.csv')

# Create boxplot of GC content at different positions
plt.figure(figsize=(10, 6))
gc_columns = ['GC1 (%)', 'GC2 (%)', 'GC3 (%)']
df_melted = df[gc_columns].melt(var_name='Position', value_name='GC Content')

sns.boxplot(data=df_melted, x='Position', y='GC Content')
plt.title('GC Content Distribution at Different Codon Positions')
plt.ylabel('GC Content (%)')
plt.xlabel('Codon Position')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('gc_content_boxplot.png', dpi=300)
plt.show()

# Create scatter plot of GC3 vs GC1
plt.figure(figsize=(10, 6))
plt.scatter(df['GC1 (%)'], df['GC3 (%)'], alpha=0.6)
plt.plot([0, 100], [0, 100], 'r--', label='GC3 = GC1')
plt.xlabel('GC1 (%)')
plt.ylabel('GC3 (%)')
plt.title('GC3 vs GC1 Content')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('gc3_vs_gc1_scatter.png', dpi=300)
plt.show()
```

## Notes

1. The tool automatically removes non-DNA characters (A, T, G, C, N) from sequences
2. GC3/GC1 ratio > 1.0 often indicates codon usage bias in the genome
3. The program warns if sequence lengths are not multiples of 3 (potential issues)
4. Overall GC content is included for reference and comparison
5. Gene names are extracted from the FASTA header by removing special characters

## Technical Support

If you encounter any issues or have suggestions, please check the analysis report displayed in the interface for warnings and ensure your FASTA file is properly formatted.
