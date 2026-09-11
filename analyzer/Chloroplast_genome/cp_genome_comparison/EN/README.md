# Chloroplast Genome Assembly Comparison Tool

## Introduction

This tool compares two assembled chloroplast genome sequences against a reference sequence to automatically determine **which assembly is more accurate**. It employs multi-dimensional scoring and outputs intuitive comparison results and visualizations.

**Key Feature**: A novel fragment position mapping method that detects **quadripartite structural errors** — including IR swap, local inversions, and junction order errors — which are difficult to identify using traditional global alignment and k-mer analysis.

## System Requirements

- **Python** 3.7 or later
- **Biopython** 1.76 or later

## Installation

```bash
pip install biopython
```

## Quick Start

1. Double-click `cp_genome_comparison.py` to launch the program
2. Select three FASTA format files:
   - **Reference Sequence**: Known chloroplast genome reference sequence
   - **Assembly 1**: Result from the first assembly
   - **Assembly 2**: Result from the second assembly
3. Click **「Start Analysis」**
4. View results: The program automatically determines which assembly is more accurate

## Interface Overview

The program has three tabs:

| Tab | Description |
|--------|------|
| 📊 Summary | Score comparison, basic statistics, K-mer analysis, **quadripartite structure analysis summary**, and final conclusion |
| 📋 Details | Detailed metrics for each dimension (length, GC, alignment, coverage) and **structural analysis details** |
| 📈 Identity Plot | Sliding window identity curves, visually showing similarity across different regions |

## Scoring System

Total Score = Base Score (max 100) + Structural Integrity Bonus (up to +5):

| Dimension | Weight | Description |
|------|:----:|------|
| Identity Score | 40% | Identity based on Biopython chunked global alignment |
| Coverage Score | 30% | K-mer coverage of the reference by the assembly |
| Length Similarity | 15% | Proximity of assembly length to reference length |
| GC Similarity | 15% | GC content difference |
| **Structural Integrity** | **+5 pts** | **Quadripartite structure correctness (IR swap, inversion, junction error detection)** |

> Structural integrity scoring rules: No errors +5.0; Global inversion (not an error) +4.5; Rearrangement/inversion penalty by severity; IR swap -2.0

## Analysis Principles

### Core Algorithms

1. **K-mer Similarity Analysis** (k=31)  
   Decompose sequences into short strings of length 31, comparing the Jaccard similarity coefficient and shared proportion of k-mer sets. This method is insensitive to sequence rearrangements.

2. **Chunked Global Alignment**  
   For large genomes (>50 kb), the program uses a chunking strategy: the reference is divided into 5000 bp chunks; each chunk is located in the assembly via k-mer matching to find the best match region, followed by a localized global alignment; results are then aggregated.

3. **Sliding Window Analysis**  
   Scans the reference with a 5000 bp window and 2000 bp step size, computing identity for each window via k-mer matching, generating an identity curve along the genome.

4. **Coverage Calculation**  
   Checks the proportion of reference k-mers present in the assembly, evaluating assembly completeness.

5. **Quadripartite Structure Analysis** (New in v2.0)  
   Divides the reference into 20 equal segments, using forward/reverse-complement k-mer dual matching to locate each segment's position and orientation in the assembly. Analyzes segment arrangement patterns:
   - **Global Inversion** (≥85% segments reversed): Assembly is the reverse-complement of the reference, not a structural error
   - **Position Breaks**: Segments not arranged in reference order, indicating rearrangement or IR swap
   - **Local Inversion**: Consecutive segments match in reverse-complement orientation
   - **Fold-back Pattern**: Hallmark of IR swap — segments in the assembly first increase then fold back

   Additionally uses 8 equally-spaced checkpoint k-mer fingerprints to verify sequence order and confirm error types.

### Structural Error Types

| Error Type | Description | Detection Method |
|----------|------|----------|
| IR Swap | IRa/IRb regions are swapped in the assembly | Position fold-back + multiple breaks + local inversions |
| Local Inversion | Part of the assembly is in reverse-complement orientation | Consecutive segment orientation reversal |
| Junction Error | Incorrect LSC/SSC/IR boundary connection order | Checkpoint k-mer order verification |
| Global Inversion | Entire assembly is the reverse-complement of the reference | ≥85% segments reversed (not an error) |

### Output Results

- Total scores and rankings for assemblies (including structural scores)
- Detailed comparison statistics
- Sliding window identity distribution plot
- Quadripartite structure analysis details
- Exportable text report

## File Formats

Supported FASTA extensions: `.fasta` `.fa` `.fna` `.fas` `.txt`

The program automatically handles encoding issues (UTF-8 / GBK) and supports FASTA files with Chinese description lines.

## Notes

- Chloroplast genomes are typically 120–180 kb; the program automatically uses memory-optimized chunked alignment for large genomes
- The alignment process may take several minutes; please be patient
- If scores are close, examine the sliding window identity plot and structural analysis results for large-scale rearrangements, inversions, or IR swaps
- Quadripartite structure analysis is only enabled for sequences longer than 30 kb

## References

1. Daniel H, Lin CS, Yu M, Chang WJ. Chloroplast genomes: diversity, evolution, and applications in genetic engineering. *Genome Biology*, 2016, 17: 134.
2. Jansen RK, Ruhlman TA. Plastid Genomes of Seed Plants. In: *Genomics of Chloroplasts and Mitochondria*. Springer, 2012: 103–126.
3. Ondov BD, Treangen TJ, et al. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology*, 2016, 17: 132.
4. Gurevich A, Saveliev V, et al. QUAST: quality assessment tool for genome assemblies. *Bioinformatics*, 2013, 29(8): 1072–1075.
5. Cock PJA, Antao T, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 2009, 25(11): 1422–1423.

## License

This tool is freely available for research purposes. If used in publications, please cite the references above.
