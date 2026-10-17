# GC / ENC Codon Usage Bias Analyzer

A GUI tool for analyzing codon usage bias in chloroplast (or other) genome CDS sequences. Provide an input folder and the tool batch-computes GC content, ENC and related metrics for every gene of every species, then summarizes everything into three tables.

## Overview

- **Input**: a folder in which each FASTA file is treated as one **species** (species name = file name), and each sequence record inside a file is treated as one **gene**.
- **Automatic gene-name recognition**: non-gene-name parts of the header (coordinates, strand direction, `join`, etc.) are removed automatically. For example:
  - `>psbD_[33312:34373](+)|NC01.fasta|NC01.fasta` → gene name `psbD`
  - `>rps16_join{[5754:5793](-),[4686:4894](-)}|...` → gene name `rps16`
- **Background analysis**: the analysis runs in a background thread, so the window shows live logs and a progress bar without freezing.

## Computed metrics

| Metric | Description |
|---|---|
| GC1 (%) / GC2 (%) / GC3 (%) | GC content at codon positions 1 / 2 / 3 |
| GC3s (%) | GC content at the 3rd position of **synonymous codons** (all degenerate codon families; Met, Trp and stop codons excluded) — identical to CodonW's GC3s |
| GC12 (%) | Average GC content of codon positions 1 and 2 |
| GC total (%) | Overall GC content of the sequence |
| ENC | Effective number of codons, computed exactly like **CodonW 1.4.4** (see notes below) |
| GC12–GC3 slope | Least-squares regression slope of GC12 on GC3 |
| R² | Coefficient of determination of that regression |

Notes:

- In the **species summary table**, slope / R² come from the regression of GC12 on GC3 across the genes of that species.
- In the **gene summary table**, slope / R² come from the regression of GC12 on GC3 across species for that gene.
- At least 3 points are required for a regression; otherwise these two columns are left empty.

### CodonW compatibility

All metrics are computed with the exact conventions of CodonW 1.4.4
(verified gene-by-gene against its output, 77/77 genes identical):

- A **terminal stop codon** (TAA/TAG/TGA) is excluded from all statistics;
  internal stop codons, if present, are kept.
- **GC3s** is calculated over all synonymous (degenerate) codons, not only
  fourfold-degenerate codon families.
- **ENC** follows the original CodonW implementation (`enc()` in
  `codon_idx.c`) rather than a plain textbook Wright (1990) formula:
  - amino-acid families observed only once (n ≤ 1) are excluded from the
    per-degeneracy-class homozygosity average;
  - families whose homozygosity is ≈ 0 (perfectly even codon usage) are also
    excluded;
  - if the 3-fold class (Ile) has no valid family, the mean of the 2-fold
    and 4-fold class averages is used;
  - if any other degeneracy class has no valid family, ENC is left empty
    (CodonW prints `*****` in that case);
  - the result is clamped only at the upper bound of 61.

## Output

Results are saved in the `results` subfolder of the selected folder:

| File | Content |
|---|---|
| `1_all_species_all_genes_detail.csv` | One row per species–gene combination with the raw computed values |
| `2_species_averages.csv` | One row per species; each metric is the **mean ± standard deviation** over all genes of that species |
| `3_gene_averages.csv` | One row per gene; each metric is the **mean ± standard deviation** over all species for that gene |

Example of the summary format:

```
Species    GC1 (%)      GC2 (%)      GC3 (%)      GC3s (%)     GC12 (%)     ENC          GC12–GC3 slope  R²
NC01       46.0 ± 6.6   39.9 ± 6.3   28.9 ± 4.8   27.3 ± 7.7   43.0 ± 5.0   47.5 ± 7.9   0.005           0.000
```

The CSV files use UTF-8 with BOM encoding and open directly in Excel without garbled characters.
If `openpyxl` is installed, the tool additionally writes a combined Excel workbook `GC_ENC_summary.xlsx` (three worksheets).

## Usage

### Requirements

- Python 3.8 or newer (standard library only; `tkinter` ships with Python)
- Optional: `openpyxl` (for xlsx output) — install with: `pip install openpyxl`

### Running

```bash
python GC_ENC_Analyzer.py
```

or simply double-click the file (provided `.py` is associated with Python).

Steps:

1. Click **Browse...** and select the folder containing the per-species FASTA files (e.g. `2各自共有`);
2. Click **Start analysis** and wait until the log says "Done";
3. Click **Open results folder**, or navigate manually to the `results` subfolder of the selected folder.

> Tip: to display Latin binomials (e.g. `A. officinalis`) in the tables, simply rename the FASTA files to the corresponding species names — the program uses the file name as the species name.

## Input data format

Standard FASTA format, for example:

```
>psbD_[33312:34373](+)|NC01.fasta|NC01.fasta
ATGACTATAGCTATTGGTAAATTTGCCAAAGAAGAAAATGATTTATTTGATATTATGGATGAT...
>rps16_join{[5754:5793](-),[4686:4894](-)}|NC01.fasta|NC01.fasta
ATGGTAAAACTTCGTTTGAAACGATGTGGTAGAAAGCAACGAGCGATCTATCGAATCGTTGCA...
```

Requirements:

- Sequences should be CDS (ideally starting with a start codon and having a length that is a multiple of 3; any trailing incomplete codon is removed automatically);
- Non-ACGT characters are ignored automatically;
- Records with fewer than 3 codons are skipped with a warning in the log;
- Supported extensions: `.fasta`, `.fa`, `.faa`, `.fna`, `.fst`.

## Files

```
GC_ENC_Analyzer/
├── GC_ENC_Analyzer.py   # Main program (GUI + all computation logic)
└── README.md            # This documentation
```

## Reference

- Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23–29.
