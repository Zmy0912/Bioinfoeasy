# Chloroplast Genome Assembly Comparison Tool

Compares two assembled chloroplast genome sequences against a reference
sequence to automatically determine **which assembly is more accurate**, using
multi-dimensional scoring and visualizations.

It includes a novel fragment position mapping method that detects
**quadripartite structural errors** — including IR swap, local inversions, and
junction order errors — which are difficult to identify with traditional global
alignment or k-mer analysis.

## Contents

| Directory | Description | Files |
|-----------|-------------|-------|
| [`EN/`](EN/) | English version | `cp_genome_comparison.py`, `README.md`, `Explanation of the technical principles.md` |
| [`CN/`](CN/) | Chinese version | `cp_genome_comparison.py`, `README.md`, `技术原理说明.md` |

## Requirements

- Python 3.7 or later
- Biopython 1.76 or later

```bash
pip install biopython
```

## Usage

1. Launch the program:

   ```bash
   cd EN
   python cp_genome_comparison.py
   ```

2. Select the **reference** sequence and the two **assembly** FASTA files.
3. Click **Start Analysis** and read the result tabs (Summary / Details /
   Identity Plot).

See `EN/README.md` for scoring details and `Explanation of the technical
principles.md` for the method and references.

Back to [Chloroplast Genome Analysis](../README.md) · [project root](../../../README.md).
