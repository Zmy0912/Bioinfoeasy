# Bioinfoeasy: A Software Package for Simplifying Bioinformatics Analysis

[![GitHub release](https://img.shields.io/github/v/release/Zmy0912/Bioinfoeasy)](https://github.com/Zmy0912/Bioinfoeasy/releases)
[![GitHub license](https://img.shields.io/github/license/Zmy0912/Bioinfoeasy)](https://github.com/Zmy0912/Bioinfoeasy/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19305267.svg)](https://doi.org/10.5281/zenodo.19305267)

> **Official update v0.1.1** | Bioinfoeasy: A software package that makes bioinformatics work easier.


---

## Acknowledgments

<p style="font-size: 28px; font-weight: bold;">山高水长，莫失莫忘</p>

> *(Though mountains are high and waters long, never forget each other)*

I will forever remember the guidance from the one who illuminated my path during those early years. My best wishes to you.

---

## Welcome everyone to use it and cite my Zenodo link: [https://doi.org/10.5281/zenodo.19305267](https://doi.org/10.5281/zenodo.19305267)

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Modules](#modules)
- [Usage](#usage)
- [Citation](#citation)
- [Contact & Support](#contact--support)

## Overview

Bioinfoeasy provides a comprehensive collection of simple and practical bioinformatics tools designed to simplify biological sequence data processing, analysis, and visualization. This software package is developed to make bioinformatics work easier for researchers and students.

## Features

- **🌐 Multi-language Support**: Developed using Python, R, and other programming languages
- **💻 Command Line Tools**: All tools provided as command-line utilities for batch processing and automation
- **🌍 Bilingual Interface**: English and Chinese language support for most tools
- **📊 Comprehensive Analysis**: Genome annotation, sequence processing, repeat analysis, codon usage, and visualization
- **🎨 Visualization Tools**: Heatmaps, Pi plots, codon usage plots, and time scale generators
- **🔄 Format Conversion**: GenBank to mVISTA, GFF3 to BED, and various sequence format conversions
- **🛠️ Utility Tools**: Batch renaming, file listing, and text processing utilities

## Project Structure

```
Bioinfoeasy/
├── .github/                              # GitHub community files
│   └── ISSUE_TEMPLATE/                   # Issue templates
│       ├── bug_report.md
│       └── feature_request.md
├── analyzer/                             # Analysis tools
│   ├── Chloroplast_genome/              # Chloroplast genome analysis
│   │   ├── Codon/                       # Codon usage analysis
│   │   │   ├── codonw_analyzer/         # CodonW-based analysis
│   │   │   ├── GC_ANALYZER/             # GC content analysis
│   │   │   ├── gc12_gc3_regression_analyzer/ # GC12-GC3 regression
│   │   │   ├── GC3ENCanalyzer/          # GC3-ENC multi-species analysis
│   │   │   └── rscu_kruskal_wallis/     # RSCU statistical analysis
│   │   └── Quadrants_analyzer/          # Quadrant analysis
│   ├── MISAtool/                        # SSR motif analysis
│   └── sequence_repeat_analyzer/        # Repeat sequence analysis
├── other_tools/                          # General utility tools
│   ├── batch_renamer_tool/              # Batch file renaming (CN/EN)
│   └── file_lister/                     # Folder file listing tool (CN/EN)
├── Plots/                               # Visualization tools
│   ├── Codon_plot/                      # Codon usage and GC analysis
│   │   ├── draw_codon_usage/            # Codon usage plotting (v0.1.0, v0.1.1)
│   │   ├── GC12_GC3/                    # GC12-GC3 analysis plotter
│   │   └── GC3_ENC/                     # GC3-ENC plotter
│   ├── Heatmap/                         # Heatmap visualization (CLI)
│   ├── HeatmapGUI/                      # GUI-based heatmap tool
│   ├── PiPlot/                          # Pi gene analysis plots
│   │   ├── pi_gene_analyzer/            # Pi gene analyzer
│   │   └── Plot/                        # Pi plot (v1.0.0, v1.1.0, v2.0.0)
│   └── Time Scale Generator/            # Time scale visualization (CN/EN)
├── Seqtool/                             # Sequence processing tools
│   ├── BED-Converter/                   # BED format conversion
│   │   ├── gene_to_0/                  # Gene to 0-based BED (CN/EN)
│   │   └── gff3_to_bed_converter/      # GFF3 to BED (CN/EN)
│   ├── comparision/                     # Sequence comparison tools
│   ├── fasta_joiner/                    # FASTA file joining
│   ├── fasta_validator/                 # FASTA file validation
│   ├── genbank_to_mvista_converter/     # GenBank to mVISTA (v1.0.0, v1.1.0)
│   ├── other_tools/                     # Miscellaneous sequence tools
│   └── Seqtool_CN/                      # Chinese version of Seqtool
├── Seqinfo/                             # Genome annotation analysis
│   ├── genome_annotation_analyzer.py
│   ├── README.md
│   └── requirements.txt
├── CITATION.cff                         # Citation information
├── CODE_OF_CONDUCT.md                   # Code of conduct
├── CONTRIBUTING.md                      # Contribution guidelines
├── LICENSE                              # License file
├── README.md                            # This file
└── SECURITY.md                          # Security policy
```

## Installation

### Prerequisites

- Python 3.7 or higher
- pip package manager

### Basic Installation

Most tools have individual `requirements.txt` files. Install dependencies for a specific tool:

```bash
cd <tool_directory>
pip install -r requirements.txt
```

### Common Dependencies

For most tools, you will need:

```bash
pip install biopython pandas matplotlib numpy openpyxl seaborn
```

## Modules

### 1. Analyzer (`analyzer/`)

**Chloroplast Genome Analysis**
- Codon usage analysis (CodonW, GC content, GC12-GC3 regression, GC3-ENC, RSCU statistics)
- Quadrant analysis of chloroplast genomes

**MISA Tool**
- SSR (Simple Sequence Repeat) motif identification and analysis

**Sequence Repeat Analyzer**
- Comprehensive repeat sequence analysis

### 2. Plots (`Plots/`)

**Codon Plot**
- Codon usage plots with multiple versions
- GC12 vs GC3 analysis
- GC3 vs ENC (Effective Number of Codons) plots
- Codon usage bias visualization

**Heatmap**
- Data clustering and visualization (CLI version)
- GUI version available for interactive use

**Pi Plot**
- Pi gene analysis visualization
- Multiple versions (v1.0.0, v1.1.0, v2.0.0) with GUI support

**Time Scale Generator**
- Phylogenetic time scale visualization
- Available in English and Chinese

### 3. Seqtool (`Seqtool/`)

**Format Converters**
- Gene to 0-based BED converter (CN/EN)
- GFF3 to BED conversion with GUI (CN/EN)
- GenBank to mVISTA format converter (v1.0.0, v1.1.0)

**Sequence Processing**
- FASTA validator for sequence integrity checking
- FASTA joiner for merging multiple files
- Sequence extraction, renaming, and removal
- Complement sequence generation
- Batch sequence processing
- FASTA cleaner with GUI

**Sequence Comparison**
- FASTA classifier and comparison tools

### 4. Seqinfo (`Seqinfo/`)

**Genome Annotation Analyzer**
Extract comprehensive information from genome annotation files:
- Basic information (ID, species, length)
- Gene statistics (CDS, rRNA, tRNA, ncRNA)
- GC content and topology
- Taxonomy and metadata

Supported formats: GenBank (.gb, .gbk), GFF/GFF3, FASTA

### 5. Other Tools (`other_tools/`)

**Batch Renamer Tool**
- Batch file renaming with pattern matching
- Available in English and Chinese

**File Lister**
- Browser-based folder file listing with table view
- Sort, search, copy, and export CSV
- Available in English and Chinese

## Usage

Each module has its own README file with detailed instructions. Here are some common examples:

### Analyze Genome Annotation

```bash
cd Seqinfo
python genome_annotation_analyzer.py chloroplast.gb
```

### Generate Codon Usage Plot

```bash
cd Plots/Codon_plot/draw_codon_usage/0.1.1
python draw_codon_usage_en.py -i input.fasta -o output.pdf
```

### Convert GFF3 to BED

```bash
cd Seqtool/BED-Converter/gff3_to_bed_converter/EN
python gff3_to_bed_gui.py
```

### Generate Time Scale

```bash
cd "Plots/Time Scale Generator/EN"
python time_scale_generator.py
```

### Batch Rename Files

```bash
cd other_tools/batch_renamer_tool/EN
python batch_renamer.py
```

### List Files in a Folder

Open `other_tools/file_lister/file_lister.html` in a browser, or double-click `文件列表查看器.html` for the Chinese version.

## Citation

If you use Bioinfoeasy in your research, please cite:

```bibtex
@software{zhang_mingyuan_2026_18545517,
  author       = {Zhang, Mingyuan},
  title        = {Bioinfoeasy: A Software Package for Simplifying Bioinformatics Analysis},
  month        = apr,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.18545517},
  url          = {https://doi.org/10.5281/zenodo.19305267}
}
```

Or use the CITATION.cff file available in the repository.

## Contact & Support

- **Author**: Mingyuan Zhang (Zmy0912)
- **Email**: myzhang0726@foxmail.com
- **Repository**: [https://github.com/Zmy0912/Bioinfoeasy](https://github.com/Zmy0912/Bioinfoeasy)
- **Issues**: [Submit an issue on GitHub](https://github.com/Zmy0912/Bioinfoeasy/issues)

For any questions or suggestions, please contact the author or raise them in Issues.

---

**Best wishes,**

**Zmy0912**

**April 1, 2026**
