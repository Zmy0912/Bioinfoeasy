# Bioinfoeasy: A Software Package for Simplifying Bioinformatics Analysis

[![GitHub release](https://img.shields.io/github/v/release/Zmy0912/Bioinfoeasy)](https://github.com/Zmy0912/Bioinfoeasy/releases)
[![GitHub license](https://img.shields.io/github/license/Zmy0912/Bioinfoeasy)](https://github.com/Zmy0912/Bioinfoeasy/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19305267.svg)](https://doi.org/10.5281/zenodo.19305267)

> **Official update v0.1.0** | Bioinfoeasy: A software package that makes bioinformatics work easier.
>
> **Great update!** You can now use the software to carry out a complete codon preference analysis process!

---

## Acknowledgments

<p style="font-size: 20px;">**山高水长，莫失莫忘**</p>

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
- **📊 Comprehensive Analysis**: Genome annotation, sequence processing, repeat analysis, and visualization
- **🎨 Visualization Tools**: Heatmaps, Pi plots, codon usage plots, and time scale generators
- **🔄 Format Conversion**: GenBank to mVISTA, GFF3 to BED, and various sequence format conversions

## Project Structure

```
Bioinfoeasy/
├── analyzer/                    # Analysis tools
│   ├── Chloroplast_genome/     # Chloroplast genome analysis
│   ├── MISAtool/               # SSR motif analysis
│   └── sequence_repeat_analyzer/ # Repeat sequence analysis
├── Plots/                      # Visualization tools
│   ├── Codon_plot/             # Codon usage and GC analysis
│   ├── Heatmap/                # Heatmap visualization
│   ├── HeatmapGUI/             # GUI-based heatmap tool
│   ├── PiPlot/                 # Pi gene analysis plots
│   └── Time Scale Generator/   # Time scale visualization
├── Seqtool/                    # Sequence processing tools
│   ├── BED-Converter/          # BED format conversion
│   ├── comparision/            # Sequence comparison tools
│   ├── fasta_joiner/           # FASTA file joining
│   ├── fasta_validator/        # FASTA file validation
│   ├── genbank_to_mvista_converter/ # GenBank to mVISTA conversion
│   ├── other_tools/            # Miscellaneous sequence tools
│   └── Seqtool_CN/             # Chinese version of Seqtool
├── Seqinfo/                    # Genome annotation analysis
│   └── genome_annotation_analyzer.py
├── txt_tools/                  # Text processing tools
│   └── batch_renamer_tool/     # Batch file renaming
└── CITATION.cff                # Citation information
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
- Quadrant analysis of chloroplast genomes
- GC content analysis
- Codon usage analysis with CodonW

**MISA Tool**
- SSR (Simple Sequence Repeat) motif identification and analysis

**Sequence Repeat Analyzer**
- Comprehensive repeat sequence analysis

### 2. Plots (`Plots/`)

**Codon Plot**
- GC12 vs GC3 analysis
- GC3 vs ENC (Effective Number of Codons) plots
- Codon usage bias visualization

**Heatmap**
- Data clustering and visualization
- GUI version available for interactive use

**Pi Plot**
- Pi gene analysis visualization
- Multiple versions with GUI support

**Time Scale Generator**
- Phylogenetic time scale visualization
- Available in English and Chinese

### 3. Seqtool (`Seqtool/`)

**Format Converters**
- BED format converter (gene to 0-based coordinates)
- GFF3 to BED conversion (with GUI)
- GenBank to mVISTA format converter

**Sequence Processing**
- FASTA validator for sequence integrity checking
- FASTA joiner for merging multiple files
- Sequence extraction and manipulation tools
- Sequence renaming and removal
- Complement sequence generation
- Batch sequence processing

**Sequence Comparison**
- FASTA classifier
- Sequence comparison tools

### 4. Seqinfo (`Seqinfo/`)

**Genome Annotation Analyzer**
Extract comprehensive information from genome annotation files:
- Basic information (ID, species, length)
- Gene statistics (CDS, rRNA, tRNA, ncRNA)
- GC content and topology
- Taxonomy and metadata

Supported formats: GenBank (.gb, .gbk), GFF/GFF3, FASTA

### 5. Text Tools (`txt_tools/`)

**Batch Renamer Tool**
- Batch file renaming with pattern matching
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
cd Plots/Codon_plot/draw_codon_usage/EN
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
