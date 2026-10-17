<div align="center">

# Bioinfoeasy

**A Software Package for Simplifying Bioinformatics Analysis**

[![GitHub release](https://img.shields.io/github/v/release/Zmy0912/Bioinfoeasy)](https://github.com/Zmy0912/Bioinfoeasy/releases)
[![GitHub license](https://img.shields.io/github/license/Zmy0912/Bioinfoeasy)](https://github.com/Zmy0912/Bioinfoeasy/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.22712158.svg)](https://doi.org/10.5281/zenodo.22712158)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/Zmy0912/Bioinfoeasy/pulls)

</div>

> **Official update v0.1.3** | Bioinfoeasy: A software package that makes bioinformatics work easier.

---

## What's new in v0.1.3

This release is centred on **codon usage bias and repeat-sequence analysis**,
and most of it grew directly out of the peer-review process of an ongoing
study. The headline additions are:

- **KaKs Analyzer** — filters Ka/Ks results and imputes missing values
  (`.xlsx` and `.kaks` auto-detected).
- **ENC-GC3s Analyzer** — quantifies how many genes fall below / on / above
  Wright's neutral expectation curve.
- **GC / ENC Analyzer** — batch GC1/GC2/GC3/GC3s/GC12 and ENC (CodonW 1.4.4
  compatible) with GC12–GC3 regression.
- **COA plotters** — supplementary figures for CodonW correspondence
  analysis, with correlation statistics against GC3s / GC / GC12.
- **SSR Region Classifier** and **MISA SSR Analyzer** — SSR classification
  into CDS / Intron / IGS, plus publication-ready SSR figures.
- **Repeat Analyzer**, **GC12-GC3 Regression Table**,
  **FASTA Header Renamer** and **CodonW2RSCUPlot**.

See the [changelog](CHANGELOG.md) for the complete list.

---

## Acknowledgments

<p style="font-size: 28px; font-weight: bold;">山高水长，莫失莫忘</p>

> *(Though mountains are high and waters long, never forget each other)*

I will forever remember the guidance from the one who illuminated my path during those early years. My best wishes to you.

---

## Table of Contents

- [What's new in v0.1.3](#whats-new-in-v013)
- [Overview](#overview)
- [Features](#features)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Modules](#modules)
- [Usage](#usage)
- [Conventions](#conventions)
- [Citation](#citation)
- [Contact & Support](#contact--support)

## Overview

Bioinfoeasy provides a comprehensive collection of simple and practical
bioinformatics tools designed to simplify biological sequence data processing,
analysis, and visualization. This software package is developed to make
bioinformatics work easier for researchers and students.

Most tools are standalone Python scripts that run entirely on your local
machine — no data is uploaded to any server.

## Features

- **🌐 Multi-language Support** — Developed using Python, R, and other
  programming languages.
- **💻 Command Line & GUI Tools** — Both batch-friendly scripts and
  point-and-click graphical interfaces.
- **🌍 Bilingual Interface** — English (`EN`) and Chinese (`CN`) versions for
  most tools.
- **📊 Comprehensive Analysis** — Genome annotation, sequence processing,
  repeat analysis, codon usage bias, Ka/Ks selection pressure, SSR
  classification, and statistics.
- **🧬 Codon Usage Bias Suite** — ENC / GC3s / GC12–GC3 neutrality,
  correspondence analysis (COA), RSCU aggregation, and Wright's neutral
  expectation test.
- **🎨 Visualization Tools** — Heatmaps, Pi plots, codon usage plots, GC12-GC3 /
  GC3-ENC / PR2 / COA plots, SSR figures, and a time scale generator.
- **🔄 Format Conversion** — GenBank to mVISTA, GFF3 to BED, CodonW output to
  RSCU-Plot input, and various sequence format conversions.
- **🛠️ Utility Tools** — Batch renaming, FASTA header management, file listing,
  and text processing utilities.
- **🔬 Review-Ready Outputs** — Several tools were written to answer specific
  reviewer requests and export tables/figures straight into a manuscript.

## Repository Structure

```text
Bioinfoeasy/
├── .github/                                  # GitHub community files
│   ├── ISSUE_TEMPLATE/
│   │   ├── bug_report.md
│   │   ├── feature_request.md
│   │   └── config.yml
│   └── PULL_REQUEST_TEMPLATE.md
├── analyzer/                                 # Analysis tools
│   ├── Chloroplast_genome/
│   │   ├── Codon/                            # Codon usage bias analysis
│   │   │   ├── codonw_analyzer/
│   │   │   │   └── blk_analyzer/{CN,EN}      # CodonW (blk) statistics
│   │   │   ├── gc_analyzer/                  # GC content analysis
│   │   │   ├── gc12_gc3_regression_analyzer/ # GC12-GC3 regression
│   │   │   ├── gc3enc_analyzer/              # GC3-ENC multi-species analysis
│   │   │   └── rscu_kruskal_wallis/{CN,EN}   # RSCU Kruskal-Wallis test
│   │   ├── cp_genome_comparison/{CN,EN}      # Assembly comparison
│   │   └── quadrants_analyzer/               # Four-region (quadrant) analysis
│   ├── ENC_GC3s_Analyzer/                    # Wright neutral expectation test
│   ├── GC_ENC_Analyzer/                      # GC1/2/3, GC3s, GC12 and ENC
│   ├── KaKs_Analyzer/                        # Ka/Ks filtering & imputation
│   ├── RepeatAnalyzer/                       # Direct & palindromic repeats
│   ├── SSR_Region_Classifier/                # SSR → CDS / Intron / IGS
│   ├── misa_tool/                            # MISA SSR toolkit
│   │   ├── misa_parser/                      # MISA result parser
│   │   └── misa_ssr_viewer/                  # SSR figures (heatmap/pie/bar)
│   └── sequence_repeat_analyzer/             # Repeat sequence analysis
├── Plots/                                    # Visualization tools
│   ├── Codon_plot/
│   │   ├── draw_codon_usage/{v0.1.0,v0.1.1}/{CN,EN}
│   │   ├── GC12_GC3/{v1.0.0,v2.0.0,v3.0.0}
│   │   │   └── gc12_gc3_regression_table/    # Regression statistics table
│   │   ├── GC3_ENC/{v1.0.0,v2.0.0,v3.0.0}
│   │   ├── PR2_plot/{v1.0.0,v2.0.0,v3.0.0,v4.0.0,v5.0.0}
│   │   ├── coa_plotter/                      # COA supplementary figures
│   │   └── coa_simple_plotter/               # COA main-panel-only variant
│   ├── Heatmap/                              # Heatmap visualization (CLI)
│   ├── heatmap_gui/                          # Heatmap tool (GUI)
│   ├── pi_plot/
│   │   ├── pi_gene_analyzer/                 # Pi gene analyzer
│   │   └── Plot/{v1.0.0,v1.1.0,v2.0.0}       # Pi plot + GUI (v2.0.0)
│   └── time_scale_generator/{CN,EN}          # Phylogenetic time scale generator
├── Seqinfo/                                  # Genome annotation analysis
│   ├── v1.0.0/
│   └── v3.0.0/
├── Seqtool/                                  # Sequence processing tools
│   ├── bed_converter/
│   │   ├── gene_to_0/{CN,EN}                 # Gene → 0-based BED
│   │   └── gff3_to_bed_converter/{CN,EN}     # GFF3 → BED (GUI)
│   ├── comparison/                           # FASTA classification & comparison
│   ├── fasta_header_renamer/                 # FASTA header ⇄ filename renaming
│   ├── fasta_joiner/{v1.0.0,v2.0.0}          # Merge FASTA files
│   ├── fasta_validator/                      # FASTA validation & fixing
│   ├── genbank_to_mvista_converter/{v1.0.0,v1.1.0}
│   └── other_tools/{CN,EN}                   # Misc. sequence utilities
├── other_tools/                              # General utility tools
│   ├── batch_renamer_tool/{CN,EN}            # Batch file renaming
│   ├── codonw2rscuplot/                      # CodonW RSCU → RSCU-Plot input
│   └── file_lister/                          # Folder file listing (HTML)
├── archive/                                  # Superseded tool versions (see its README)
├── .editorconfig
├── .gitattributes
├── .gitignore
├── CHANGELOG.md                              # Version history
├── CITATION.cff                              # Citation metadata
├── CODE_OF_CONDUCT.md
├── CONTRIBUTING.md
├── LICENSE                                   # MIT
├── README.md                                 # This file
├── SECURITY.md
└── requirements.txt                          # Aggregated dependencies
```

## Installation

### Prerequisites

- Python 3.7 or higher (3.9+ recommended)
- `pip` package manager

### Install all common dependencies

```bash
git clone https://github.com/Zmy0912/Bioinfoeasy.git
cd Bioinfoeasy
pip install -r requirements.txt
```

### Install for a single tool

Most tools ship their own `requirements.txt` with the minimal dependency set:

```bash
cd <tool_directory>
pip install -r requirements.txt
```

## Modules

### 1. Analyzer (`analyzer/`)

| Sub-module | Description |
|------------|-------------|
| `Chloroplast_genome/Codon/codonw_analyzer` | CodonW (blk) codon usage bias statistics (CN/EN) |
| `Chloroplast_genome/Codon/gc_analyzer` | GC content analysis |
| `Chloroplast_genome/Codon/gc12_gc3_regression_analyzer` | GC12-GC3 neutral/selection regression |
| `Chloroplast_genome/Codon/gc3enc_analyzer` | Multi-species ENC-GC3 analysis |
| `Chloroplast_genome/Codon/rscu_kruskal_wallis` | RSCU Kruskal-Wallis statistical test (CN/EN) |
| `Chloroplast_genome/cp_genome_comparison` | Compare two chloroplast genome assemblies (CN/EN) |
| `Chloroplast_genome/quadrants_analyzer` | Four-region (quadrant) chloroplast analysis |
| `ENC_GC3s_Analyzer` | Quantifies genes below / on / above Wright's ENC-GC3s neutral curve |
| `GC_ENC_Analyzer` | Batch GC1/GC2/GC3/GC3s/GC12 and ENC (CodonW 1.4.4 compatible) with GC12–GC3 regression |
| `KaKs_Analyzer` | Ka/Ks result filtering, grouping and missing-value imputation (`.xlsx` / `.kaks`) |
| `RepeatAnalyzer` | Direct and palindromic repeat statistics from Vmatch `repfind` output |
| `SSR_Region_Classifier` | Classify every SSR into CDS / Intron / IGS using MISA + GenBank files |
| `misa_tool` | MISA SSR toolkit — result parser (`misa_parser`) and figure viewer (`misa_ssr_viewer`) |
| `sequence_repeat_analyzer` | Repeat sequence analysis |

### 2. Plots (`Plots/`)

| Sub-module | Description |
|------------|-------------|
| `Codon_plot/draw_codon_usage` | Codon usage bias plots (v0.1.0, v0.1.1; CN/EN) |
| `Codon_plot/GC12_GC3` | GC12 vs GC3 neutrality plots (v1.0.0 – v3.0.0) |
| `Codon_plot/GC12_GC3/gc12_gc3_regression_table` | GC12–GC3 regression statistics as a sortable table (batch, CSV export) |
| `Codon_plot/GC3_ENC` | GC3 vs ENC plots (v1.0.0 – v3.0.0) |
| `Codon_plot/PR2_plot` | PR2 (parity rule 2) plots (v1.0.0 – v5.0.0) |
| `Codon_plot/coa_plotter` | CodonW correspondence analysis (COA) supplementary figures + correlation tests |
| `Codon_plot/coa_simple_plotter` | Streamlined COA variant — main axis panel only, batch/multi-species |
| `Heatmap` | Heatmap visualization (CLI) |
| `heatmap_gui` | Heatmap visualization (GUI) |
| `pi_plot/pi_gene_analyzer` | Gene nucleotide polymorphism (Pi) analyzer |
| `pi_plot/Plot` | Pi sliding-window gene annotation plots, incl. GUI (v2.0.0) |
| `time_scale_generator` | Phylogenetic time scale generator (CN/EN) |

### 3. Seqtool (`Seqtool/`)

| Sub-module | Description |
|------------|-------------|
| `bed_converter/gene_to_0` | Gene to 0-based BED coordinate converter (CN/EN) |
| `bed_converter/gff3_to_bed_converter` | GFF3 to BED converter with GUI (CN/EN) |
| `comparison` | FASTA sample classification and sequence comparison |
| `fasta_header_renamer` | Batch FASTA header ⇄ filename renaming with a multi-sequence editor (CN/EN) |
| `fasta_joiner` | Merge multiple FASTA files (v1.0.0, v2.0.0) |
| `fasta_validator` | FASTA format validation and fixing |
| `genbank_to_mvista_converter` | GenBank/GFF3 to mVISTA converter (v1.0.0, v1.1.0) |
| `other_tools` | Sequence extraction, renaming, removal, complement, cleaner GUI, etc. (CN/EN) |

### 4. Seqinfo (`Seqinfo/`)

**Genome Annotation Analyzer** (`v1.0.0`, `v3.0.0`) — extracts comprehensive
information from genome annotation files:

- Basic information (ID, species, length)
- Gene statistics (CDS, rRNA, tRNA, ncRNA)
- GC content and topology
- Taxonomy and metadata

Supported formats: GenBank (`.gb`, `.gbk`), GFF/GFF3, FASTA.

### 5. Other Tools (`other_tools/`)

| Sub-module | Description |
|------------|-------------|
| `batch_renamer_tool` | Batch file renaming with pattern matching (CN/EN) |
| `codonw2rscuplot` | Aggregate per-species CodonW/MEGA RSCU output into RSCU-Plot input |
| `file_lister` | Browser-based folder file listing with sort/search/CSV export (CN/EN) |

## Usage

Each module has its own `README.md` with detailed instructions. Common examples:

```bash
# Analyze a genome annotation file
cd Seqinfo/v3.0.0
python genome_annotation_analyzer.py chloroplast.gb

# Generate a codon usage plot (English UI)
cd Plots/Codon_plot/draw_codon_usage/v0.1.1/EN
python draw_codon_usage_en.py -i input.fasta -o output.pdf

# Convert GFF3 to BED (GUI)
cd Seqtool/bed_converter/gff3_to_bed_converter/EN
python gff3_to_bed_gui.py

# Generate a time scale
cd Plots/time_scale_generator/EN
python time_scale_generator.py

# Batch rename files
cd other_tools/batch_renamer_tool/EN
python batch_renamer.py

# Filter and impute Ka/Ks results (GUI, .xlsx or .kaks)
cd analyzer/KaKs_Analyzer
python kaks_analyzer.py

# Test genes against Wright's ENC-GC3s neutral expectation curve (GUI)
cd analyzer/ENC_GC3s_Analyzer
python ENC_GC3s_Analyzer.py

# Batch GC / ENC codon usage bias statistics for many species (GUI)
cd analyzer/GC_ENC_Analyzer
python GC_ENC_Analyzer.py

# COA supplementary figures from CodonW output (GUI)
cd Plots/Codon_plot/coa_plotter
python coa_plot_plotter.py

# Classify SSRs into CDS / Intron / IGS (GUI)
cd analyzer/SSR_Region_Classifier
python ssr_region_classifier.py
```

Newer modules place their Chinese documentation under a `中文/` sub-directory
rather than `CN/`; both are supported and the tool scripts are identical.

For the browser-based file lister, open
`other_tools/file_lister/file_lister.html` (English) or
`other_tools/file_lister/文件列表查看器.html` (Chinese) directly in a browser.

## Conventions

To keep the repository consistent, all modules follow these rules:

- **Language directories** are named `EN` (English) and `CN` (Chinese).
- **Version directories** use the form `vMAJOR.MINOR.PATCH` (e.g. `v1.0.0`).
- **Tool directory names** use lowercase `snake_case`.
- **Every module** provides a `README.md` describing what it does and how to run it.
- Superseded versions are moved to `archive/` instead of being kept inline.

## Citation

If you use Bioinfoeasy in your research, please cite:

```bibtex
@software{zhang_mingyuan_2026_bioinfoeasy,
  author    = {Zhang, Mingyuan},
  title     = {Bioinfoeasy: A Software Package for Simplifying Bioinformatics Analysis},
  month     = sep,
  year      = 2026,
  publisher = {Zenodo},
  version   = {v0.1.3},
  doi       = {10.5281/zenodo.22712158},
  url       = {https://doi.org/10.5281/zenodo.22712158}
}
```

Version-specific DOIs are listed on the Zenodo record. See also the
[`CITATION.cff`](CITATION.cff) file in this repository.

## Contact & Support

- **Author**: Mingyuan Zhang (Zmy0912)
- **Email**: myzhang0726@foxmail.com
- **Repository**: [https://github.com/Zmy0912/Bioinfoeasy](https://github.com/Zmy0912/Bioinfoeasy)
- **Issues**: [Submit an issue on GitHub](https://github.com/Zmy0912/Bioinfoeasy/issues)

For any questions or suggestions, please contact the author or raise them in Issues.

---

**Best wishes,**

**Zmy0912**

**October 17, 2026**
