# FASTA Sequence Processing Tool Suite

A suite of desktop applications for processing FASTA format biological sequence files, featuring two core tools: sequence comparison and sample classification.

---

## 📦 Tool Overview

### 1. FASTA Sequence Comparison Tool
For comparing multiple FASTA files to identify common and unique sequences.

**Main Features:**
- Compare multiple FASTA files simultaneously
- Identify sequences common to all files
- Identify sequences unique to each file
- Provide detailed statistical information
- Support sequence viewing and export

### 2. FASTA Sample Classification Tool
Classify sequences in mixed FASTA files by sample names.

**Main Features:**
- Automatically recognize sample names from sequence headers
- Organize sequences by sample into separate files
- Support selective sample export
- Provide sample statistics and preview features

---

## 🚀 Quick Start

### System Requirements
- Python 3.6 or higher
- Windows / macOS / Linux
- No additional dependencies required

### Installation Steps

1. **Check Python Environment**
   ```bash
   python --version
   ```

2. **Launch Tools**
   - Double-click `fasta_compare.bat` to launch sequence comparison tool
   - Double-click `fasta_classifier.bat` to launch sample classification tool

---

## 📚 Documentation Navigation

### Quick Start
- **[Quick Start Guide](./quick_start.md)** - 5-minute quick start with complete examples

### Detailed Documentation
- **[User Manual](./user_manual.txt)** - Complete usage manual with detailed feature explanations
- **[FASTA Sequence Comparison Tool README](./fasta_compare_README.md)** - Detailed documentation for sequence comparison tool
- **[FASTA Sample Classification Tool README](./fasta_classifier_README.md)** - Detailed documentation for sample classification tool

---

## 📁 File List

```
comparison/
├── README.md                        # This file (main project documentation)
├── quick_start.md                   # Quick start guide
├── user_manual.txt                  # Complete user manual
├── fasta_compare_README.md          # Sequence comparison tool documentation
├── fasta_classifier_README.md       # Sample classification tool documentation
├── fasta_compare.py                 # Sequence comparison tool source code
├── fasta_classifier.py              # Sample classification tool source code
├── fasta_compare.bat                # Sequence comparison tool launch script
└── fasta_classifier.bat             # Sample classification tool launch script
```

---

## 💡 Usage Scenarios

### Scenario 1: Multi-Species Sequence Comparison
**Goal**: Find genes common and unique to multiple species

**Workflow**:
1. Use sequence comparison tool to compare CDS files from multiple species
2. View common sequences (conserved genes)
3. View unique sequences (species-specific genes)
4. Export results for further analysis

### Scenario 2: Sample Sequence Classification
**Goal**: Classify mixed sequences by sample name

**Workflow**:
1. Use sample classification tool to load mixed FASTA file
2. Program automatically recognizes sample names
3. Export all samples or selectively export some samples
4. Get sample-organized sequence files

### Scenario 3: Complete Analysis Workflow
**Goal**: From multiple species files to sample-organized common genes

**Workflow**:
1. Use sequence comparison tool to compare multiple species files
2. Export common sequences (genes common to all species)
3. Use sample classification tool to load exported common sequences
4. Classify by sample name and export
5. Get sample-organized common gene sequences

---

## 🔧 Core Features

### FASTA Sequence Comparison Tool
| Feature | Description |
|---------|-------------|
| Multi-file loading | Load multiple FASTA files simultaneously |
| Common sequence identification | Find sequences common to all files |
| Unique sequence identification | Find sequences unique to each file |
| Statistics | Provide detailed statistical data |
| Sequence viewing | View detailed information of any sequence |
| Result export | Export common or unique sequences |

### FASTA Sample Classification Tool
| Feature | Description |
|---------|-------------|
| Auto recognition | Automatically extract sample names from headers |
| Sample classification | Organize sequences by sample name |
| Selective export | Support exporting partial or all samples |
| Sample statistics | Display sequence count and length for each sample |
| Sequence preview | Preview sample content before export |

---

## 📊 Data Flow Diagram

```
┌─────────────────────┐
│  Multiple FASTA files│
│  (species1.fasta)  │
│  (species2.fasta)  │
│  (species3.fasta)  │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  FASTA Sequence     │
│  Comparison Tool    │
│  - Compare seqs     │
│  - Identify common/ │
│    unique           │
└──────────┬──────────┘
           │
           ▼
    ┌──────┴──────┐
    │             │
    ▼             ▼
┌────────┐   ┌──────────┐
│Common  │   │Unique    │
│Sequences│  │Sequences │
└───┬────┘   └────┬─────┘
    │             │
    ▼             ▼
┌─────────────────────────┐
│  FASTA Sample          │
│  Classification Tool   │
│  - Classify by sample  │
│  - Export sample files │
└──────────┬──────────────┘
           │
           ▼
┌─────────────────────────┐
│  Sample-organized       │
│  sequence files         │
│  (A10.fasta, A11.fasta) │
└─────────────────────────┘
```

---

## 🎯 Typical Applications

### Bioinformatics Research
- Comparative genomics
- Population genetics analysis
- Phylogenetic analysis
- Gene family studies

### Data Preprocessing
- FASTA file formatting
- Sequence quality control
- Data organization and archiving

### Teaching and Learning
- Bioinformatics courses
- Sequence analysis teaching
- Data processing practice

---

## ❓ FAQ

### Q1: What's the difference between the two tools?
**A**: The sequence comparison tool compares multiple files to find common and unique sequences; the sample classification tool classifies sequences in a mixed file by sample.

### Q2: Which FASTA formats are supported?
**A**: Supports standard FASTA format, header starts with ">", sequence starts after new line. Supports .fasta and .fa extensions.

### Q3: Is comparison based on sequence content or header?
**A**: The sequence comparison tool compares based on header names, not sequence content.

### Q4: Can I use it on Mac or Linux?
**A**: Yes, but need to run Python scripts directly:
```bash
python fasta_compare.py
python fasta_classifier.py
```

### Q5: How to check Python version?
**A**: Enter in command line: `python --version`

---

## 🔍 Sequence Name Recognition Rules

### Sequence Comparison Tool
Program extracts sequence names from headers, removing position information:
- `rps12_join{...}|A10.fasta` → `rps12`
- `psbA_[76:1137](-)` → `psbA`
- `gene_name[some_info]` → `gene_name`

### Sample Classification Tool
Program extracts sample names from end of headers:
- `>gene_name|A10.fasta` → `A10`
- `>sequence|sample1` → `sample1`
- `>name A10.fasta` → `A10`

---

## 📝 Version Information

- **FASTA Sequence Comparison Tool**: v1.0.0
- **FASTA Sample Classification Tool**: v1.0.0
- **Last Updated**: 2025

---

## 🤝 Technical Support

For questions or suggestions, please:
1. Check detailed documentation
2. Review troubleshooting section in user manual
3. Check FAQ

---

## 📄 License

This tool is open source software, free to use and modify.

---

## 🙏 Acknowledgments

Thanks to all users and developers using this tool!

---

**Get Started:** Check [Quick Start Guide](./quick_start.md) to start immediately!

**Detailed Docs:** Check [User Manual](./user_manual.txt) for more details.
