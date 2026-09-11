# FASTA Tools Quick Start Guide

Welcome to the FASTA sequence processing tool suite! This guide will help you quickly get started with these two tools.

---

## Tool Overview

### 1. FASTA Sequence Comparison Tool
For comparing multiple FASTA files to find common and unique sequences.

### 2. FASTA Sample Classification Tool
Classify mixed FASTA files by sample names into multiple files.

---

## Quick Start

### Step 1: Check Python Environment
Open command line and enter:
```bash
python --version
```
Ensure version is 3.6 or higher. If Python is not installed, download and install from https://python.org

### Step 2: Launch Tools
Double-click the corresponding .bat file:
- `fasta_compare.bat` - Launch sequence comparison tool
- `fasta_classifier.bat` - Launch sample classification tool

---

## Tool 1: FASTA Sequence Comparison Tool - 5-Minute Quick Start

### Scenario: Compare CDS sequences from 3 species

**1. Prepare Data**
```
CDS/
  ├── species1.fasta
  ├── species2.fasta
  └── species3.fasta
```

**2. Load Files**
1. Launch the program
2. Click "Select Folder"
3. Select `CDS/` folder
4. Program automatically loads all .fasta files

**3. Start Comparison**
1. Click "Start Comparison"
2. Wait a few seconds

**4. View Results**

**Common Sequences** (genes shared by all species):
- Switch to "Common Sequences" tab
- View sequences common to all species (e.g., rbcL, matK)
- Double-click any sequence to view details

**Unique Sequences** (genes unique to each species):
- Switch to "Unique Sequences" tab
- View sequences unique to each species

**Statistics**:
- Switch to "Statistics" tab
- View detailed statistical distribution

**5. Export Results**
- Click "Export Common Sequences" → Save as `common.fasta`
- Click "Export Unique Sequences" → Save to `output/` folder

---

## Tool 2: FASTA Sample Classification Tool - 5-Minute Quick Start

### Scenario: Classify mixed samples by sample name

**1. Prepare Data**
You have a `mixed_sequences.fasta` file containing sequences from multiple samples:
```
>rps12|A10.fasta
ATGCGT...
>cemA|A11.fasta
ATGCGT...
>psbA|A12.fasta
ATGCGT...
```

**2. Load File**
1. Launch the program
2. Click "Select File"
3. Select `mixed_sequences.fasta`

**3. View Samples**
- Left side shows all recognized samples: A10, A11, A12
- Shows sequence count for each sample

**4. Export Samples**

**Export All Samples**:
1. Click "Export All Samples"
2. Select save folder
3. Program creates: `A10.fasta`, `A11.fasta`, `A12.fasta`

**Export Selected Samples**:
1. Check only needed samples (e.g., A10 and A11)
2. Click "Export Selected Samples"
3. Program only creates: `A10.fasta`, `A11.fasta`

**5. Preview Sequences**
- Click sample name (e.g., A10)
- Right side shows all sequences for that sample
- Click sequence to view complete content

---

## Complete Workflow Example

### From multiple species files to sample-organized sequences

**Step 1: Use sequence comparison tool**
```
Input: FASTA files from 3 species
  ├── species1.fasta
  ├── species2.fasta
  └── species3.fasta

Operation: Compare files, export common sequences

Output: common.fasta (sequences common to all species)
```

**Step 2: Use sample classification tool**
```
Input: common.fasta (sequences from multiple samples)
  Format: >gene_name|sample.fasta

Operation: Classify by sample name

Output: Sample-organized sequence files
  ├── A10.fasta
  ├── A11.fasta
  └── A12.fasta
```

**Result**: Now you have sample-organized common gene sequences for further analysis.

---

## Common Tasks Quick Reference

| Task | Tool | Steps |
|------|------|-------|
| Compare multiple FASTA files | Sequence Comparison Tool | Select folder → Start comparison → View results |
| Find common sequences | Sequence Comparison Tool | After comparison, switch to "Common Sequences" tab |
| Find unique sequences | Sequence Comparison Tool | After comparison, switch to "Unique Sequences" tab |
| Classify by samples | Sample Classification Tool | Select file → View samples → Export |
| Export specific samples | Sample Classification Tool | Check samples → Export selected samples |
| View sequence details | Either tool | Use "Sequence Viewer" function |
| Export results | Either tool | Click "Export" button |

---

## UI Element Quick Reference

### Sequence Comparison Tool Interface
```
┌─────────────────────────────────────────┐
│ [Select Folder] [Start] [Export] [Clear] │  ← Toolbar
├─────────────────────────────────────────┤
│ Loaded Files:                            │  ← File list
│ ┌─────────────────────────────────────┐ │
│ │ species1.fasta                      │ │
│ │ species2.fasta                      │ │
│ └─────────────────────────────────────┘ │
├─────────────────────────────────────────┤
│ [Common] [Unique] [Stats] [View]        │  ← Tabs
│ ┌─────────────────────────────────────┐ │
│ │ Seq Name    Count    Files          │ │
│ │ rbcL        3         sp1,sp2,sp3  │ │
│ └─────────────────────────────────────┘ │
├─────────────────────────────────────────┤
│ Select folder containing FASTA files    │  ← Status bar
└─────────────────────────────────────────┘
```

### Sample Classification Tool Interface
```
┌─────────────────────────────────────────┐
│ [Select File] [Export] [Select All] [Clear] │  ← Toolbar
├─────────────────────────────────────────┤
│ File: mixed_sequences.fasta             │  ← File info
├────────────────┬────────────────────────┤
│ Sample List     │ Sequence Details       │  ← Split view
│ ☐ A10 (50)     │ Sequences for sample:  │
│ ☐ A11 (45)     │ ├ rps12                │
│ ☐ A12 (52)     │ ├ cemA                 │
│                │ └ psbA                 │
│                │ Sequence content:       │
│                │ >rps12|A10.fasta       │
│                │ ATGCGT...              │
└────────────────┴────────────────────────┘
```

---

## Quick Operation Tips

### Sequence Comparison Tool
- Double-click sequence → View details
- Enter key → Confirm action
- Esc key → Cancel action (for some functions)

### Sample Classification Tool
- Click checkbox → Toggle selection state
- Click sample name → Preview sequences
- Select All button → Select all samples with one click

---

## Common Issues Quick Resolution

**Problem: Program won't start**
→ Check Python installation: `python --version`

**Problem: Can't find FASTA files**
→ Ensure file extension is .fasta or .fa

**Problem: Display garbled text**
→ Convert file encoding to UTF-8

**Problem: Can't recognize sample names**
→ Check header format, ensure sample identifiers are present

**Problem: Slow processing**
→ Wait for program to complete, large files require more time

---

## Output File Description

### Sequence Comparison Tool Output
- `common.fasta` - Sequences common to all files
- `filename_unique.fasta` - Unique sequences for each file

### Sample Classification Tool Output
- `sample_name.fasta` - Sequence file for each sample

---

## Next Steps

- Check detailed `user_manual.txt` for complete functionality
- Check `fasta_compare_README.md` for sequence comparison tool details
- Check `fasta_classifier_README.md` for sample classification tool details
- Practice to familiarize yourself with the tools

---

## Need Help?

Encountering problems?
1. Check troubleshooting section in user manual
2. Check FAQ
3. Verify file format is correct
4. Validate Python environment

---

**Happy using!**

---

*Version: v1.0.0 | Updated: 2025*
