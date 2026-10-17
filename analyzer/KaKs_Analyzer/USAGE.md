# User Guide

This guide walks you through every step of using the KaKs analysis tool.

---

## 1. Installation & Launch

### 1.1 Install Dependencies

```bash
pip install pandas openpyxl
```

> Windows users who see `tkinter not found` should reinstall Python
> from [python.org](https://www.python.org/downloads/) and make sure
> "tcl/tk and IDLE" is checked.

### 1.2 Launch the Program

```bash
cd <your_project_dir>
python kaks_analyzer.py
```

You can also double-click `kaks_analyzer.py` (provided Windows
associates `.py` with Python).

---

## 2. Interface Overview

```
+-------------------------------------------------------------+
| KaKs Result Analysis & Filtering Tool v2.0                  |
+-------------------------------------------------------------+
| (1) Select Input File                                       |
|   [<your_project_dir>\example_kaks.xlsx] [Browse...] |
|                                          [Load] [Example]   |
+-------------------------------------------------------------+
| (2) Filter Parameters                                       |
|   Missing ratio:         Columns to impute:                 |
|   o 20% (drop >=20%)   [ ] Impute Ka (NaN -> median)        |
|   o 30% (drop >=30%)   [ ] Impute Ks (NaN -> median)        |
|   * 50% (drop >=50%)   [x] Impute Ka/Ks (NaN/0 -> median)   |
|   o 70% (drop >=70%)   * Ka/Ks == 0 treated as missing      |
+-------------------------------------------------------------+
| (3) Output Directory                                        |
|   [<your_project_dir>]              [Browse...]    |
+-------------------------------------------------------------+
| [Start Analysis] [Clear Log] [Open Output Folder]  [About]  |
| ooooooooooooooooooooooooooooooooooooooooooooooooooooo 100%  |
+-------------------------------------------------------------+
| [Run Log] [Report Preview] [Data Preview]                   |
| +---------------------------------------------------------+ |
| | 22:15:32 Loaded: 6888 rows, 8 genes, 41 species         | |
| | 22:15:33 Data preview updated.                          | |
| | 22:15:40 Starting | threshold=50% | impute=['Ka/Ks']    | |
| | 22:15:41 Done!                                          | |
| |   Total genes: 8  Kept: 5  Removed: 3                   | |
| |   Generated files:                                      | |
| |     - example_kaks_filtered_50pct_*.xlsx                | |
| |     - example_kaks_filtered_50pct_*_report.xlsx         | |
| |     - example_kaks_filtered_50pct_*_summary.txt         | |
| +---------------------------------------------------------+ |
+-------------------------------------------------------------+
```

---

## 3. Step-by-step

### Step 1: Choose a File

#### Option A: Use Example (recommended for first use)

Click **"Use Example"**; the program fills in the path to
`example_kaks.xlsx` and loads it automatically.

#### Option B: Pick Manually

Click **"Browse..."** and select one of:

- a `.xlsx` file — Excel organized sheet (6 columns: Sp1, Sp2, Gene,
  Ka, Ks, Ka/Ks)
- a `.kaks` file — raw tab-separated output from KaKs Calculator
  (22 columns, first column `Sequence`)

#### After loading

- The data is loaded and deduplicated automatically.
- The **Data Preview** tab shows the first 500 rows of the input.
- The **Run Log** tab prints:
  `Loaded successfully: X rows, Y genes, Z species.`

### Step 2: Set Filter Parameters

#### Missing-ratio Threshold (pick one)

| Option | Keep condition |
|--------|----------------|
| 20%    | Ka/Ks missing ratio < 20% |
| 30%    | Ka/Ks missing ratio < 30% |
| **50%** (default) | Ka/Ks missing ratio < 50% |
| 70%    | Ka/Ks missing ratio < 70% |

#### Median Imputation Columns (multi-select)

- ☑ **Impute Ka/Ks** (default): replaces missing entries with the
  median of valid Ka/Ks values from the same gene's other species
- ☐ Impute Ka: replaces NaN in Ka with the median of valid Ka values
- ☐ Impute Ks: replaces NaN in Ks with the median of valid Ks values

> **Tip:** Usually ticking only Ka/Ks is sufficient. Ka and Ks columns
> rarely contain zeros, so imputation is less meaningful for them.

### Step 3: Choose Output Directory

Defaults to the folder containing the program. Click **"Browse..."** to
pick a different one.

### Step 4: Start Analysis

Click **"Start Analysis"**; the program runs the following steps in a
background thread:

1. Load and deduplicate the data
2. Compute the per-gene missing ratio
3. Drop genes over the threshold
4. Impute the remaining missing values
5. Write the three output files

The progress bar updates in real time; when finished the log prints a
result summary.

### Step 5: Inspect Results

- **Run Log** — processing steps and final summary
- **Report Preview** — per-gene statistics table
- **Output Folder** — three files (xlsx data + xlsx report + txt
  summary)

---

## 4. Output Files Explained

### 4.1 Data File `<name>_filtered_<threshold>pct_<timestamp>.xlsx`

Filtered and imputed data. Same layout as the input:

| Species1 | Species2 | Gene | Ka | Ks | Ka/Ks |
|----------|----------|------|-----|-----|-------|
| A.acutifolius | A.adscendens | accD | 0.019657 | 0.009183 | 2.140640 |
| A.acutifolius | A.cooperi    | accD | 0.031616 | 0.014280 | 2.213960 |
| ...      | ...      | ...  | ... | ... | ...   |

> Note: previously-missing (NaN or 0) cells now contain the imputed
> median.

### 4.2 Report File `<name>_filtered_<threshold>pct_<timestamp>_report.xlsx`

One row per gene:

| Gene | Total | MissingCount | MissingKa | MissingKs | MissingRate | FilledCount | MedianKa | MedianKs | MedianKaKs | Status |
|------|-------|--------------|-----------|-----------|-------------|-------------|----------|----------|------------|--------|
| accD | 861 | 45 | 22 | 35 | 5.23% | 816 | 0.0221 | 0.01678 | 1.252 | Kept |
| atpH | 861 | 861 | 861 | 669 | 100.00% | 0 | NA | 0.01757 | NA | Removed |
| ...  | ... | ...   | ...   | ...   | ...         | ... | ...    | ...    | ...        | ...    |

### 4.3 Text Summary `<name>_filtered_<threshold>pct_<timestamp>_summary.txt`

A human-readable text report with the following structure:

```
================================================================================
KaKs Analysis & Filtering Report
================================================================================
Source file       : example_kaks.xlsx
Generated at      : 2026-09-21 21:30:00
Missing ratio     : < 50%
Missing rule      : Ka/Ks NaN or == 0 (Ka/Ks==0 is treated as a degenerate
                    value produced when Ka fails to compute)
Imputation        : Median of valid values from the same gene's other species
--------------------------------------------------------------------------------
Total genes       : 8
Kept genes        : 5
Removed genes     : 3
--------------------------------------------------------------------------------
Removed genes (3):
   - atpF
   - atpH
   - atpI
--------------------------------------------------------------------------------
Per-gene statistics:
Gene       Total MissKa/Ks MissKa MissKs    Rate +Filled   MedKa   MedKs MedKa/Ks    Status
accD         861        45     22     35   5.23%     816   0.0221  0.0168     1.252      Kept
...
================================================================================
Field descriptions:
  Total       : Number of rows for the gene (i.e. species pairs).
  MissKa/Ks   : Missing count for the Ka/Ks column (NaN or == 0).
  MissKa      : Missing count for the Ka column (NaN).
  MissKs      : Missing count for the Ks column (NaN).
  Rate        : Ka/Ks missing ratio = MissKa/Ks / Total.
  +Filled     : Number of valid rows after imputation.
  MedKa/MedKs/MedKa/Ks : Median used for imputation (computed over valid values).
  Status      : Kept or Removed.
```

---

## 5. Input File Requirements

### 5.1 Excel Organized Sheet (`.xlsx`)

Columns (in this order):

| # | Column     | Type   | Description |
|---|------------|--------|-------------|
| 1 | Species1   | text   | Species 1 (reference species) |
| 2 | Species2   | text   | Species 2 (compared species) |
| 3 | Gene       | text   | Gene name |
| 4 | Ka         | number | Non-synonymous substitution rate |
| 5 | Ks         | number | Synonymous substitution rate |
| 6 | Ka/Ks      | number | Ratio of Ka to Ks |

> If the first column contains rows with `Sequence` / NaN / None (e.g.
> summary lines), the tool automatically skips them.

### 5.2 KaKs Raw Output (`.kaks`)

A 22-column tab-separated file produced by KaKs_calculator. The first
column `Sequence` uses the following format:

```
Species1<id>&Species2<id>|gene[_copy2][_mafft]
```

Example:

```
Asparagus_acutifolius_176986014997886&Asparagus_filicinus_17698599752654|accD_mafft
Asparagus_acutifolius_176986014997886&Asparagus_filicinus_17698599752654|ndhB_copy2_mafft
```

The tool automatically extracts `Species1` / `Species2` / `Gene` and
normalizes gene names:

- `ndhB_copy2_mafft` → `ndhB`
- `accD_mafft` → `accD`

---

## 6. Advanced Tips

### 6.1 Command-Line / Batch Use

The `KaksProcessor` class in `kaks_analyzer.py` can be imported and
called directly:

```python
from kaks_analyzer import KaksProcessor

df = KaksProcessor.load_file("data.xlsx")
result = KaksProcessor.analyze(df, threshold=0.30, fill_columns=["Ka/Ks"])
files = KaksProcessor.export(result, "output/", threshold=0.30,
                              source_file="data.xlsx")
print(result["total_genes"], result["kept_genes"])
```

### 6.2 Running Multiple Times

Every "Start Analysis" click re-runs the analysis on the currently
loaded data **without modifying the source file**. To switch datasets,
click **"Browse..."** and pick a different file.

### 6.3 Generating Reports for All Four Thresholds

Click "Start Analysis" four times in a row, once per threshold. Output
filenames include a timestamp, so previous outputs are not overwritten.

---

## 7. Troubleshooting

| Problem | Solution |
|---------|----------|
| `ModuleNotFoundError: No module named 'pandas'` on launch | Run `pip install pandas openpyxl` |
| `tkinter not found` on launch | Reinstall Python and tick the tcl/tk component |
| "Missing 'Sequence' column" while loading | Check that the first line of the .kaks file is `Sequence\tMethod\t...` |
| All `MedianKa/Ks` in the report are 0 | Confirm the dataset actually contains valid Ka/Ks values (>0) |
| Program hangs / no response | Inspect the raw file for unusual characters in a text editor first |

---

If you have questions or suggestions, please inspect the source files
and re-run the tool.