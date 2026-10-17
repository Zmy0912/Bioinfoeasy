# KaKs Result Analysis & Filtering Tool

> A graphical tool for filtering KaKs calculation results and imputing
> missing values. It auto-detects both the organized `.xlsx` spreadsheet
> and the original `.kaks` text output produced by KaKs Calculator,
> groups rows by gene, drops genes whose missing-value ratio exceeds a
> user-selected threshold, and imputes missing values with the per-gene
> median of valid values from other species.

---

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Core Algorithm](#core-algorithm)
- [Output Files](#output-files)
- [Project Structure](#project-structure)
- [FAQ](#faq)
- [Changelog](#changelog)

---

## Features

- ✅ **Dual format auto-detection**: supports both `.xlsx` (6-column
  organized sheet) and `.kaks` (KaKs Calculator's original 22-column
  tab-separated output)
- ✅ **Graphical interface**: tkinter-based GUI, no command line required
- ✅ **4 preset thresholds**: 20% / 30% / 50% / 70% missing-ratio cutoffs
- ✅ **Smart missing-value detection**: treats `Ka/Ks == 0` as missing
  (KaKs software emits 0 when Ka computation fails; the value is not
  biologically meaningful)
- ✅ **Median imputation**: uses the per-gene median of **valid values
  (>0)** so the imputed result is no longer 0
- ✅ **Automatic deduplication**: handles duplicate
  `(Species1, Species2, Gene)` rows by preferring the one with a valid
  Ka/Ks value
- ✅ **Three outputs per analysis**: filtered data (xlsx) + per-gene
  report (xlsx) + concise text report (txt)

---

## Quick Start

### Requirements

- Python 3.8+
- Dependencies: `pandas`, `openpyxl` (tkinter ships with the standard
  library)

### Install Dependencies

```bash
pip install pandas openpyxl
```

### Run the Program

```bash
cd <your_project_dir>
python kaks_analyzer.py
```

Once launched, click **"Use Example"** to load the bundled
`example_kaks.xlsx` and try the tool immediately.

---

## Usage

### 1. Select an Input File

The tool supports two formats, automatically detected by extension:

| Extension       | Format                  | Description |
|-----------------|-------------------------|-------------|
| `.xlsx` / `.xls` | Organized spreadsheet   | Columns: `Species1, Species2, Gene, Ka, Ks, Ka/Ks` |
| `.kaks` / `.tsv` / `.txt` | KaKs raw output | 22 tab-separated columns; the first column `Sequence` contains `Sp1&Sp2\|gene_mafft` |

Click **"Browse..."** to pick a file, or **"Use Example"** to load the
bundled example.

### 2. Set Filter Parameters

- **Missing-ratio threshold** (4 radio buttons):
  - `20%` — genes whose Ka/Ks missing ratio is ≥ 20% are dropped
  - `30%` — genes whose Ka/Ks missing ratio is ≥ 30% are dropped
  - `50%` — genes whose Ka/Ks missing ratio is ≥ 50% are dropped
  - `70%` — genes whose Ka/Ks missing ratio is ≥ 70% are dropped
- **Columns to impute with median**: tick `Ka` / `Ks` / `Ka/Ks`. By
  default only `Ka/Ks` is ticked.

> **Why is Ka/Ks == 0 treated as missing?**
> In KaKs Calculator, when Ka cannot be computed the software outputs
> 0 as a placeholder. Such values have no biological meaning and
> severely drag the median downward, so they are treated as missing.
> The median is computed only over `Ka/Ks > 0` valid values.

### 3. Choose Output Directory

Click **"Browse..."** to select where to save results; defaults to the
folder containing the program.

### 4. Click **"Start Analysis"**

The analysis runs in a background thread while the progress bar updates
in real time. When finished:

- **Run Log** tab: processing steps and result summary
- **Report Preview** tab: tabular view of per-gene statistics
- **Data Preview** tab: first 500 rows of the raw input

### 5. View the Output Files

When processing completes, the log prints the full paths of three output
files. Click **"Open Output Folder"** to inspect them.

---

## Core Algorithm

### Missing-value Detection

| Column   | Missing definition |
|----------|---------------------|
| `Ka/Ks`  | NaN **or** == 0 (Ka computation failed) |
| `Ka`     | NaN |
| `Ks`     | NaN |

### Missing Ratio

$$
\text{MissingRatio} = \frac{\text{Missing Ka/Ks count}}{\text{Total rows}}
$$

### Median Imputation

For each kept gene, the median of every target column is computed using
**valid values only**:

| Column   | Valid-value definition |
|----------|------------------------|
| `Ka/Ks`  | `> 0` |
| `Ka`     | not NaN |
| `Ks`     | not NaN |

Missing entries (NaN or 0) in the column are then replaced by that
median.

### Automatic Deduplication

Rows are deduplicated by `(Species1, Species2, Gene)`, preferring the
row whose Ka/Ks is a valid value (>0). For example, `rps7` appears in 82
rows in the raw data because of the `_copy2` variant and is collapsed to
41 rows automatically.

### Gene-name Normalization

Original `.kaks` gene names look like `ndhB_copy2_mafft`.
After normalization they become `ndhB` (both the `_mafft` and `_copy\d+`
suffixes are stripped).

---

## Output Files

Each run writes three files:

| File | Contents |
|------|----------|
| `<name>_filtered_<threshold>pct_<timestamp>.xlsx` | Filtered + imputed data |
| `<name>_filtered_<threshold>pct_<timestamp>_report.xlsx` | Per-gene statistics |
| `<name>_filtered_<threshold>pct_<timestamp>_summary.txt` | Concise text report (with field descriptions) |

### Report Fields

| Field | Description |
|-------|-------------|
| `Gene` | Gene name |
| `Total` | Number of rows for the gene (i.e. species pairs) |
| `MissKa/Ks` | Missing count for Ka/Ks column (NaN or == 0) |
| `MissKa` | Missing count for Ka column (NaN) |
| `MissKs` | Missing count for Ks column (NaN) |
| `Rate` | Ka/Ks missing ratio = MissKa/Ks / Total |
| `+Filled` | Number of valid rows after imputation |
| `MedKa` | Median of Ka column (over valid values) |
| `MedKs` | Median of Ks column (over valid values) |
| `MedKa/Ks` | Median of Ka/Ks column (over valid values) |
| `Status` | Kept or Removed |

---

## Project Structure

```
<your_project_dir>\
├── kaks_analyzer.py          # Main program (GUI + core logic)
├── example_kaks.xlsx         # Example file (organized sheet)
├── example_kaks.kaks         # Example file (raw KaKs format)
├── make_example.py           # Script to regenerate the example files
├── USAGE.md                   # Detailed user guide
├── CHANGELOG.md              # Version history
├── requirements.txt          # Python dependencies
├── merge.axt.kaks            # Raw KaKs output (user-provided)
├── 0921shaixuan.xlsx         # Organized xlsx (user-provided)
└── test_*.py / inspect_*.py  # Test scripts (optional)
```

---

## FAQ

### Q1: Why is Ka/Ks == 0 treated as missing?

**A:** KaKs Calculator outputs 0 as a placeholder when Ka cannot be
computed; the value has no biological meaning. If 0 is treated as a
valid value, the median is dragged way down (for `accD` it drops from
1.252 to 0.27) and the imputation result is still 0 — defeating the
purpose of imputation.

### Q2: Why does the example dataset contain only 8 genes?

**A:** To keep the demo quick and readable, the example files contain
only the first 8 genes, but they cover a variety of missing ratios
(5% – 100%). To test against the full 79 genes, load the raw file
`merge.axt.kaks` or `0921shaixuan.xlsx`.

### Q3: How should I choose the threshold?

**A:** It depends on how strict your downstream analysis needs to be:

- **20%** — Strict; the kept data has high completeness (best for
  downstream analyses such as phylogenetic reconstruction)
- **30%** — Moderate (commonly used in KaKs studies)
- **50%** — Lenient; keeps more genes
- **70%** — Most lenient; drops only the genes that are almost
  completely missing

### Q4: Median imputation vs. deletion?

**A:** Median imputation introduces some statistical bias (it
underestimates variance) but preserves more genes for downstream
analysis. If the sample size is large enough, dropping missing genes
without imputation may be preferable.

### Q5: How large a file can the tool handle?

**A:** Tested with a 70,000-row `.kaks` file: loading takes about 1 s
and processing about 0.5 s. For files much larger than that
(millions of rows), consider chunked processing.

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md).

---

## Dependencies

```
pandas >= 1.3
openpyxl >= 3.0
```

The GUI uses Python's standard library `tkinter`, no extra install
needed.