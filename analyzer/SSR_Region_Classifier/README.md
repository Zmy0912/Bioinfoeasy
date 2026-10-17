# SSR Region Classifier

A graphical (GUI) tool that classifies every SSR into **CDS / Intron / IGS** based on its coordinate in the genome, using a MISA result file together with the matching GenBank (GB) annotation files. The outputs are ready-to-plot tables for Excel, Origin, R or Python.

## Features

- **Intuitive GUI** built on Tkinter — no command line required.
- **Automatic species matching**: the first column of the MISA file (`ID`, e.g. `NC01`) is matched against the basename of each GB file (`NC01.gb`).
- **Region classification rules**:
  - **CDS** — the SSR (or its midpoint) lies inside a CDS exon.
  - **Intron** — the SSR lies inside an intron inferred from the `gene` span minus its CDS exons.
  - **IGS** — anything else is treated as intergenic spacer.
- **Multiple output formats**: detail table, long-form summary, wide-form summary, region × SSR-type matrix — all directly usable for plotting.

## Usage

### 1. Launch

Double-click `Run_analysis.bat`, or run from the command line:

```bash
python ssr_region_classifier.py
```

### 2. In the GUI

1. Click **Browse...** to select the **MISA file** (e.g. `0142.MISA`).
2. Select the **GB folder** (the sub-directory named after the MISA file, containing `NC01.gb`, `NC02.gb`, ...).
3. Select the **Output directory** (a default is suggested automatically).
4. Click **Run analysis**.

### 3. Inspect the results

When finished, click **Open output folder** to view the generated files.

## Output files

| File | Description |
| --- | --- |
| `SSR_with_region.txt` | MISA contents plus two new columns: `Region` (CDS / Intron / IGS) and `RegionDetail` (the interval it falls into). |
| `SSR_region_summary.tsv` | Long form: `Species, Region, Count, Percent` — ideal for ggplot2 stacked bars. |
| `SSR_region_summary_wide.tsv` | Wide form: `Species, CDS, Intron, IGS, Total`. |
| `SSR_type_by_region.tsv` | Rows = region, columns = SSR type (`p1 / p2 / ... / c`) — convenient for grouped bar charts. |
| `log.txt` | Full run log. |

## Input format

### MISA file
- Tab-separated (`.misa` or `.txt`); the header must contain at least
  `ID, SSR nr., SSR type, SSR, size, start, end`.
- The first column `ID` (e.g. `NC01`) must match the basename of the GB file (`NC01.gb`).

### GB folder
- Each GB file must contain `gene` and `CDS` features in its FEATURES table.
- Complex locations such as `complement(join(...))` are supported.

## Dependencies

- Python 3.8+
- Only the standard library (`tkinter`, `csv`, `re`, `threading`, ...) is used — **no extra packages required**.

## Directory layout (example)

```
0919\
├─ 0142.MISA                ← input MISA file
├─ 0142\                    ← GB folder
│   ├─ NC01.gb
│   ├─ NC02.gb
│   └─ ...
└─ SSR_Region_Classifier\   ← this program
    ├─ ssr_region_classifier.py
    └─ Run_analysis.bat
```