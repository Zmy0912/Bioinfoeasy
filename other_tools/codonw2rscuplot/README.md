# CodonW2RSCUPlot

> Aggregate the RSCU output of **codonW** (or MEGA) into the input format
> required by the [**Codon Usage Visualization / RSCU-Plot**](https://pcg-lab.shinyapps.io/RSCU-Plot/)
> web application.

A lightweight, GUI-based Python tool with **zero third-party dependencies**.

---

## 1. Overview

The [RSCU-Plot](https://pcg-lab.shinyapps.io/RSCU-Plot/) web app only
handles **visualisation** (heatmaps, clustering, PCA, etc.); it does NOT
compute RSCU itself. It expects a **wide-table CSV** with the structure:

| Col 1 | Col 2 | Col 3 ... |
|-------|-------|-----------|
| Amino Acid | Codon | RSCU value of each species |

But in practice codonW gives us **one file per species** (or per sequence),
which cannot be uploaded directly. This tool bridges the gap:

```
RSCU files per species (.blk / .out / .txt / .csv)
                  |
                  |  CodonW2RSCUPlot (this tool)
                  v
             Wide-table CSV
                  |
                  v
  RSCU-Plot (https://pcg-lab.shinyapps.io/RSCU-Plot)
```

Features:

- Point-and-click GUI (tkinter) -- no command-line skills required.
- Zero third-party dependencies (Python standard library only).
- Auto-detects both **MEGA `.blk`** and **codonW `.out` / `.txt` / `.csv`** formats.
- Auto-converts DNA letter `T` to RNA letter `U`, so files from different
  sources can be merged uniformly.
- File name (without extension) is used as the species column name.
- Rows are ordered by the standard codon table (matching the official
  template: Phe -> Leu -> Ile -> ... -> Gly).
- Empty / unparsable files are reported clearly and kept as empty columns.
- Runs in a background thread -- the UI stays responsive.

---

## 2. Requirements

- **Python >= 3.7**
- Built-in `tkinter` (bundled with official installers on Windows / macOS;
  on Linux you may need `sudo apt install python3-tk`).
- No `pip install` needed.

---

## 3. Usage

### 3.1 Prepare input files

Put one RSCU file per species into the same folder, e.g.:

```
my_rscu_results/
├── Species_A.blk
├── Species_B.blk
├── Species_C.out
└── Species_D.txt
```

The file name (without extension) becomes the species column name in the
output CSV.

### 3.2 Run the GUI

```bash
python codonw2rscuplot.py
```

A window appears. Steps:

1. Click **Browse...** to pick the folder containing the RSCU files.
2. Set the output CSV path (default: `rscu_input_for_plot.csv` in the
   same folder).
3. (Optional) Uncheck *Include stop codons (TER)* if you want to drop them.
4. Click **Convert**. Progress and the final summary appear in the log
   area below.

### 3.3 Upload to RSCU-Plot

Open <https://pcg-lab.shinyapps.io/RSCU-Plot/> and in the *Data Upload*
panel:

- File type: **CSV**, separator: **Comma**.
- Tick *Header row*.
- Upload the CSV generated above.
- Switch to the *Data Preview* tab to verify, then go to *Visualisation*.

---

## 4. Input format

### 4.1 MEGA-style (`.blk`)

Four codon entries per line; the first row of each column has the amino
acid label, the second does not:

```
Phe UUU  791 1.27 Ser UCU  469 1.69 Tyr UAU  656 1.61 Cys UGU  197 1.54
    UUC  450 0.73     UCC  274 0.99     UAC  159 0.39     UGC   59 0.46
Leu UUA  754 1.96     UCA  335 1.21 TER UAA   36 1.40 TER UGA   20 0.78
    UUG  476 1.24     UCG  161 0.58     UAG   21 0.82 Trp UGG  393 1.00
...
22178 codons in NC01_joined_from (used Universal Genetic code)
```

### 4.2 codonW-style (`.out` / `.txt` / `.csv`)

One codon entry per line: `AminoAcid  Codon  Count  RSCU`.

```
Phe  UUU  791  1.27
Phe  UUC  450  0.73
Leu  UUA  754  1.96
...
```

### 4.3 General rules

| Field | Allowed values |
|-------|---------------------|
| Amino acid | `Phe`, `Leu`, ..., `TER` (case-insensitive; first-letter capital is enough) |
| Codon | `UUU`, `TTT`, ... (`T` is auto-normalised to `U`) |
| Count | any integer |
| RSCU | any float (e.g. `1.27`) |

Codons not listed in `STANDARD_ORDER` (e.g. non-standard codes) are
automatically appended at the end of the CSV -- no data is lost.

---

## 5. Output format

UTF-8 with BOM CSV (opens directly in Excel without garbled characters).
Example:

```csv
"Amino Acid","Codon","Species_A","Species_B","Species_C"
"Phe","UUU",1.27,1.31,1.25
"Phe","UUC",0.73,0.69,0.75
"Leu","UUA",1.96,2.00,1.92
"Leu","UUG",1.24,1.20,1.28
...
```

Columns 1 and 2 are fixed as `Amino Acid` and `Codon`; remaining columns
are the RSCU values of each species.

---

## 6. Command-line / programmatic usage

The GUI is just a wrapper around `merge_rscu_files`, which can be imported
and used in any script:

```python
from codonw2rscuplot import merge_rscu_files

n_species, n_codons, empty = merge_rscu_files(
    input_dir="my_rscu_results",
    output_csv="rscu_input_for_plot.csv",
    include_ter=True,
)
print(n_species, n_codons, empty)
```

---

## 7. FAQ

**Q: I get "No RSCU result files found in the selected folder".**
A: The tool recognises four extensions: `.blk`, `.out`, `.txt`, `.csv`. Make
sure that:
- The files really live in that folder.
- Their extension is one of the four.
- They are not empty (>= 1 byte and contain at least one recognisable
  codon entry).

**Q: One species column is entirely empty.**
A: That usually means the corresponding file is empty, has an unusual
encoding, or contains no recognisable codon entries. The species name
will be listed in the warning log -- regenerate the RSCU result for that
species.

**Q: Should I include stop codons?**
A: Yes -- the RSCU-Plot app explicitly asks for stop-codon rows. Only
uncheck *Include stop codons (TER)* when working strictly on coding
regions.

**Q: Can I customise the codon order?**
A: The order is currently fixed to the standard table. If you really need
a custom order, edit the `STANDARD_ORDER` constant near the top of the
source file.

**Q: Can this tool aggregate the per-gene RSCU rows inside a codonW `.out`?**
A: No. This tool assumes each input file already represents the RSCU of
**one species**. Please use codonW's *Codon usage* / *Correspondence
analysis* modules to first produce a per-species summary, and only then
let this tool merge them horizontally.

---

## 8. Credits

- Visualisation is provided by
  [**Codon Usage Visualization (RSCU-Plot)**](https://pcg-lab.shinyapps.io/RSCU-Plot/)
  (developed by Aboozar Soorni, pcg-lab.iut.ac.ir).
  If you use RSCU-Plot, please cite the relevant papers as required by
  the app:
  - <https://doi.org/10.1038/s41598-025-29339-x>
  - <https://doi.org/10.1007/s00425-025-04875-y>
  - <https://doi.org/10.1186/s12864-025-11729-0>
  - <https://doi.org/10.1007/s11033-025-10789-5>
- RSCU computation can be done with
  [**codonW**](http://codonw.sourceforge.net/) or
  [**MEGA**](https://www.megasoftware.net/).

---

## 9. Layout

```
codonw2rscuplot/
├── codonw2rscuplot.py   # main program
└── README.md            # this file
```