# Repeat Analyzer

A GUI tool for counting and summarizing **direct (forward) repeats** and **palindromic (reverse) repeats** from Vmatch `repfind` output files.

## Overview

- Each input folder contains one plain-text result file per **species**; the file name (extension and trailing `repeats` removed) is used as the species name (e.g. `NC01repeats.txt` → `NC01`).
- The input format is the standard `repfind` output, one repeat per line:

  ```
  length1  pos1  P/D  length2  pos2  distance  E-value
  ```

  where **P = Palindromic repeat** and **D = Direct repeat**. Lines starting with `#` are ignored.

- The tool parses every `.txt` file in the selected folder and displays per-species statistics in a table, with an optional CSV export.

## Requirements

- Python 3.8+ (tested with Python 3.13)
- No third-party packages required — the tool uses only the Python standard library (Tkinter is bundled with the official Python installer on Windows).

## Usage

Run the program:

```
python repeat_analyzer.py
```

Then:

1. Click **Select Folder** and choose the folder containing your repeat result files.
2. Optionally set the filters:
   - **Min length** — only count repeats with an arm length ≥ this value (default 10, matching Vmatch's `-l 10`).
   - **Max E-value** — only count repeats with an E-value ≤ this value (leave blank for no limit, e.g. `1e-5`).
3. Click **Analyze**. The table shows one row per species plus a highlighted **TOTAL** summary row.
4. Click any column header to sort the table (ascending/descending toggle).
5. Click **Export CSV** to save the table. The CSV is written with a UTF-8 BOM, so it opens correctly in Microsoft Excel.

## Output Columns

| Column | Description |
|---|---|
| Species | Species name derived from the file name |
| Direct (D) Repeats | Number of forward repeats |
| Palindromic (P) Repeats | Number of reverse/palindromic repeats |
| Total | Total number of repeats |
| Total Length (bp) | Sum of repeat arm lengths (bp) |
| Mean Length (bp) | Average repeat length |
| Min Length (bp) | Shortest repeat |
| Max Length (bp) | Longest repeat |

A **TOTAL** row at the bottom aggregates all species.

## Notes

- Only `.txt` files in the selected folder are parsed.
- Malformed lines (not matching the 7-field `repfind` layout) are skipped automatically.
- The default example dataset was generated with `-best 50`, so each species shows up to 50 repeats; adjust the filters to subset them.
