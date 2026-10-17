# Other Tools

General-purpose utility tools that are not specific to sequence analysis.

## Contents

| Module | Description |
|--------|-------------|
| [`batch_renamer_tool/`](batch_renamer_tool/) | Batch file renaming with pattern matching (CN/EN) |
| [`codonw2rscuplot/`](codonw2rscuplot/) | Aggregate per-species CodonW/MEGA RSCU output into the wide-table CSV expected by the [RSCU-Plot](https://pcg-lab.shinyapps.io/RSCU-Plot/) web app; zero third-party dependencies |
| [`file_lister/`](file_lister/) | Browser-based folder file listing with sort, search, copy and CSV export (CN/EN) |

## Usage

```bash
# Batch renamer (English)
cd batch_renamer_tool/EN
python batch_renamer.py

# Build an RSCU-Plot input table from CodonW output
cd codonw2rscuplot
python codonw2rscuplot.py
```

For the file lister, open `file_lister/file_lister.html` (English) or
`file_lister/文件列表查看器.html` (Chinese) directly in a browser — no Python
required.

Back to the [project root](../README.md).
