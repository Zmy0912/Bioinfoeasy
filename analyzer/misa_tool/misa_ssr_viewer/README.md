# MISA SSR Analyzer

A cross-platform desktop GUI application for visualizing the results of
[MISA](https://webblast.ipk-gatersleben.de/misa/) (MIcroSAtellite) searches.
It turns a `.misa` result file into three publication-style figures and also
exports the underlying data tables for use in GraphPad Prism, Excel, TBtools,
R, Origin, etc.

![Heatmap](_smoke_heatmap.png)
![Pie chart](_smoke_pie.png)
![Stacked bar](_smoke_stacked.png)

---

## Features

| Panel | What it shows |
|------:|---------------|
| **B. Heatmap** | SSR motif presence (blue = 1) / absence (white = 0) across all species. Motifs are ordered by repeat-unit length (1 - 10+ bp) and grouped in the upper subplot. Species are ordered by **Ward hierarchical clustering** on the binary motif matrix. |
| **C. Pie chart** | Percentage distribution of SSR categories (Mono-, Di-, Tri-, Tetra-, Penta-, Hexa-nucleotide, Compound *c* / *c**) across all species, with leader-line labels to avoid overlap. |
| **D. Stacked bar** | Per-species SSR category distribution and abundance. Species are sorted by total SSR count (descending) so the most SSR-rich species appear first. |

Additional capabilities:

- One-click export of all underlying data tables (tab-separated, UTF-8-BOM)
  so they open cleanly in Excel and other tools.
- Save individual or all figures as PNG, PDF or SVG (300 dpi by default).
- Embedded matplotlib navigation toolbar (zoom, pan, save, etc.).

---

## Requirements

- Python **≥ 3.8**
- `tkinter` (bundled with most official CPython installers; on Linux you may
  need `sudo apt install python3-tk`)
- `matplotlib` **≥ 3.5**
- `numpy`

Install the Python dependencies with:

```bash
pip install matplotlib numpy
```

---

## Run

```bash
python misa_ssr_viewer.py
```

A window opens with three tabs (Heatmap, Pie chart, Stacked bar). Use the
top toolbar to pick a `.misa` file, choose an output directory, and press
**Run Analysis**.

### Input file

A MISA result file (default extension `.misa`) with a tab-separated header
`ID   SSR nr.   SSR type   SSR   size   start   end` and one SSR per row.

### Exported data tables

When you press **Export Tables** (or `File → Export data tables`), the
following files are written to the chosen output directory:

| File | Contents |
|------|----------|
| `B1_heatmap_presence_absence_matrix.txt` | Species × motif binary matrix, species in **Ward cluster order**, motifs in (length, alphabetical) order. |
| `B1b_heatmap_matrix_original_species_order.txt` | Same matrix but species in their original numeric order (useful for software with its own clustering). |
| `B2_motif_info.txt` | Motif order, repeat-unit length, species count, frequency (%). |
| `B3_motif_length_group_count.txt` | Number of motifs per repeat-unit-length group (used for the upper subplot). |
| `B4_species_ward_cluster_order.txt` | Final Ward leaf order with each species’ total SSR count and distinct motif count. |
| `C_pie_SSR_category_percentage.txt` | All eight MISA categories (`p1`–`p6`, `c`, `c*`) with counts and percentages. |
| `C_pie_SSR_category_percentage_merged.txt` | Same data but with `c` and `c*` merged into one *Compound* slice. |
| `D1_stackedbar_species_category_counts.txt` | Wide format: species (descending total) × category counts + total. |
| `D2_stackedbar_long_format.txt` | Long format (Species, MISA_Type, Category, Count, SpeciesTotal) suitable for ggplot2 / Origin. |

---

## Workflow

1. **Run MISA** on your FASTA file:
   ```
   misa.pl input.fasta output.misa
   ```
2. Launch the analyzer:
   ```
   python misa_ssr_viewer.py
   ```
3. Click **Browse...** next to *MISA file* and select `output.misa`.
4. (Optional) Click **Browse...** next to *Output dir* and pick a folder for
   the exported tables and figures.
5. Press **Run Analysis** — the three tabs populate immediately.
6. Press **Export Tables** to write the eight `.txt` files.
7. Use **Plots → Save current figure as PNG...** (or *Save all figures as
   PNG...*) to export the figures as 300 dpi images.

---

## Method notes

- **Motif extraction** uses the regular expression
  `\(([ACGTN]+)\)(\d+)`, which captures every repeat unit inside a (compound)
  SSR string. Complex compound SSRs (`c*`) whose interruption contains
  degeneracies such as `(A<TT>TTT...)` contribute only their well-formed
  motifs; this matches how MISA reports them.
- **Ward clustering** is computed on the binary motif matrix using squared
  Euclidean distance (= number of differing motifs) and Ward’s
  minimum-variance linkage, implemented from scratch with the Lance–Williams
  update formula. No `scipy` is required.
- **Heatmap colours**: a discrete two-colour map (`white` = 0, `#1f77b4` = 1)
  as described in the figure caption. Columns are lightly shaded by their
  repeat-unit-length group for visual grouping.
- **Pie chart labels** are placed outside the pie with thin leader lines and
  a minimum vertical spacing of 0.20 (data units, radius = 1) so even the
  smallest slices do not overlap.

---

## Project layout

```
.
├── README.md
├── misa_ssr_viewer.py        # the GUI application
├── process_misa.py           # original CLI batch script (still usable)
├── 01_20.fasta.misa          # example input
└── SSR_plot_data/            # tables produced by process_misa.py
```

`process_misa.py` is kept for batch / scripted use; the GUI version is
recommended for interactive exploration.

---

## Licence & citation

If you use this tool in a publication, please cite the original MISA paper:

> Beier S., Thiel T., Münch T., Scholz U., Mascher M. (2017).
> MISA-web: a web server for microsatellite prediction.
> *Bioinformatics* 33(16): 2583–2585. https://doi.org/10.1093/bioinformatics/btx198