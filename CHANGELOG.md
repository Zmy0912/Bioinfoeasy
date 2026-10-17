# Changelog

All notable changes to **Bioinfoeasy** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.3] - 2026-10-17

### Added
- `analyzer/KaKs_Analyzer/` — GUI tool for filtering Ka/Ks results and
  imputing missing values. Auto-detects both the organized `.xlsx` sheet and
  the raw `.kaks` output of KaKs Calculator, groups rows by gene and drops
  genes whose missing-value ratio exceeds a user-selected threshold.
- `analyzer/ENC_GC3s_Analyzer/` — quantifies the proportion of genes falling
  below, on or above Wright's (1990) neutral expectation curve in an ENC-GC3s
  plot, across a whole folder of CodonW `.out` files.
- `analyzer/GC_ENC_Analyzer/` — batch codon usage bias statistics
  (GC1/GC2/GC3/GC3s/GC12, with ENC computed exactly as in CodonW 1.4.4) plus
  GC12–GC3 regression, automatic gene-name recognition and background
  threading with a progress bar.
- `analyzer/SSR_Region_Classifier/` — classifies every SSR into CDS / Intron /
  IGS from a MISA result file together with the matching GenBank annotations.
- `analyzer/misa_tool/misa_ssr_viewer/` — turns a `.misa` file into
  publication-style heatmap, pie-chart and stacked-bar figures, with
  one-click export of the underlying tables.
- `analyzer/RepeatAnalyzer/` — counts and summarizes direct and palindromic
  repeats from Vmatch `repfind` output; standard library only.
- `Plots/Codon_plot/coa_plotter/` — supplementary figures for CodonW
  correspondence analysis, including correlation tests against nucleotide
  composition parameters (GC3s, GC, GC12).
- `Plots/Codon_plot/coa_simple_plotter/` — streamlined COA variant that keeps
  only the orthogonal-axis main panel while retaining batch input.
- `Plots/Codon_plot/GC12_GC3/gc12_gc3_regression_table/` — GC12–GC3 regression
  statistics for one or many species as a sortable table with CSV export.
- `Seqtool/fasta_header_renamer/` — batch FASTA header ⇄ filename renaming,
  with recursive scanning and a multi-sequence editor.
- `other_tools/codonw2rscuplot/` — aggregates per-species CodonW/MEGA RSCU
  output into the wide-table CSV expected by the RSCU-Plot web app; zero
  third-party dependencies.

### Changed
- `analyzer/misa_tool/`: `misa_parser` moved into its own sub-directory so
  that it sits alongside the new `misa_ssr_viewer`.
- Root `README.md`: refreshed the feature list, repository structure, module
  tables, usage examples and citation metadata for v0.1.3.
- Updated the module index READMEs under `analyzer/`, `Plots/`,
  `Plots/Codon_plot/`, `Seqtool/` and `other_tools/` to cover the new tools.
- Clarified in `README.md` that the tools introduced in this release keep
  their Chinese documentation under a `中文/` directory, while the `CN`/`EN`
  convention established in v0.1.2 is retained for all pre-existing modules.

## [0.1.2] - 2026-09-12

### Added
- Root-level metadata: `.gitignore`, `.gitattributes`, `.editorconfig`,
  `requirements.txt` (aggregated dependencies) and this `CHANGELOG.md`.
- Index README files for every top-level module (`analyzer/`, `Plots/`,
  `Seqtool/`, `Seqinfo/`, `other_tools/`) and their sub-modules.
- `archive/` directory to hold superseded tool versions, with its own README.
- GitHub community files: pull request template and issue template config.

### Changed
- **Repository-wide structure normalization**
  - Unified all version directories to `vMAJOR.MINOR.PATCH`
    (e.g. `v1.0` → `v1.0.0`, `0.1.0` → `v0.1.0`).
  - Unified all language directories to `CN` / `EN`
    (replaced the legacy `中文` directory names).
  - Standardized directory names to lowercase `snake_case`
    (`GC_ANALYZER` → `gc_analyzer`, `GC3ENCanalyzer` → `gc3enc_analyzer`,
    `Quadrants_analyzer` → `quadrants_analyzer`, `MISAtool` → `misa_tool`,
    `PiPlot` → `pi_plot`, `HeatmapGUI` → `heatmap_gui`,
    `Time Scale Generator` → `time_scale_generator`,
    `BED-Converter` → `bed_converter`).
  - Fixed the misspelled `Seqtool/comparision` → `Seqtool/comparison`.
  - Merged the duplicated `Seqtool/Seqtool_CN` module into
    `Seqtool/other_tools/CN` and moved the English scripts to
    `Seqtool/other_tools/EN`.
  - Split the mixed `cp_genome_comparison/` directory into clean
    `EN/` and `CN/` sub-directories.
  - Flattened redundant nested directories in `batch_renamer_tool` and
    `time_scale_generator`, and removed an empty `GC12_GC3/v3.0` directory.
- Consolidated all superseded versions into `archive/codon_plot_legacy/`.
- Corrected the repository identity to **Bioinfoeasy** across
  `README.md`, `CITATION.cff` and `CONTRIBUTING.md`.
- Fixed invalidation in `CITATION.cff` (`cff-version` was incorrectly set to
  the software version) and synchronized the version and DOI.

## [0.1.1] - 2026-09-11

### Changed
- Major update to the codon-usage plotting module (`Plots/Codon_plot/`),
  including the GC12-GC3, GC3-ENC, PR2-Plot and codon-usage plotters.

## [0.1.0] - 2026-04-01

### Added
- First public release of Bioinfoeasy.
- `analyzer/` – chloroplast genome and repeat analysis
  (CodonW, GC content, GC12-GC3 regression, GC3-ENC, RSCU Kruskal-Wallis,
  quadrants analysis, MISA SSR parser, sequence repeat analyzer).
- `Plots/` – visualization tools (codon plots, heatmaps and heatmap GUI,
  Pi plots, time scale generator).
- `Seqtool/` – sequence processing and format conversion
  (FASTA validator/joiner, BED converters, GenBank→mVISTA converter,
  FASTA comparison and miscellaneous sequence utilities).
- `Seqinfo/` – genome annotation file analyzer.
- `other_tools/` – batch renamer and file lister.
- Project documentation: `README.md`, `LICENSE` (MIT), `CITATION.cff`,
  `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `SECURITY.md`.

[Unreleased]: https://github.com/Zmy0912/Bioinfoeasy/compare/v0.1.3...HEAD
[0.1.3]: https://github.com/Zmy0912/Bioinfoeasy/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/Zmy0912/Bioinfoeasy/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/Zmy0912/Bioinfoeasy/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/Zmy0912/Bioinfoeasy/releases/tag/v0.1.0
