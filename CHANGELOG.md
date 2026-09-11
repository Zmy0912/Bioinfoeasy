# Changelog

All notable changes to **Bioinfoeasy** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/Zmy0912/Bioinfoeasy/compare/v0.1.2...HEAD
[0.1.2]: https://github.com/Zmy0912/Bioinfoeasy/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/Zmy0912/Bioinfoeasy/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/Zmy0912/Bioinfoeasy/releases/tag/v0.1.0
