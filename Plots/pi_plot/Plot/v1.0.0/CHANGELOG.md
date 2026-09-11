# Changelog

This document records all important changes to the Pi Sliding Window Gene Labeling Program.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned Features
- Support for batch analysis of multiple samples
- Add interactive chart options
- Support CSV format for Pi files
- Add statistical test options

## [1.0.0] - 2026-01-03

### Added
- Complete Pi sliding window analysis functionality
- Automatic significant gene identification algorithm
- Two significance criteria:
  - Relative significance (based on Z-score)
  - Absolute significance (based on Pi threshold)
- Intelligent label layout system:
  - Automatic label overlap avoidance
  - Adjacent gene labels assigned different heights
  - Support for up to 50 iterations for optimization
- High-quality chart generation:
  - 300 DPI resolution
  - Customizable font size
  - Multiple colors for different types of significant genes
- Detailed analysis report generation:
  - PNG visualization chart
  - Text format detailed information file
- Flexible command-line parameter interface
- Complete documentation system:
  - README.md (project introduction)
  - USAGE_GUIDE.md (detailed usage guide)
  - CHANGELOG.md (update log)
  - requirements.txt (dependency list)

### Technical Features
- Pandas-based data processing
- Matplotlib visualization
- NumPy numerical computation support
- Chinese font display support
- Automatic Y-axis scaling to accommodate labels

### Gene Labeling Features
- Support for extracting gene information from GFF3 files
- Calculate gene Z-score relative to local region
- Sort genes by significance level
- Intelligent label position algorithm:
  - Group adjacent genes within 15000bp
  - Assign genes within groups to different height layers
  - Collision detection to avoid overlaps
- Label color coding:
  - Yellow: Relative significance only
  - Light blue: Absolute significance only
  - Orange: Both criteria satisfied

### Parameter System
- `--pi-file`: Pi data file path (required)
- `--gff3-file`: GFF3 gene annotation file path (required)
- `--output`: Output file name
- `--gene-fontsize`: Gene label font size
- `--max-genes`: Maximum number of genes to label
- `--window-size`: Local comparison window size
- `--significance-std`: Relative significance criterion
- `--absolute-threshold`: Absolute Pi threshold

### Documentation
- Detailed README.md documentation
- Comprehensive usage guide (USAGE_GUIDE.md)
- Dependency list (requirements.txt)
- Command-line help information
- Inline code documentation strings

## Version Information

### Version Number Format

Version numbers use the format `MAJOR.MINOR.PATCH`:
- **MAJOR**: Incompatible API changes
- **MINOR**: Backwards-compatible functionality additions
- **PATCH**: Backwards-compatible bug fixes

### Change Types

- **Added**: New features
- **Changed**: Changes in existing functionality
- **Deprecated**: Soon-to-be removed features
- **Removed**: Removed features
- **Fixed**: Bug fixes
- **Security**: Security-related fixes

---

[Unreleased]: https://github.com/MingyuanZhang/pi_plot/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/MingyuanZhang/pi_plot/releases/tag/v1.0.0
