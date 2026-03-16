# Changelog

This document records all important changes to the Pi Sliding Window Gene Annotation Program.

Format based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
version numbers follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Planned Features
- Support batch analysis of multiple samples
- Add interactive chart options
- Support CSV format Pi files
- Add statistical test options

## [1.1.0] - 2026-01-04

### Added
- Added `--lower-threshold` parameter to support detection of lowest Pi value genes
- Support multiple significance type combinations (L+R+A, R+A, L+R, L+A, etc.)
- Added 7 significance type markers, covering all possible combinations

### Enhanced
- **Gene sorting logic refactor**: Prioritize display of genes meeting multiple conditions
  - L+R+A (3 conditions) > R+A, L+R, L+A (2 conditions) > R, A, L (1 condition)
  - Within same priority level, sort by Z-score from high to low
- **Parameter processing logic optimization**: Enable relative threshold mode only when `--significance-std` is explicitly specified
- **Color coding system update**: Assign unique colors for each significance type

### Fixed
- Fixed issue where `max_genes` limit was not applied in "lower threshold + absolute threshold" mode

### Documentation
- Updated README.md, added lower threshold parameter description
- Updated USAGE_GUIDE.md, added new usage scenarios and parameter combination descriptions
- Added gene sorting priority table
- Updated color coding descriptions

## [1.0.0] - 2026-01-03

### Added
- Complete Pi sliding window analysis functionality
- Automatic algorithm for identifying significant genes
- Two significance criteria:
  - Relative significance (based on Z-score)
  - Absolute significance (based on Pi threshold)
- Smart label layout system:
  - Automatically avoid label overlaps
  - Assign different heights to adjacent gene labels
  - Support up to 50 iterations of optimization
- High-quality chart generation:
  - 300 DPI resolution
  - Customizable font size
  - Multiple colors to annotate different types of significant genes
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
- Data processing based on pandas
- Visualization using matplotlib
- NumPy numerical calculation support
- Support Chinese font display
- Automatic Y-axis scaling to accommodate labels

### Gene Annotation Features
- Support extracting gene information from GFF3 files
- Calculate Z-score of gene relative to local region
- Sort genes by significance level
- Smart label position algorithm:
  - Group adjacent genes within 15000bp
  - Assign genes in group to different height layers
  - Collision detection to avoid overlaps
- Label color coding (v1.1.0 and later):
  - Yellow: Relative only [R]
  - Light blue: Absolute only [A]
  - Light green: Lower only [L]
  - Orange: Meet both relative and absolute [R+A]
  - Light yellow: Meet lower and relative [L+R]
  - Light cyan: Meet lower and absolute [L+A]
  - Light green: Meet all three [L+R+A]

### Parameter System
- `--pi-file`: Pi data file path (required)
- `--gff3-file`: GFF3 gene annotation file path (required)
- `--output`: Output file name
- `--gene-fontsize`: Gene annotation font size
- `--max-genes`: Maximum number of genes to annotate
- `--window-size`: Local comparison window size
- `--significance-std`: Relative significance standard
- `--absolute-threshold`: Absolute Pi threshold
- `--lower-threshold`: Lower Pi threshold (added in v1.1.0)

### Documentation
- Detailed README.md documentation
- Comprehensive usage guide (USAGE_GUIDE.md)
- Dependency list (requirements.txt)
- Command-line help information
- Inline code documentation strings

## Version Notes

### Version Number Format

Version numbers use `MAJOR.MINOR.PATCH` format:
- **MAJOR**: Incompatible API changes
- **MINOR**: Backwards-compatible functionality additions
- **PATCH**: Backwards-compatible bug fixes

### Change Types

- **Added**: New features
- **Changed**: Changes to existing functionality
- **Deprecated**: Soon-to-be removed features
- **Removed**: Removed features
- **Fixed**: Bug fixes
- **Security**: Security-related fixes

---

## Version History
[1.1.0]: https://github.com/Zmy0912/Bioinfo/tree/master/PiPlot/v1.1.0
[1.0.0]: https://github.com/Zmy0912/Bioinfo/tree/master/PiPlot/v1.0.0

---
END
