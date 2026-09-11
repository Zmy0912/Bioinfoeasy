# BED/GFF3 Coordinate Converter (gene_to_0) - Project Structure

## Overview

This tool (directory `gene_to_0`) is a graphical tool for converting BED and GFF3 file coordinates by normalizing gene coordinates to start from 0.

## Project Structure

```
Seqtool/bed_converter/gene_to_0/EN/
├── bed_converter.py          # Main program file
├── README.md                 # Project overview (English)
├── QUICK_START.md            # Quick start guide (5-minute setup)
├── USER_GUIDE.md             # Detailed user manual
├── FAQ.md                    # Frequently asked questions
├── CHANGELOG.md              # Version history and updates
├── PROJECT_STRUCTURE.md      # This file
└── examples/
    ├── README.md             # Example files guide
    ├── example.bed           # BED format example
    └── example.gff3         # GFF3 format example
```

## File Descriptions

### Core Program
- **bed_converter.py**: Main Python application with Tkinter GUI

### Documentation
- **README.md**: Project introduction, features, installation, and basic usage
- **QUICK_START.md**: Step-by-step guide for quick start
- **USER_GUIDE.md**: Comprehensive user manual with detailed explanations
- **FAQ.md**: Solutions to common problems and issues
- **CHANGELOG.md**: Version history, new features, and improvements
- **PROJECT_STRUCTURE.md**: Project structure overview

### Examples
- **examples/README.md**: Guide for using example files
- **examples/example.bed**: Sample BED file for testing
- **examples/example.gff3**: Sample GFF3 file for testing

## Features

### Supported Formats
- **BED Format**: 0-based coordinate system
- **GFF3 Format**: 1-based coordinate system (maintains standard)

### Key Capabilities
- Graphical user interface
- File preview before conversion
- Coordinate normalization
- Preserve comments and attributes
- Statistics display
- Format validation

## Requirements

- Python 3.6 or higher
- tkinter (usually included with Python)

## Installation

No installation required. Simply download and run:

```bash
python bed_converter.py
```

## Usage Summary

1. Select file type (BED/GFF3)
2. Browse and select input file
3. Load and preview the file
4. Convert coordinates
5. Preview results
6. Save converted file

## File Formats

### BED Format
```
geneID  start  end  featureType  [phase]
```
- 0-based coordinates
- Columns: Gene ID, Start, End, Feature Type, Phase (optional)

### GFF3 Format
```
seqid  source  type  start  end  score  strand  phase  attributes
```
- 1-based coordinates
- Attributes must include `gene_id`
- Standard GFF3 format

## Documentation Quick Links

- **Quick Start**: See `QUICK_START.md` for immediate setup
- **Detailed Guide**: See `USER_GUIDE.md` for complete instructions
- **Common Issues**: See `FAQ.md` for troubleshooting
- **Version Info**: See `CHANGELOG.md` for updates
- **Examples**: See `examples/README.md` for test files

## Development

### Version
- Current: v1.1.0
- Release Date: 2026-03-08

### Language
- Python 3.x
- English interface and documentation

### License
Provided as-is for educational and research purposes.

## Support

For questions or issues:
1. Check `FAQ.md` for common problems
2. Review `USER_GUIDE.md` for detailed instructions
3. Test with `examples/` files
4. Verify input file format

## Contributing

This is a standalone tool. Contributions welcome:
- Bug fixes
- Feature suggestions
- Documentation improvements
- Format support additions

## Acknowledgments

Designed for genomic coordinate normalization tasks in bioinformatics research.
