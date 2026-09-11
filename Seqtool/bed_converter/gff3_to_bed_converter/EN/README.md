# GFF3 to BED Converter

A graphical tool for converting GFF3 format to BED format.

## Features

- Graphical user interface, easy to use
- Select multiple feature types for conversion (gene, mRNA, exon, CDS)
- Real-time conversion log display
- Auto-generated BED file
- Preview output before saving
- Ensures consistent gene names across all features

## File Description

- `gff3_to_bed_gui.py` - Main program file
- `start_converter.bat` - Windows startup script
- `requirements.txt` - Dependency list (uses Python standard library mainly)
- `README.md` - Detailed documentation
- `快速开始.md` - Quick start guide (Chinese)

## Usage

### Method 1: Use Startup Script (Recommended)

Double-click the `start_converter.bat` file to launch the program.

### Method 2: Use Python Command

Open command line, navigate to the program directory, and run:

```bash
python gff3_to_bed_gui.py
```

## Steps

1. Click "Browse" to select GFF3 input file
2. The program will automatically parse the file and display detected gene information
3. The program will automatically set BED output file path (can also be modified manually)
4. Check the feature types to convert
5. (Optional) Click "Preview Output" to check conversion results
6. Click "Start Conversion" or save from preview window

## Key Features

- **Guaranteed Gene Name Consistency**: Ensures all features of the same gene use the same gene name
- **Smart Gene Recognition**: Automatically extracts gene IDs from GFF3 attributes field
- **Preview Function**: Preview output content before saving to avoid unexpected results
- **Real-time Log**: Displays conversion progress and statistics

## Format Description

### GFF3 Format

GFF3 files are tab-delimited text files with 9 columns:

1. seqid - Sequence identifier
2. source - Source
3. type - Feature type (gene, mRNA, exon, CDS, etc.)
4. start - Start position (1-based)
5. end - End position
6. score - Score
7. strand - Strand direction (+/-.)
8. phase - Phase
9. attributes - Attribute information

### BED Format

BED files are tab-delimited text files. This program generates BED format with 4-5 columns:

1. geneID/transcriptID - Gene ID or transcript ID (format: geneID or geneID/transcriptID)
2. start - Start position (1-based, keeping GFF3 original coordinates)
3. end - End position
4. featureType - Feature type (gene, mRNA, exon, CDS, etc.)
5. phase - Phase (optional column, only CDS features have values)

**Note**: For gene features, only geneID is shown. For other feature types (mRNA, exon, CDS, etc.), geneID is shown.

## Important Notes

- BED output uses 1-based coordinate system, keeping GFF3 original coordinates
- Ensure GFF3 file format is correct
- BED output file will be overwritten, please choose output path carefully
- gene type features only show geneID
- Other feature types (mRNA, exon, CDS, etc.) use geneID format
- phase column only has actual values for CDS features, others are empty

## Example

Example files are provided in the program directory:
- `基因成员.gff3` - GFF3 format example input
- `BED示例.txt` - BED format example output

## System Requirements

- Python 3.6 or higher
- Windows operating system
- No additional dependencies required

## License

This tool is free software and can be freely used and modified.
