# BED/GFF3 Coordinate Converter

A tool to convert BED and GFF3 file coordinates by normalizing gene coordinates to start from 0.

## Features

- **Dual Format Support**: Process both BED and GFF3 file formats
- **Coordinate Normalization**: Convert gene internal coordinates to start from 0
- **Preview Function**: View conversion results before saving
- **GUI Interface**: User-friendly graphical interface
- **Format Preservation**: Maintain file structure and comments

## Requirements

- Python 3.6 or higher
- tkinter (usually included with Python)

## Installation

No installation required. Simply download the `bed_converter.py` file and run it.

```bash
python bed_converter.py
```

## Usage

1. **Launch the Program**
   - Double-click `bed_converter.py`
   - Or run from command line: `python bed_converter.py`

2. **Select File Format**
   - Choose "BED Format" or "GFF3 Format" radio button

3. **Load Input File**
   - Click "Browse" to select your input file
   - Or enter the file path directly
   - Click "Load File" to preview the original content

4. **Convert Coordinates**
   - Click "Convert" to process the file
   - View the results in the preview area
   - Check the statistics for gene and feature counts

5. **Save Results**
   - Click "Save File" to save the converted data
   - Choose the destination and filename
   - BED files default to `.bed` extension
   - GFF3 files default to `.gff3` extension

6. **Clear and Restart**
   - Click "Clear" to reset and process another file

## File Formats

### BED Format

```
gene_id  start  end  feature_type  phase(optional)
```

- Column 1: Gene/transcript ID
- Column 2: Start position (0-based)
- Column 3: End position (0-based)
- Column 4: Feature type (exon, CDS, gene, mRNA, etc.)
- Column 5: Phase (optional: 0, 1, 2, or .)

### GFF3 Format

```
seqid  source  type  start  end  score  strand  phase  attributes
```

- Column 1: Sequence ID
- Column 2: Source
- Column 3: Feature type (gene, mRNA, exon, CDS, etc.)
- Column 4: Start position (1-based)
- Column 5: End position (1-based)
- Column 6: Score
- Column 7: Strand (+, -, or .)
- Column 8: Phase (0, 1, 2, or .)
- Column 9: Attributes (key=value pairs separated by ;)

## Conversion Logic

For each gene:

1. Identify the gene's minimum start position (gene_start)
2. Subtract gene_start from all feature coordinates within the gene
3. Preserve relative distances between features

### BED Example

**Before:**
```
AoNCED1  61386543  61389501  gene
AoNCED1  61388987  61389501  exon
AoNCED1  61387535  61387538  exon
```

**After:**
```
AoNCED1  0  2958  gene
AoNCED1  2444  2958  exon
AoNCED1  992  995  exon
```

### GFF3 Example

**Before:**
```
chr1  Ensembl  gene  61386544  61389501  .  +  .  gene_id=AoNCED1
chr1  Ensembl  exon  61388988  61389501  .  +  .  gene_id=AoNCED1
```

**After:**
```
chr1  Ensembl  gene  1  2958  .  +  .  gene_id=AoNCED1
chr1  Ensembl  exon  2445  2958  .  +  .  gene_id=AoNCED1
```

## Notes

- Comment lines (starting with #) are preserved
- For GFF3 files, ensure attributes contain `gene_id` property
- The tool maintains the original file structure and all attributes
- Always preview results before saving
- GFF3 files use 1-based coordinate system (standard GFF3)

## Troubleshooting

**Q: GFF3 conversion produces no data**
A: Ensure your GFF3 file contains `gene_id` in the attributes field.

**Q: Conversion fails with error**
A: Check that the input file format matches the selected type (BED/GFF3).

**Q: Can I process multiple files at once?**
A: Currently, the tool processes one file at a time. Click "Clear" between files.

## License

This tool is provided as-is for educational and research purposes.

## Version History

- v1.1: Added GFF3 format support
- v1.0: Initial release with BED format support

## Contact

For issues or questions, please refer to the usage documentation or check your input file format.
