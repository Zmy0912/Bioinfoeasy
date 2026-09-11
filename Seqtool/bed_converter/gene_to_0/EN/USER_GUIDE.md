# BED/GFF3 Coordinate Converter - User Guide

## Features

This tool supports coordinate conversion for two file formats:
- **BED Format**: Converts gene internal coordinates in BED files to start from 0
- **GFF3 Format**: Converts gene internal coordinates in GFF3 files to start from 0

The conversion maintains the relative distances between features within each gene.

## Usage Steps

### General Steps

1. **Launch the Program**
   - Double-click `bed_converter.py`
   - Or run from command line: `python bed_converter.py`

2. **Select File Format**
   - Choose "BED Format" or "GFF3 Format" radio button

3. **Select Input File**
   - Click "Browse" button to select the file to convert
   - Or enter the file path directly in the input box

4. **Load File**
   - Click "Load File" button to view the original file content
   - The file content will be displayed in the preview area

5. **Convert Coordinates**
   - Click "Convert" button to start conversion
   - The program will automatically identify gene IDs and convert coordinates
   - Conversion results will be displayed in the preview area in real-time

6. **Save File**
   - Click "Save File" button to save the converted results
   - Choose the save location and filename
   - BED format default filename is "converted.bed"
   - GFF3 format default filename is "converted.gff3"

7. **Clear**
   - Click "Clear" button to clear all data and start over

## BED File Format

Input BED file format:
```
geneID  start  end  featureType  phase(optional)
```

- Column 1: Gene/transcript ID
- Column 2: Start position (0-based)
- Column 3: End position (0-based)
- Column 4: Feature type (exon, CDS, gene, mRNA, etc.)
- Column 5: Phase (optional: 0, 1, 2, or .)

## GFF3 File Format

Input GFF3 file format (standard GFF3):
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

The program extracts `gene_id` from the attributes field to identify genes.

## Conversion Logic

For each gene:
1. Find the minimum start position of the gene (gene_start)
2. Subtract gene_start from all feature coordinates within the gene
3. Maintain relative distances between features

### BED Format Example

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

### GFF3 Format Example

**Before:**
```
chr1  Ensembl  gene  61386544  61389501  .  +  .  gene_id=AoNCED1
chr1  Ensembl  exon  61388988  61389501  .  +  .  gene_id=AoNCED1
chr1  Ensembl  exon  61387536  61387539  .  +  .  gene_id=AoNCED1
```

**After (maintaining GFF3 format):**
```
chr1  Ensembl  gene  1  2958  .  +  .  gene_id=AoNCED1
chr1  Ensembl  exon  2445  2958  .  +  .  gene_id=AoNCED1
chr1  Ensembl  exon  993  996  .  +  .  gene_id=AoNCED1
```

## Notes

### General Notes
- The program preserves comment lines (lines starting with #) from the original file
- All features belonging to the same gene are processed together
- Preview the original file before conversion to ensure correct format
- It is recommended to preview conversion results before saving

### BED Format Notes
- BED format uses 0-based coordinate system
- Converted coordinates start directly from 0

### GFF3 Format Notes
- GFF3 uses 1-based coordinate system
- The converted file maintains the GFF3 format and 1-based coordinate system
- Ensure the GFF3 file contains `gene_id` property in the attributes field
- The program preserves all original GFF3 file structure, comment lines, and attributes

## Runtime Environment

- Python 3.x
- Requires tkinter (usually installed with Python)

## Troubleshooting

**Q: Why do GFF3 converted coordinates still start from 1?**
A: This is normal. The GFF3 standard uses a 1-based coordinate system. The program uses 0-based for internal calculations but converts to standard 1-based format for output to ensure compatibility.

**Q: My GFF3 file shows no data after conversion?**
A: Please check if the file contains `gene_id` in the attributes field. The program identifies genes through gene_id.

**Q: Can I process multiple files at once?**
A: Currently, the program processes one file at a time. Click "Clear" between files.

For other issues, please check if the input file format meets the requirements.
