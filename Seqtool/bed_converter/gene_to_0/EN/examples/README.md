# Example Files Guide

This folder contains example files for testing the BED/GFF3 coordinate converter.

## File List

### 1. example.bed
- BED format example file
- Contains 3 genes (AoNCED1, AoNCED2, AoNCED3)
- Includes gene, mRNA, exon, CDS features
- Includes phase information
- Contains Chinese comments

### 2. example.gff3
- Standard GFF3 format example file
- Contains 3 genes (AoNCED1, AoNCED2, AoNCED3)
- Located on different chromosomes (chr1, chr2)
- Includes complete GFF3 header information
- Contains gene, mRNA, exon, CDS features
- Includes both positive and negative strand genes
- Attributes field contains gene_id

## Using Examples

### Testing BED Format Conversion

1. Launch `bed_converter.py`
2. Select the "BED Format" radio button
3. Browse and select `examples/example.bed`
4. Click "Load File" to view original content
5. Click "Convert" to perform coordinate conversion
6. View the results in the preview area
7. Check the statistics (gene count and feature count)
8. Click "Save File" to save results

### Testing GFF3 Format Conversion

1. Launch `bed_converter.py`
2. Select the "GFF3 Format" radio button
3. Browse and select `examples/example.gff3`
4. Click "Load File" to view original content
5. Click "Convert" to perform coordinate conversion
6. View the results in the preview area
7. Check the statistics (gene count and feature count)
8. Note that GFF3 format maintains 1-based coordinate system
9. Click "Save File" to save results

## Expected Results

### BED Format After Conversion

The first feature of each gene should start from coordinate 0. For example AoNCED1:

```bed
AoNCED1	0	2958	gene
AoNCED1	2444	2958	mRNA
AoNCED1	2444	2958	exon
AoNCED1	992	995	exon
```

### GFF3 Format After Conversion

The first feature of each gene should start from coordinate 1 (GFF3 standard). For example AoNCED1:

```gff3
chr1	Ensembl	gene	1	2958	.	+	.	ID=gene:AO001;gene_id=AoNCED1;Name=NCED1
chr1	Ensembl	mRNA	1	2958	.	+	.	ID=transcript:AO001.1;Parent=gene:AO001;gene_id=AoNCED1
chr1	Ensembl	exon	2445	2958	.	+	.	ID=exon:AO001.1.1;Parent=transcript:AO001.1;gene_id=AoNCED1
chr1	Ensembl	exon	993	996	.	+	.	ID=exon:AO001.1.2;Parent=transcript:AO001.1;gene_id=AoNCED1
```

## Verification Points

1. **Coordinate Range**: Check if coordinates within each gene start from 0 (BED) or 1 (GFF3)
2. **Relative Distances**: Verify that relative distances between features remain unchanged
3. **Feature Count**: Confirm that converted feature count matches the original file
4. **Attribute Preservation**: All attributes in GFF3 files should be fully preserved
5. **Comment Lines**: Comment lines should be fully preserved
6. **Statistics**: Gene count and feature count should display correctly

## Format Validation

### BED Format Validation
- At least 4 columns
- Coordinates are non-negative integers
- Features with the same gene ID are grouped together

### GFF3 Format Validation
- Contains 9 columns
- Start position < End position
- Attributes field contains gene_id
- Coordinates are positive integers (1-based)

## Troubleshooting

If conversion fails or results are incorrect:

1. Check if the file format meets requirements
2. Confirm you selected the correct file type (BED/GFF3)
3. Compare your file format with the BED example
4. Check if the GFF3 file contains the gene_id attribute
5. Refer to the FAQ for common issues

## Custom Examples

You can create your own test files based on these examples:

- Modify gene IDs
- Adjust coordinate ranges
- Add more feature types
- Update attribute fields

Ensure the file format remains consistent with the example files.
