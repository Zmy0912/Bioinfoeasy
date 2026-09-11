# Frequently Asked Questions (FAQ)

## General Issues

### 1. Program won't start

**Problem**: Double-clicking `bed_converter.py` does nothing

**Solutions**:
- Ensure Python 3.6 or higher is installed
- Run from command line: `python bed_converter.py` to see error messages
- Check if tkinter is installed: `python -m tkinter` (should open a test window)

### 2. Display issues in the interface

**Problem**: Garbled text or incomplete display in the interface

**Solutions**:
- Ensure the system uses UTF-8 encoding
- Try setting encoding in command line: `chcp 65001` (Windows)
- Use a font that supports Chinese characters

## BED Format Issues

### 3. Negative coordinates after BED conversion

**Problem**: Some coordinates become negative after conversion

**Reason**: For the same gene ID, some features have start positions less than the gene's minimum start position

**Solution**: Check if all features with the same gene ID in the original file actually belong to the same gene

### 4. Fewer lines after conversion

**Problem**: Output file has fewer lines than input file

**Reason**: Some lines don't meet BED format requirements and were skipped

**Solution**: Check the original file format and ensure each line has at least 4 columns

## GFF3 Format Issues

### 5. No data after GFF3 conversion

**Problem**: Preview area is empty after clicking convert

**Reason**: Program cannot extract gene_id from the file

**Solutions**:
- Ensure the GFF3 file's attributes field contains `gene_id=xxx`
- Check if the file format is standard GFF3 (tab-separated)
- Example: `gene_id=AoNCED1;transcript_id=xxx`

### 6. Some lines not converted in GFF3 file

**Problem**: Some data lines disappear after conversion

**Reason**: These lines don't have the gene_id attribute

**Solutions**:
- Ensure all lines that need conversion have the gene_id attribute
- Check if the attributes field format is correct (separated by semicolons)

### 7. GFF3 coordinates still 1-based after conversion

**Problem**: Want coordinates to start from 0, but results still start from 1

**Note**: This is normal. The GFF3 standard uses a 1-based coordinate system

**Explanation**:
- The program uses 0-based for internal calculations
- Output is converted to standard 1-based format for compatibility
- If you need 0-based coordinates, consider using BED format

## File Operation Issues

### 8. Save file failed

**Problem**: Error message when clicking save

**Possible Causes**:
- Destination path doesn't exist
- No write permission
- File is being used by another program

**Solutions**:
- Choose a different save location
- Close other programs that might be using the file
- Check file permissions

### 9. Loading large files is slow

**Problem**: Program responds slowly when loading large files

**Solutions**:
- Large files (>100MB) may take some time to load
- Be patient and don't click repeatedly
- You can test functionality with a small file first

## Performance Issues

### 10. Processing takes too long

**Problem**: It takes a long time to process large amounts of data

**Suggestions**:
- This is normal and depends on data volume
- You can test with a small file to verify functionality
- Ensure your computer has sufficient memory

## Format Issues

### 11. Which GFF3 feature types are supported?

**Supported Feature Types**:
- gene (genes)
- mRNA / transcript (transcripts)
- exon (exons)
- CDS (coding sequences)
- UTR (untranslated regions)
- Other standard GFF3 feature types

### 12. What are BED file format requirements?

**Required Columns**:
- Column 1: Gene ID
- Column 2: Start position (integer)
- Column 3: End position (integer)
- Column 4: Feature type

**Optional Columns**:
- Column 5: Phase information

## Error Messages

### 13. "Index out of range" error

**Reason**: Some lines have insufficient columns

**Solution**: Check the file format and ensure each line has at least the required number of columns

### 14. "Cannot convert to integer" error

**Reason**: Start or end positions are not valid numbers

**Solution**: Check if the coordinate columns contain pure numbers

## Other Issues

### 15. Can I batch process multiple files?

**Current Version**: Batch processing is not supported

**Suggestions**:
- Use scripts or write loops to process multiple files
- Or wait for batch processing feature in future versions

### 16. Can converted data be used in other software?

**Compatibility**:
- BED format: Compatible with most genomic analysis tools
- GFF3 format: Conforms to standard GFF3 specifications, highly compatible

**Recommendation**: Confirm the target software's format requirements before use

### 17. Does the program modify original files?

**Answer**: No

**Explanation**:
- The program only reads original files
- Conversion results are saved as new files
- Original files remain unchanged

### 18. How to verify conversion results are correct?

**Verification Methods**:
1. Check if the first feature of each gene starts from 0 (BED) or 1 (GFF3)
2. Compare feature distances in original and converted files
3. Manually calculate and verify with small samples
4. Check statistics for gene and feature counts

## Need More Help?

If you encounter issues not listed above:
1. Check if the input file format is correct
2. Review the user guide documentation
3. Try testing with example files
4. Check Python and library versions

## Tips

- It's recommended to back up original files
- Test with example files before processing important data
- Preview results before saving
- Be aware of differences between BED and GFF3 formats
