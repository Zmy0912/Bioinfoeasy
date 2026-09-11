# Quick Start Guide

## Get Started in 5 Minutes

### Step 1: Launch the Program

Double-click `bed_converter.py`, or enter the following in the command line:

```bash
python bed_converter.py
```

### Step 2: Select File Format

- **BED Format**: If your file format is `gene_id  start  end  type  phase`, select "BED Format"
- **GFF3 Format**: If your file is in standard GFF3 format, select "GFF3 Format"

### Step 3: Open File

1. Click the "Browse" button
2. Select your file
3. Click "Load File" to view the content

### Step 4: Convert Coordinates

1. Click the "Convert" button
2. Wait for the processing to complete
3. View the results in the preview area

### Step 5: Save Results

1. Verify the preview results are correct
2. Click "Save File"
3. Choose the save location and filename
4. Conversion complete!

## Test with Examples

Want to experience the functionality quickly? Use the example files:

### BED Example Test
1. Select "BED Format"
2. Open `examples/example.bed`
3. Click convert, view results

### GFF3 Example Test
1. Select "GFF3 Format"
2. Open `examples/example.gff3`
3. Click convert, view results

## Key Tips

- ✅ Click "Load File" to preview before conversion
- ✅ Click "Save File" after verifying preview results
- ✅ Original files are never modified
- ✅ Comment lines are automatically preserved
- ✅ Relative distances between features remain unchanged

## Format Quick Reference

### BED Format
```
GeneID  Start  End  FeatureType  [Phase]
```

### GFF3 Format
```
Chromosome  Source  Type  Start  End  Score  Strand  Phase  Attributes
```
Attributes must contain `gene_id=xxx`

## Having Issues?

- Refer to `USER_GUIDE.md` for detailed features
- Refer to `FAQ.md` for common problems
- Check if the input file format is correct

## Start Using Now!

You now understand the basic operations and can start converting your files!
