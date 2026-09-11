# Quick Start Guide

## Launch Program

Double-click `start_converter.bat` file to launch the program.

## Convert Files

1. Click "Browse" button to select `基因成员.gff3` file
2. The program will automatically set output file to `基因成员.bed`
3. Keep default options (convert all feature types)
4. Click "Start Conversion" button
5. After conversion completes, view the generated BED file

## View Results

After conversion, you can find the generated `.bed` file in the program directory. Open it with a text editor to view.

## Customize Options

- To convert only specific types (e.g., only CDS), uncheck other options
- You can manually modify the name and save location of BED output file

## Common Questions

**Q: Program won't start?**
A: Make sure Python 3.6 or higher is installed

**Q: Converted file is empty?**
A: Check if GFF3 file format is correct and contains selected feature types

**Q: Chinese characters display incorrectly?**
A: Ensure GFF3 file is opened with UTF-8 encoding

For more help, please refer to README.md file.
