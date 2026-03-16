# Batch Renamer Tool

This is a graphical user interface tool for batch replacing names in files.

## Features

- Support loading TXT format name mapping files (first column is original name, second column is new name, separated by tab)
- Support selecting any type of target file for name replacement
- Provide split-screen preview, view original content and replaced content simultaneously
- Display replacement statistics
- Simple and intuitive graphical user interface

## Usage

### 1. Run the Program

```bash
python batch_renamer.py
```

### 2. Operation Steps

1. **Select Name Mapping File**: Click the "Browse..." button to select a TXT file containing the correspondence between original names and new names
   - File format requirement: two columns per line, separated by tab (TAB)
   - First column: original name
   - Second column: new name

2. **Select Target File**: Click the "Browse..." button to select the file that needs name replacement

3. **Select Output File**: Click the "Browse..." button to select or enter the output file path

4. **Load Name Mapping**: Click the "Load Name Mapping" button, the program will automatically parse the mapping file

5. **Preview Replacement**: Click the "Preview Replacement" button, the left side shows the original content, the right side shows the content after replacement

6. **Execute Replacement**: After confirming the preview is correct, click the "Execute Replacement" button, the program will save the result to the output file

## Example

### Name Mapping File Example (new text document.txt)

```
AccessionNumber    ScientificName
OQ851342.1    Asparagus officinalis
NC_086747.1    Asparagus subscandens
...
```

### Notes

- The name mapping file must use tab (TAB) to separate the two columns of data
- The program replaces names from longest to shortest based on old name length to avoid partial matching issues
- The replacement operation is irreversible, it is recommended to preview and confirm before executing
- Supports all text format files

## System Requirements

- Python 3.6+
- tkinter (usually installed with Python)

## Dependencies

This program only uses Python standard library, no additional dependencies need to be installed.
