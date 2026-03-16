# FASTA Format Validator and Fixer

A powerful tool for FASTA format validation, error highlighting, and automatic fixing with a graphical user interface.

![Python](https://img.shields.io/badge/Python-3.x-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Version](https://img.shields.io/badge/Version-3.0-orange.svg)

## Features

- ✅ **FASTA Format Validation** - Check if files follow FASTA format specifications
- ✅ **Error Position Highlighting** - Mark invalid characters with `^` symbols
- ✅ **Automatic Fixing** - One-click fix for common format errors
- ✅ **Dual Validation Modes** - Standard mode and strict ATCG mode
- ✅ **Graphical Interface** - Intuitive four-window comparison display
- ✅ **Detailed Error Reports** - Include line numbers, error types, and specific descriptions

## Interface Overview

The program features four text boxes (top to bottom):

1. **Original File Content** - Display original FASTA file content
2. **Error Position Markers** - Display error markers (`^` indicates invalid characters)
3. **Fixed Content** - Display fixed FASTA file content
4. **Validation Results/Error Information** - Display detailed error list

## Installation

### System Requirements

- Python 3.x
- Windows operating system
- Supported file encodings: UTF-8, GBK, GB2312, Latin-1

### Quick Start

Double-click `launch_fasta_validator.bat` to launch the program

Or use command line:

```bash
python fasta_validator.py
```

## Usage

### Basic Operations

1. Click "Select File" to load a FASTA file
2. Choose validation mode (optional):
   - **Standard Mode**: Allows all IUPAC nucleotide codes (ACGTURYSWKMBDHVN)
   - **Strict ATCG Mode**: Only allows A, C, G, T letters
   - Toggle "Strict ATCG Mode" checkbox to switch modes
3. Enable error highlighting (optional):
   - Check "Highlight Error Positions" checkbox
   - Error positions will be displayed in the "Error Position Markers" text box
4. Click "Validate Format" to check file format
   - Error information will appear in the "Validation Results/Error Information" text box
   - Error position markers will appear in the "Error Position Markers" text box
5. If errors found, click "Fix File"
   - Fixed content will appear in the "Fixed Content" text box
   - Fixing is based on current selected mode
6. Review fixed result and click "Save Fixed" if satisfied

### Error Position Highlighting

When error highlighting is enabled, the program marks invalid characters with `^` in the sequence:

```
>seq1_with_errors
Line 2: ATCGATCG12345
            ^^^^^
            ^ Invalid characters: 12345

Line 3: ATCG@#$%^&*()
            ^^^^^^^^^
            ^ Invalid characters: #$%&()*@^
```

### Color Coding

- 🔵 **Blue**: Sequence header lines
- ⚫ **Black**: Original sequence lines
- 🔴 **Red**: Error marker lines (`^` symbol)
- ⚪ **Gray**: Informational text

## Validation Modes

### Standard Mode

- **Valid characters**: A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N, -
- **Use case**: General FASTA validation, supports all IUPAC nucleotide codes
- **Features**: Allows degenerate bases and uncertainty symbols

### Strict ATCG Mode

- **Valid characters**: A, C, G, T, -
- **Use case**: Standard DNA sequence validation only
- **Features**: Only allows four standard base letters and gap symbols

**Example Comparison:**

| Sequence Content | Standard Mode | Strict ATCG Mode |
|-----------------|---------------|-----------------|
| `ATCGNNNNURYS` | ✅ Valid | ❌ Invalid |
| `ATCGATCGATCG` | ✅ Valid | ✅ Valid |
| `ATCG12345` | ❌ Invalid | ❌ Invalid |

## Supported File Formats

.fasta, .fa, .fas, .fna, .ffn, .faa, .frn, .txt

## Common Errors and Fixes

### 1. Line 1 must start with '>' (header line)

- **Error example**: First line doesn't have `>` symbol
- **Fix**: Automatically add header

### 2. Sequence data without header line

- **Error example**: Sequence data found without a header
- **Fix**: Automatically add header line

### 3. Header line missing sequence ID

- **Error example**: `>` with nothing after it
- **Fix**: Automatically generate default ID (e.g., sequence_1, sequence_2...)

### 4. Contains invalid characters

- **Standard mode error example**: Sequence contains numbers, special symbols
- **Strict mode error example**: Sequence contains N, U, R, etc. (non-ATCG)
- **Fix**: Automatically remove non-sequence characters based on current mode
- **Display format**:
  - Standard mode: `Line X: Contains invalid characters 'XXX' (expected ACGTURYSWKMBDHVN)`
  - Strict mode: `Line X: Contains invalid characters 'XXX' (expected ATCG)`

### 5. No sequence headers found

- **Error example**: Entire file has no lines starting with `>`
- **Fix**: Add headers for each sequence segment

### 6. File is empty

- **Error message**: File is empty
- **Cannot be fixed automatically**

## Test Files

The project includes test files for functionality demonstration:

- `correct.fasta` - Correct FASTA format file
- `error1.fasta` - Missing `>` symbol error
- `error2.fasta` - Contains invalid characters
- `error3.fasta` - Empty header error
- `error_markers_test.fasta` - Test error position highlighting
- `non_atcg.fasta` - Non-ATCG character test

## Project Structure

```
.
├── fasta_validator.py          # Main program
├── 启动FASTA验证器.bat          # Quick launch script
├── README.md                   # Project overview (this file)
├── CHANGELOG.md                # Version history
├── LICENSE.md                  # MIT License
├── 使用说明.md                 # Detailed usage guide
├── QUICK_START.md              # Quick start guide
├── correct.fasta               # Test file: correct format
├── error1.fasta               # Test file: missing >
├── error2.fasta               # Test file: invalid characters
├── error3.fasta               # Test file: empty header
├── error_markers_test.fasta   # Test file: error markers
└── non_atcg.fasta             # Test file: non-ATCG chars
```

## Important Notes

1. Please backup original files before fixing
2. Please manually review fixed results before saving
3. For complex format errors, manual editing is recommended
4. The program preserves all valid sequence information

## Technical Details

### Error Marking Algorithm

1. Iterate through each line character by character
2. Check if character is in valid character set
3. Valid characters: mark with space
4. Invalid characters: mark with `^`
5. Collect list of invalid characters
6. Display marker line and invalid character information

### Performance Optimizations

- Only calculate error positions when checkbox is enabled
- Cache error position information, no recalculation
- Use async scrolling for text boxes to prevent UI freezing

## Version History

### v3.0 (2026-02-03)

- ✨ New error position highlighting feature
- ✨ New strict ATCG validation mode
- 🎨 Optimized interface layout, added four text boxes
- 🐛 Fixed error information window display issues
- 📝 Updated documentation

### v2.0

- ✨ New strict ATCG validation mode
- 🎨 Improved interface layout
- 🐛 Fixed validation logic issues

### v1.0

- 🎉 Initial release
- ✨ Basic FASTA validation features
- ✨ Error detection and fixing
- ✨ Graphical user interface

## Contributing

Issues and improvement suggestions are welcome!

**When reporting issues, please include**:
- Version number
- Operating system
- Problem description
- Steps to reproduce
- Relevant files or screenshots

## License

MIT License

## Authors

Zmy0912

## Acknowledgments

Thanks to all researchers and developers who have contributed to FASTA format specifications and bioinformatics tools.

---

**Quick Start**: Double-click `launch_fasta_validator.bat` to get started!

**Detailed Guide**: See `instructions.md` for more details.
