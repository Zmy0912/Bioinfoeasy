# Quick Start Guide

## 30-Second Quick Start

### 1. Launch Program

Double-click `launch_fasta_validator.bat`

Or run:
```bash
python fasta_validator.py
```

### 2. Select File

Click "Select File" button, choose a FASTA file

### 3. Validate Format

Click "Validate Format" to check for errors

### 4. Fix Errors

Click "Fix File" to automatically fix errors

### 5. Save Result

Click "Save Fixed" to save the result

## Main Features

| Feature | Description |
|---------|-------------|
| ✅ Validate | Check FASTA format |
| ✅ Highlight | Mark errors with `^` |
| ✅ Fix | Auto-fix format errors |
| ✅ Dual Mode | Standard & Strict ATCG |

## Test Files

Use the provided test files to explore features:

```bash
# Correct format
correct.fasta

# With errors
error1.fasta       # Missing >
error2.fasta       # Invalid characters
error3.fasta       # Empty header
```

## Interface Overview

The program includes four text boxes:

1. **Original File Content** - Show original file
2. **Error Position Markers** - Mark errors with `^`
3. **Fixed Content** - Show fixed result
4. **Validation Results** - Show error list

## Validation Modes

### Standard Mode
Allows all IUPAC nucleotide codes (ACGTURYSWKMBDHVN)

### Strict ATCG Mode
Only allows A, C, G, T letters

## Need Help?

Check detailed documentation:
- `README.md` - Project homepage
- `instructions.md` - Detailed guide (in English)

## System Requirements

- Python 3.x
- Windows operating system

---

**Version**: 3.0 | **Last Updated**: 2026-02-03
