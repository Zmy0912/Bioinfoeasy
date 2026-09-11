# GenBank to mVISTA Converter

Converts GenBank/GFF3 annotation files into the input format expected by
mVISTA.

## Contents

| Version | Files |
|---------|-------|
| [`v1.0.0/`](v1.0.0/) | `genbank_to_mvista.py`, `README.md`, `requirements.txt` |
| [`v1.1.0/`](v1.1.0/) | `genbank_to_mvista.py`, `README.md`, `requirements.txt` |

`v1.1.0` is the current release. Note that the dependency set was slimmed down
between versions.

## Requirements

- Python 3.7+
- `biopython`, `PyQt5` (v1.1.0); `pandas`, `openpyxl`, `bcbio-gff` (v1.0.0)

```bash
pip install -r v1.1.0/requirements.txt
```

## Usage

```bash
cd v1.1.0
python genbank_to_mvista.py
```

Back to [Seqtool](../README.md) · [project root](../../README.md).
