# Seqinfo — Genome Annotation Analyzer

Extracts comprehensive information from genome annotation files.

## Contents

| Version | Description |
|---------|-------------|
| [`v1.0.0/`](v1.0.0/) | Initial release |
| [`v3.0.0/`](v3.0.0/) | Current release |

Each version contains `genome_annotation_analyzer.py`, `README.md` and
`requirements.txt`.

## What it extracts

- Basic information (ID, species, length)
- Gene statistics (CDS, rRNA, tRNA, ncRNA)
- GC content and topology
- Taxonomy and metadata

## Supported formats

- GenBank (`.gb`, `.gbk`)
- GFF / GFF3
- FASTA

## Requirements

```bash
cd v3.0.0
pip install -r requirements.txt
```

- `biopython>=1.81`
- `pandas>=1.5.0`
- `openpyxl>=3.0.0`
- `bcbio-gff>=0.6.0` (optional — the tool falls back to a simplified parser)

## Usage

```bash
cd v3.0.0
python genome_annotation_analyzer.py chloroplast.gb
```

Back to the [project root](../README.md).
