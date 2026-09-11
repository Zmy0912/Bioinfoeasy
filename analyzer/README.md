# Analyzer

Analysis tools for chloroplast genomes and repetitive sequences.

## Contents

| Module | Description |
|--------|-------------|
| [`Chloroplast_genome/`](Chloroplast_genome/) | Chloroplast genome analysis: codon usage bias, assembly comparison, four-region analysis |
| [`misa_tool/`](misa_tool/) | MISA SSR (Simple Sequence Repeat) motif parser |
| [`sequence_repeat_analyzer/`](sequence_repeat_analyzer/) | Comprehensive repeat sequence analysis |

## Quick Start

```bash
# MISA SSR parsing
cd misa_tool
pip install -r requirements.txt
python misa_parser.py

# Sequence repeat analysis
cd sequence_repeat_analyzer
pip install -r requirements.txt
python sequence_repeat_analyzer.py
```

Each module has its own `README.md` with detailed usage instructions.

Back to the [project root](../README.md).
