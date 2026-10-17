# Analyzer

Analysis tools for chloroplast genomes and repetitive sequences.

## Contents

| Module | Description |
|--------|-------------|
| [`Chloroplast_genome/`](Chloroplast_genome/) | Chloroplast genome analysis: codon usage bias, assembly comparison, four-region analysis |
| [`ENC_GC3s_Analyzer/`](ENC_GC3s_Analyzer/) | Quantifies genes below / on / above Wright's ENC-GC3s neutral expectation curve |
| [`GC_ENC_Analyzer/`](GC_ENC_Analyzer/) | Batch GC1/GC2/GC3/GC3s/GC12 and ENC (CodonW 1.4.4 compatible) with GC12–GC3 regression |
| [`KaKs_Analyzer/`](KaKs_Analyzer/) | Ka/Ks filtering, grouping and missing-value imputation (`.xlsx` / `.kaks`) |
| [`RepeatAnalyzer/`](RepeatAnalyzer/) | Direct and palindromic repeat statistics from Vmatch `repfind` output |
| [`SSR_Region_Classifier/`](SSR_Region_Classifier/) | Classify every SSR into CDS / Intron / IGS using MISA + GenBank files |
| [`misa_tool/`](misa_tool/) | MISA SSR toolkit — result parser (`misa_parser`) and figure viewer (`misa_ssr_viewer`) |
| [`sequence_repeat_analyzer/`](sequence_repeat_analyzer/) | Comprehensive repeat sequence analysis |

## Quick Start

```bash
# MISA SSR parsing
cd misa_tool/misa_parser
pip install -r requirements.txt
python misa_parser.py

# Sequence repeat analysis
cd sequence_repeat_analyzer
pip install -r requirements.txt
python sequence_repeat_analyzer.py

# Ka/Ks filtering and imputation
cd KaKs_Analyzer
python kaks_analyzer.py

# ENC-GC3s neutral expectation test
cd ENC_GC3s_Analyzer
python ENC_GC3s_Analyzer.py

# Batch GC / ENC codon usage bias statistics
cd GC_ENC_Analyzer
python GC_ENC_Analyzer.py

# SSR region classification
cd SSR_Region_Classifier
python ssr_region_classifier.py

# Repeat counting (standard library only)
cd RepeatAnalyzer
python repeat_analyzer.py
```

Each module has its own `README.md` with detailed usage instructions.

Back to the [project root](../README.md).
