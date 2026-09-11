# FASTA Validator — Test Files

Sample FASTA files for exercising the validator. Each file triggers a different
validation outcome.

| File | Purpose |
|------|---------|
| `correct.fasta` | A well-formed FASTA file (should pass) |
| `error1.fasta` | Malformed header / sequence |
| `error2.fasta` | Malformed header / sequence |
| `error3.fasta` | Malformed header / sequence |
| `non_atcg.fasta` | Contains characters outside A/T/C/G (and valid ambiguity codes) |
| `error_markers_test.fasta` | Sequences containing error marker characters |

See [`../README.md`](../README.md) for how to run the validator against these
files.

Back to [FASTA Validator](../README.md) · [project root](../../../README.md).
