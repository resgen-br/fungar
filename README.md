# FUNGAR: antiFUNGAl gene Resistance detection

FUNGAR is a Bash pipeline for detecting antifungal resistance genes from metagenomic sequencing reads.
It includes a species-specific database of genes and associated amino acid mutations.

## Features

- Supports paired-end or single-end reads
- DIAMOND-based protein alignment
- Detects known resistance mutations from CSV files
- Outputs a results file and a summary
- Configurable ORF length, genetic code, query coverage, and identity thresholds
