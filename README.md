# FUNGAR: antiFUNGAl gene Resistance detection

FUNGAR is a Bash pipeline for detecting antifungal resistance genes from metagenomic sequencing reads.
It includes a species-specific database of genes and associated amino acid mutations.

## Features

- Supports paired-end or single-end reads
- DIAMOND-based protein alignment
- Detects known resistance mutations from CSV files
- Outputs a results file and a summary
- Configurable ORF length, genetic code, query coverage, and identity thresholds

## Instalation

FUNGAR is written in Bash and Python 3 languages and is designed for Unix-like operating systems (built and tested on Ubuntu 24.04.3 LTS).
Make sure that DIAMOND and _pandas_ library are on PATH.

- Recommended instalation
Install DIAMOND and _pandas_ with conda:

```
conda install -c bioconda -c anaconda pandas diamond 
```

- Or create a specific environment for running FUNGAR:

```
conda create -n fungar_env -c bioconda -c anaconda pandas diamond 
```
