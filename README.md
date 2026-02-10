# FUNGAR: antiFUNGAl gene Resistance detection

FUNGAR is a Bash pipeline for detecting antifungal resistance genes from metagenomic sequencing reads.
It includes a species-specific database of genes and associated amino acid mutations.

## Features

- Supports paired-end or single-end reads
- DIAMOND-based protein alignment
- Detects known resistance mutations from CSV files
- Outputs a results file and a summary
- Configurable ORF length, genetic code, query coverage, and identity thresholds

## Installation

FUNGAR is written in Bash and Python 3 languages and is designed for Unix-like operating systems (built and tested on Ubuntu 24.04.3 LTS).
Make sure that DIAMOND and _pandas_ library are on PATH.

- Recommended instalation

Install DIAMOND and _pandas_ with conda:

```
conda install -c bioconda -c anaconda pandas diamond 
```

Or create a specific environment for running FUNGAR:

```
conda create -n fungar_env -c bioconda -c anaconda pandas diamond 
```

## Running FUNGAR
```
Usage:
$ fungar -id <sample_id> -sp <species> -o <outdir> [-d <database_dir>] [-orf <min_orf_length>] [-code <int>] [input options]

Required:
-id   Sample ID (e.g., patient123)
-sp   Species name (e.g., Candida_albicans)
-o    Output directory

Input options (choose one):
-1    Forward reads (R1.fastq.gz) - for paired-end mode
-2    Reverse reads (R2.fastq.gz) - for paired-end mode
or
-s    Single-end reads (SE.fastq.gz)

Options:
-d      Directory containing protein DB + mutation CSV
-c      CPU threads (default: 4)
-orf    Minimum ORF length (default: 50)
-code   Genetic code table for translation (default: 1)
--min-query-cover <int>   Minimum query coverage % (default: 10)
--min-pident <int>        Minimum % identity (default: 0)
--keep  Keep intermediate alignment files
```
