[![arXiv](https://img.shields.io/badge/arXiv-2602.16728-b31b1b.svg)](https://arxiv.org/abs/2602.16728)

# FUNGAR: antiFUNGAl gene Resistance detection

FUNGAR is a Bash pipeline for detecting antifungal resistance mutations directly
from metagenomic short reads. It uses DIAMOND for protein-level alignment and
species-specific mutation databases to call known resistance variants.

## Features

- Supports paired-end or single-end reads
- DIAMOND-based protein (blastx) alignment
- Detects known resistance mutations from CSV databases
- **New:** Minimum read-support threshold (`--min-support`) to filter low-confidence calls
- **New:** In paired-end mode, require both reads (mates) to align to the reference (`--require-both-mates`, default: true) to drastically reduce false positives
- **New:** `--list_species` option to quickly list available species in the database
- **New:** Outputs a self-contained **HTML report** with tables, charts, and interpretation notes
- **New:** Runtime and peak-memory tracking per alignment step
- **New:** Extended output columns: `Subject_Accession`, `Source` (clinical/agricultural), `Strand`
- Configurable ORF length, genetic code, query coverage, and identity thresholds

## Installation

FUNGAR requires **DIAMOND** and **pandas** (Python 3). Tested on Ubuntu 24.04 LTS.

```bash
# Recommended — create a dedicated conda environment
conda create -n fungar_env -c bioconda -c anaconda pandas diamond
conda activate fungar_env
```

## Running FUNGAR

```
Usage:
  fungar -id <sample_id> -sp <species> -o <outdir> [input] [options]

Required:
  -id   Sample ID (e.g. patient123)
  -sp   Species name matching a database entry (e.g. Aspergillus_fumigatus)
  -o    Output directory

Input (choose one):
  -1    Forward reads (R1.fastq.gz) — paired-end mode
  -2    Reverse reads (R2.fastq.gz) — paired-end mode
  -s    Single-end reads (SE.fastq.gz)

Options:
  -d      Directory containing protein DB (.dmnd) + mutation CSV
  -c      CPU threads (default: 4)
  -orf    Minimum ORF length in nt (default: 10)
  -code   Genetic code table (default: 1; standard)
  --min-query-cover  Minimum query coverage % (default: 10)
  --min-pident       Minimum % identity (default: 95)
  --min-support      Minimum reads supporting a mutation to report it (default: 1)
  --require-both-mates  Require both mates to align to reference (true/false, default: true)
  --keep             Keep intermediate alignment files
  --list_species     List available species in the database and exit
  -h, --help         Show this help
```

### Example

```bash
fungar \
    -id sample01 -sp Aspergillus_fumigatus \
    -1 R1.fastq.gz -2 R2.fastq.gz \
    -o results/ -d database/ \
    --min-support 2
```

## Output Files

| File | Description |
|------|-------------|
| `<id>_results.csv` | Per-read mutation hits (before deduplication). New columns: `Subject_Accession`, `Source`, `Strand` |
| `<id>_summary.csv` | Deduplicated mutations with `Support_Reads` count, filtered by `--min-support` |
| `<id>_report.html` | Self-contained HTML report with charts and interpretation notes |
| `<id>_runtime.txt` | Wall-clock time and peak RSS memory for each DIAMOND run |
| `FUNGAR.<id>.log`  | Full pipeline log |

## HTML Report

The HTML report (`*_report.html`) is generated automatically at the end of each run.
It includes:
- Summary statistics cards (mutations detected, genes affected, total read support)
- Resistance mutation table colour-coded by drug class (azoles, echinocandins, polyenes, QoI)
- Bar chart of read support per mutation
- Collapsed per-read alignment table
- Run parameters footer

The file is fully self-contained (no external dependencies at view-time) and can be
shared, printed, or saved as PDF via the browser.

## Benchmarking

`fungar_benchmark.py` provides quantitative performance evaluation:

```bash
python3 fungar_benchmark.py \
    -sp Aspergillus_fumigatus \
    -d  database/ \
    -o  benchmark_out/ \
    --fungar fungar \
    --depths 0.2 0.5 1 5 10 20 \
    --read-length 150 \
    --variable-read-lengths \
    --min-read-length 20 \
    --threads 4
```

The script:
1. Reads the species protein FASTA and mutation CSV
2. Generates synthetic FASTQ reads stochastically (using a Poisson read-count sampler based on target depth) by reverse-translating protein sequences with each known mutation planted at the target position (positive runs) and without (wild-type runs for FPR estimation)
3. Simulates variable read lengths (mimicking fastp-trimmed inputs) via `--variable-read-lengths`, which randomly trims 30% of reads to sizes between `--min-read-length` (default 20bp) and full length
4. Runs FUNGAR on every synthetic dataset
5. Calculates sensitivity, precision, F1, and false-positive rate per mutation and depth
6. Outputs `benchmark_<species>.csv` and `benchmark_<species>_report.html`

> **Note**: Benchmarks use synthetic stochastically generated reads. Real-world sensitivity will
> depend on sequencing error, read length variation, species cross-mapping, and sample complexity.

## Known Limitations and Boundary Cases

### Reads partially covering the diagnostic residue
Reads where the diagnostic amino acid position falls **outside** the aligned region
are silently filtered. This is intentional and conservative. Increasing `--min-query-cover`
reduces partial-alignment noise.

### Multiple ORFs and frameshifts
DIAMOND blastx evaluates all six reading frames and reports the best-scoring ORF per read.
FUNGAR uses the translated query (`qseq_translated`) from that ORF. Reads with internal
stop codons or frameshifts are handled by DIAMOND's internal ORF caller (`--min-orf`).

### Heterozygosity and Variant Allele Frequency
Because FUNGAR relies on translated searches (aligning translated nucleotide reads against a protein database), it abstracts away nucleotide-level variation. Consequently, synonymous heterozygous alleles are completely invisible to the pipeline. Furthermore, FUNGAR only reports reads supporting a functional mutation and does not quantify wild-type read depth at the locus. Therefore, you **cannot** use FUNGAR's read support to calculate Variant Allele Frequency (VAF) or infer true heterozygosity. Quantitative allele-frequency analysis requires mapping reads to a nucleotide reference genome and using dedicated variant-calling tools.

### Closely related species and false positives
Proteins conserved across species (e.g. β-tubulin, Cyp51) can produce cross-species
alignments. The `Subject_Accession` column in `*_results.csv` identifies the specific
database entry that matched; verify accession identity before clinical interpretation.

### Confirmed vs. candidate mutations
Mutations detected computationally are **candidates**. Biological confirmation
(culture-based MIC, whole-gene Sanger sequencing) is required before clinical use.
The HTML report flags all detections accordingly.

### Paired-end mate alignment filtering
By default in paired-end mode, FUNGAR requires **both reads (mates)** of a pair to successfully align to the reference protein database. If only one mate aligns, the read is filtered out to reduce false positives (e.g. from spurious cross-species alignments of a single highly-conserved mate). Users can deactivate this safety feature by setting `--require-both-mates false`.

## Comparison with Related Tools

| Feature | FUNGAR | ChroQueTas | ResFinder (fungi) |
|---------|--------|------------|-------------------|
| Short-read (metagenomics) | ✅ | ❌ (assembly) | ❌ |
| Protein-level alignment | ✅ | ✅ | ❌ |
| Point mutation detection | ✅ | ✅ | Partial |
| Multi-species DB | ✅ | Limited | Limited |
| HTML report | ✅ | ❌ | ✅ |
| Read-support threshold | ✅ | N/A | N/A |

FUNGAR is designed as a **complement** to assembly-based tools: it sacrifices
per-gene assembly completeness for the ability to detect resistance mutations
directly from short metagenomic reads without a reference genome.
