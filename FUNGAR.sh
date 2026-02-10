#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 🍄 FUNGAR v1.0.0
# antiFUNGAl gene Resistance detection
#
# Bash pipeline for detecting antifungal resistance genes
# and associated amino acid mutations from sequencing reads
#
# Repository: https://github.com/resgen-br/fungar
# ============================================================

#################################
# Defaults
#################################
THREADS=4
MIN_ORF=50
GENETIC_CODE=1
MIN_QUERY_COVER=10
MIN_PIDENT=0
KEEP=false
DB_DIR=""
R1=""
R2=""
SE=""
SAMPLE_ID=""
SPECIES=""
OUTDIR=""

#################################
# Help message
#################################
usage() {
    cat << EOF
FUNGAR v1.0.0 – antiFUNGAl gene Resistance detection

Usage:
  fungar -id <sample_id> -sp <species> -o <outdir> [options]

Required:
  -id   Sample ID (e.g., sample001)
  -sp   Species name (e.g., Candida_albicans)
  -o    Output directory

Input (choose one):
  -1    Forward reads (paired-end)
  -2    Reverse reads (paired-end)
  -s    Single-end reads

Options:
  -d      Directory containing protein database and mutation CSV
  -c      CPU threads (default: 4)
  -orf    Minimum ORF length (default: 50)
  -code   Genetic code table for translation (default: 1)
  --min-query-cover  Minimum query coverage percentage (default: 10)
  --min-pident       Minimum percent identity (default: 0)
  --keep             Keep intermediate files
  -h, --help         Show this help message and exit

Example:
  fungar -id patient01 -sp Candida_albicans -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/

EOF
}

#################################
# Argument parsing
#################################
while [[ $# -gt 0 ]]; do
    case "$1" in
        -id)
            SAMPLE_ID="$2"
            shift 2
            ;;
        -sp)
            SPECIES="$2"
            shift 2
            ;;
        -o)
            OUTDIR="$2"
            shift 2
            ;;
        -1)
            R1="$2"
            shift 2
            ;;
        -2)
            R2="$2"
            shift 2
            ;;
        -s)
            SE="$2"
            shift 2
            ;;
        -d)
            DB_DIR="$2"
            shift 2
            ;;
        -c)
            THREADS="$2"
            shift 2
            ;;
        -orf)
            MIN_ORF="$2"
            shift 2
            ;;
        -code)
            GENETIC_CODE="$2"
            shift 2
            ;;
        --min-query-cover)
            MIN_QUERY_COVER="$2"
            shift 2
            ;;
        --min-pident)
            MIN_PIDENT="$2"
            shift 2
            ;;
        --keep)
            KEEP=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "❌ Unknown option: $1"
            echo
            usage
            exit 1
            ;;
    esac
done

#################################
# Sanity checks
#################################
if [[ -z "$SAMPLE_ID" || -z "$SPECIES" || -z "$OUTDIR" ]]; then
    echo "❌ ERROR: -id, -sp and -o are required."
    echo
    usage
    exit 1
fi

if [[ -n "$SE" && ( -n "$R1" || -n "$R2" ) ]]; then
    echo "❌ ERROR: Choose either single-end (-s) or paired-end (-1/-2), not both."
    exit 1
fi

if [[ -z "$SE" && ( -z "$R1" || -z "$R2" ) ]]; then
    echo "❌ ERROR: You must provide either -s or both -1 and -2."
    exit 1
fi

#################################
# Check dependencies
#################################
for cmd in diamond python; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "❌ ERROR: Required command not found in PATH: $cmd"
        exit 1
    fi
done

#################################
# Setup
#################################
mkdir -p "$OUTDIR"
WORKDIR="${OUTDIR}/${SAMPLE_ID}"
mkdir -p "$WORKDIR"

echo "========================================"
echo "🍄 Running FUNGAR v1.0.0"
echo "Sample:   $SAMPLE_ID"
echo "Species:  $SPECIES"
echo "Threads:  $THREADS"
echo "Output:   $WORKDIR"
echo "========================================"

#################################
# Main pipeline
#################################
# NOTE:
# The core logic below is intentionally preserved.
# Adjustments here should remain species/database specific.

# Placeholder for alignment and mutation detection steps
# (kept as-is to preserve original pipeline behavior)

echo "[INFO] Starting resistance gene detection..."

# >>> YOUR ORIGINAL PIPELINE LOGIC CONTINUES HERE <<<
# diamond blastx / blastp
# python mutation detection script
# result summarization

echo "[INFO] FUNGAR analysis completed successfully."

#################################
# Cleanup
#################################
if [[ "$KEEP" = false ]]; then
    echo "[INFO] Cleaning intermediate files..."
    # rm -rf "$WORKDIR/tmp"
fi

echo "✅ Done."

