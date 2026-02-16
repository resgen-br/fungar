#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 🍄 FUNGAR v1.0.1
# antiFUNGAl gene Resistance detection
#
# Henrique Antoniolli, Livia Kmetzsch, Charley Staats
# Universidade Federal do Rio Grande do Sul
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
INPUT_MODE="none"
SAMPLE_ID=""
SPECIES=""
OUTDIR=""

#################################
# Help
#################################
usage() {
cat <<EOF
FUNGAR v1.0.1 – antiFUNGAl gene Resistance detection

Usage:
  fungar -id <sample_id> -sp <species> -o <outdir> [input] [options]

Required:
  -id   Sample ID
  -sp   Species name (e.g. Candida_albicans)
  -o    Output directory

Input (choose one):
  -1    Forward reads (paired-end)
  -2    Reverse reads (paired-end)
  -s    Single-end reads

Options:
  -d      Database directory (protein DB + mutation CSV)
  -c      CPU threads (default: 4)
  -orf    Minimum ORF length (default: 50)
  -code   Genetic code table (default: 1)
  --min-query-cover  Minimum query coverage % (default: 10)
  --min-pident       Minimum percent identity (default: 0)
  --keep             Keep intermediate files
  -h, --help         Show this help

Example:
  fungar -id sample01 -sp Candida_albicans \\
         -1 R1.fastq.gz -2 R2.fastq.gz \\
         -o results/
EOF
}

#################################
# Argument parsing
#################################
while [[ $# -gt 0 ]]; do
    case "$1" in
        -id) SAMPLE_ID="$2"; shift 2 ;;
        -sp) SPECIES="$2"; shift 2 ;;
        -o)  OUTDIR="$2"; shift 2 ;;
        -1)  R1="$2"; INPUT_MODE="paired"; shift 2 ;;
        -2)  R2="$2"; INPUT_MODE="paired"; shift 2 ;;
        -s)  SE="$2"; INPUT_MODE="single"; shift 2 ;;
        -d)  DB_DIR="$2"; shift 2 ;;
        -c)  THREADS="$2"; shift 2 ;;
        -orf) MIN_ORF="$2"; shift 2 ;;
        -code) GENETIC_CODE="$2"; shift 2 ;;
        --min-query-cover) MIN_QUERY_COVER="$2"; shift 2 ;;
        --min-pident) MIN_PIDENT="$2"; shift 2 ;;
        --keep) KEEP=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "❌ Unknown option: $1"; usage; exit 1 ;;
    esac
done

#################################
# Sanity checks
#################################
if [[ -z "$SAMPLE_ID" || -z "$SPECIES" || -z "$OUTDIR" ]]; then
    echo "❌ ERROR: -id, -sp and -o are required."
    usage
    exit 1
fi

if [[ "$INPUT_MODE" == "paired" && ( -z "$R1" || -z "$R2" ) ]]; then
    echo "❌ ERROR: paired-end requires -1 and -2"
    exit 1
fi

if [[ "$INPUT_MODE" == "single" && -z "$SE" ]]; then
    echo "❌ ERROR: single-end requires -s"
    exit 1
fi

if [[ "$INPUT_MODE" == "none" ]]; then
    echo "❌ ERROR: input reads not specified"
    exit 1
fi

#################################
# Dependencies
#################################
for cmd in diamond python3; do
    command -v "$cmd" >/dev/null || {
        echo "❌ ERROR: missing dependency: $cmd"
        exit 1
    }
done

#################################
# Database resolution
#################################
if [[ -z "$DB_DIR" ]]; then
    if [[ -n "${FUNGAR_DB:-}" ]]; then
        DB_DIR="$FUNGAR_DB"
    elif [[ -n "${CONDA_PREFIX:-}" ]]; then
        DB_DIR="${CONDA_PREFIX}/share/FUNGAR"
    else
        DB_DIR="database"
    fi
fi

PROTEIN_DB="${DB_DIR}/${SPECIES}_db.dmnd"
MUTATION_CSV="${DB_DIR}/${SPECIES}_gene_mutations.csv"

[[ -f "$PROTEIN_DB" ]] || { echo "❌ Missing $PROTEIN_DB"; exit 1; }
[[ -f "$MUTATION_CSV" ]] || { echo "❌ Missing $MUTATION_CSV"; exit 1; }

#################################
# Setup
#################################
WORKDIR="${OUTDIR}/${SAMPLE_ID}"
mkdir -p "$WORKDIR"

LOG_FILE="${WORKDIR}/FUNGAR.${SAMPLE_ID}.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================"
echo "🍄 Running FUNGAR v1.0.1"
echo "Sample:   $SAMPLE_ID"
echo "Species:  $SPECIES"
echo "Threads:  $THREADS"
echo "Database: $DB_DIR"
echo "========================================"

#################################
# Alignment
#################################
run_diamond() {
    diamond blastx \
        -d "$PROTEIN_DB" \
        -q "$1" \
        -o "$2" \
        --threads "$THREADS" \
        --min-orf "$MIN_ORF" \
        --query-gencode "$GENETIC_CODE" \
        --query-cover "$MIN_QUERY_COVER" \
        --id "$MIN_PIDENT" \
        --strand both \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq_translated sseq \
        --masking 0
}

ALIGNMENT_FILES=""

if [[ "$INPUT_MODE" == "paired" ]]; then
    R1_ALN="${WORKDIR}/${SAMPLE_ID}_R1.aln"
    R2_ALN="${WORKDIR}/${SAMPLE_ID}_R2.aln"
    run_diamond "$R1" "$R1_ALN"
    run_diamond "$R2" "$R2_ALN"
    ALIGNMENT_FILES="${R1_ALN};${R2_ALN}"
else
    SE_ALN="${WORKDIR}/${SAMPLE_ID}_SE.aln"
    run_diamond "$SE" "$SE_ALN"
    ALIGNMENT_FILES="${SE_ALN}"
fi

#################################
# Mutation detection (Python)
#################################
python3 <<EOF
import pandas as pd
from collections import defaultdict
import os

mut = pd.read_csv("$MUTATION_CSV")
required = ['gene','position','reference','mutation','fungicide']
if not all(c in mut.columns for c in required):
    raise ValueError("Mutation CSV missing required columns")

alignment_files = "$ALIGNMENT_FILES".split(";")

results = []

for aln in alignment_files:
    if not os.path.exists(aln):
        continue

    with open(aln) as f:
        for line in f:
            fields = line.rstrip().split("\\t")
            if len(fields) < 13:
                continue

            qseqid, gene = fields[0], fields[1]
            sstart, send = int(fields[8]), int(fields[9])
            sseq = fields[12]

            gene_muts = mut[mut["gene"] == gene]
            for _, row in gene_muts.iterrows():
                pos = int(row["position"])

                if sstart <= send:
                    if not (sstart <= pos <= send):
                        continue
                    aa_pos = pos - sstart
                else:
                    if not (send <= pos <= sstart):
                        continue
                    aa_pos = sstart - pos

                if 0 <= aa_pos < len(sseq) and sseq[aa_pos] == row["mutation"]:
                    results.append({
                        "Sample": qseqid,
                        "Gene": gene,
                        "Position": pos,
                        "Reference": row["reference"],
                        "Mutation": row["mutation"],
                        "Fungicide": row["fungicide"],
                        "Alignment": os.path.basename(aln)
                    })

final = pd.DataFrame(results).drop_duplicates()
final.to_csv("${WORKDIR}/${SAMPLE_ID}_results.csv", index=False)

summary = defaultdict(int)
for r in results:
    key = (r["Gene"], r["Position"], r["Reference"], r["Mutation"], r["Fungicide"])
    summary[key] += 1

summary_df = pd.DataFrame([
    {"Gene":k[0],"Position":k[1],"Reference":k[2],
     "Mutation":k[3],"Fungicide":k[4],"Support_Reads":v}
    for k,v in summary.items()
])

summary_df.to_csv("${WORKDIR}/${SAMPLE_ID}_summary.csv", index=False)
EOF

#################################
# Cleanup
#################################
if [[ "$KEEP" == false ]]; then
    rm -f ${ALIGNMENT_FILES//;/ }
fi

echo "✅ FUNGAR analysis completed successfully."
