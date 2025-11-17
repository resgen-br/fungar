#!/bin/bash

##################################################
# 🍄 FUNGAR v1.0.0: antiFUNGAl gene Resistance detection
# by Henrique Antoniolli, Livia Kmetzsch and Charley Staats
# Universidade Federal do Rio Grande do Sul
# Contact: staats@ufrgs.br
##################################################

# Defaults
CPUS=4
KEEP_FILES=false
INPUT_MODE="none"
MIN_ORF=50
MIN_QUERY_COVER=10
MIN_PIDENT=0
GENETIC_CODE=1
DATABASE_DIR=""

usage() {
    echo "Usage: $0 -id <sample_id> -sp <species> -o <outdir> [-d <database_dir>] [-orf <min_orf_length>] [-code <int>] [input options]"
    echo "Required:"
    echo "  -id   Sample ID (e.g., patient123)"
    echo "  -sp   Species name (e.g., Candida_albicans)"
    echo "  -o    Output directory"
    echo "Input options (choose one):"
    echo "  -1    Forward reads (R1.fastq.gz) - for paired-end mode"
    echo "  -2    Reverse reads (R2.fastq.gz) - for paired-end mode"
    echo "  -s    Single-end reads (SE.fastq.gz)"
    echo "Options:"
    echo "  -d      Directory containing protein DB + mutation CSV"
    echo "  -c      CPU threads (default: 4)"
    echo "  -orf    Minimum ORF length (default: 50)"
    echo "  -code   Genetic code table for translation (default: 1)"
    echo "  --min-query-cover <int>   Minimum query coverage % (default: 10)"
    echo "  --min-pident <int>        Minimum % identity (default: 0)"
    echo "  --keep  Keep intermediate alignment files"
    exit 1
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -id) SAMPLE="$2"; shift ;;
        -sp) SPECIES="$2"; shift ;;
        -1) R1="$2"; INPUT_MODE="paired"; shift ;;
        -2) R2="$2"; INPUT_MODE="paired"; shift ;;
        -s) SE="$2"; INPUT_MODE="single"; shift ;;
        -o) OUTDIR="$2"; shift ;;
        -d) DATABASE_DIR="$2"; shift ;;
        -c) CPUS="$2"; shift ;;
        -orf) MIN_ORF="$2"; shift ;;
        -code) GENETIC_CODE="$2"; shift ;;
        --min-query-cover) MIN_QUERY_COVER="$2"; shift ;;
        --min-pident) MIN_PIDENT="$2"; shift ;;
        --keep) KEEP_FILES=true ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

echo
echo "🍄 FUNGAR v1.0.0: pipeline for antiFUNGAl gene Resistance detection"
echo

# Validate required arguments
if [[ -z "$SAMPLE" || -z "$SPECIES" || -z "$OUTDIR" ]]; then
    echo "ERROR: Missing required arguments!"
    usage
fi

# Validate input mode
if [[ "$INPUT_MODE" == "none" ]]; then
    echo "ERROR: No input files specified. Must use either -1/-2 (paired) or -s (single)"
    usage
elif [[ "$INPUT_MODE" == "paired" && (-z "$R1" || -z "$R2") ]]; then
    echo "ERROR: Paired-end mode requires both -1 and -2 arguments"
    usage
elif [[ "$INPUT_MODE" == "single" && -z "$SE" ]]; then
    echo "ERROR: Single-end mode requires -s argument"
    usage
fi

#############################################################
# DATABASE HANDLING
#############################################################
# Priority:
#   1. -d
#   2. FUNGAR_DB environment variable
#   3. conda installation
#   4. fallback: "database"
if [ -z "$DATABASE_DIR" ]; then
    if [ -n "$FUNGAR_DB" ]; then
        DATABASE_DIR="$FUNGAR_DB"
    elif [ -n "$CONDA_PREFIX" ]; then
        DATABASE_DIR="${CONDA_PREFIX}/share/FUNGAR"
    else
        DATABASE_DIR="database"
    fi
fi

# Create output directory
mkdir -p "$OUTDIR" || { echo "ERROR: Could not create $OUTDIR"; exit 1; }

# Set log file
LOG_FILE="${OUTDIR}/FUNGAR.${SAMPLE}.log"

# Start logging
{
echo "################################################"
echo "# FUNGAR v1.0.0: antiFUNGAl gene Resistance detection #"
echo "################################################"
echo "Started: $(date)"
echo "Sample: $SAMPLE"
echo "Species: $SPECIES"
echo "Database directory: $DATABASE_DIR"
echo "Input mode: $INPUT_MODE"
if [ "$INPUT_MODE" == "paired" ]; then
    echo "- R1: $R1"
    echo "- R2: $R2"
else
    echo "- Single-end: $SE"
fi
echo "Parameters:"
echo "- CPUs: $CPUS"
echo "- Minimum ORF length: $MIN_ORF"
echo "- Genetic code: $GENETIC_CODE"
echo "- Minimum query cover: $MIN_QUERY_COVER"
echo "- Minimum identity: $MIN_PIDENT"
echo "- Keep intermediates: $KEEP_FILES"
echo ""
} > "$LOG_FILE"

# File validation
check_file() {
    if [ ! -f "$1" ]; then
        echo "ERROR: File not found - $1" | tee -a "$LOG_FILE"
        exit 1
    fi
}

echo "[INFO] Verifying input files..." | tee -a "$LOG_FILE"
if [ "$INPUT_MODE" == "paired" ]; then
    check_file "$R1"
    check_file "$R2"
else
    check_file "$SE"
fi

# Verify species-specific files
PROTEIN_DB="${DATABASE_DIR}/${SPECIES}_db.dmnd"
MUTATION_CSV="${DATABASE_DIR}/${SPECIES}_gene_mutations.csv"
check_file "$PROTEIN_DB"
check_file "$MUTATION_CSV"

# Define output files
if [ "$INPUT_MODE" == "paired" ]; then
    R1_OUT="${OUTDIR}/${SAMPLE}_R1.aln"
    R2_OUT="${OUTDIR}/${SAMPLE}_R2.aln"
else
    SE_OUT="${OUTDIR}/${SAMPLE}_SE.aln"
fi
FINAL_RESULTS="${OUTDIR}/${SAMPLE}_results.csv"
SUMMARY_RESULTS="${OUTDIR}/${SAMPLE}_summary.csv"

# DIAMOND alignment function
run_diamond() {
    local input=$1
    local output=$2
    local read=$3

    echo "[INFO] Running DIAMOND on $read reads..." | tee -a "$LOG_FILE"

    diamond blastx \
        -d "$PROTEIN_DB" \
        -q "$input" \
        -o "$output" \
        --threads "$CPUS" \
        --min-orf "$MIN_ORF" \
        --query-gencode "$GENETIC_CODE" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq_translated sseq \
        --strand both \
        --query-cover "$MIN_QUERY_COVER" \
        --id "$MIN_PIDENT" \
        --masking 0 2>> "$LOG_FILE"

    if [ $? -ne 0 ]; then
        echo "ERROR: DIAMOND failed on $read reads" | tee -a "$LOG_FILE"
        exit 1
    fi
}

# Run alignments
if [ "$INPUT_MODE" == "paired" ]; then
    run_diamond "$R1" "$R1_OUT" "R1"
    run_diamond "$R2" "$R2_OUT" "R2"
else
    run_diamond "$SE" "$SE_OUT" "single-end"
fi

# Mutation detection with Python
echo "[INFO] Detecting resistance mutations..." | tee -a "$LOG_FILE"
python3 <<EOF | tee -a "$LOG_FILE"
import pandas as pd
import os
from collections import defaultdict

mutations = pd.read_csv("$MUTATION_CSV")
required_cols = ['gene','position','reference','mutation','fungicide']
if not all(col in mutations.columns for col in required_cols):
    raise ValueError(f"CSV missing required columns: {required_cols}")

all_results = []

# Select alignments
if "$INPUT_MODE" == "paired":
    alignments = [("${R1_OUT}", "R1"), ("${R2_OUT}", "R2")]
else:
    alignments = [("${SE_OUT}", "SE")]

for aln_file, read_type in alignments:
    if not os.path.exists(aln_file):
        continue
    with open(aln_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 13:
                continue
            qseqid, sseqid = fields[0], fields[1]
            sstart, send = int(fields[8]), int(fields[9])
            sseq = fields[12]
            
            gene_muts = mutations[mutations['gene'] == sseqid]
            for _, row in gene_muts.iterrows():
                pos = int(row['position'])
                # Handle forward/reverse alignments
                if sstart <= send:
                    if sstart <= pos <= send:
                        aa_pos = pos - sstart
                    else:
                        continue
                else:  # reverse
                    if send <= pos <= sstart:
                        aa_pos = sstart - pos
                    else:
                        continue
                if 0 <= aa_pos < len(sseq) and sseq[aa_pos] == row['mutation']:
                    all_results.append({
                        'Sample': qseqid,
                        'Gene': sseqid,
                        'Position': pos,
                        'Reference': row['reference'],
                        'Mutation': row['mutation'],
                        'Fungicide': row['fungicide'],
                        'Read': read_type
                    })

# Save final results
if all_results:
    df = pd.DataFrame(all_results).drop_duplicates(
        subset=['Sample','Gene','Position','Mutation'],
        keep='first'
    )
    df.to_csv("$FINAL_RESULTS", index=False)
else:
    pd.DataFrame(columns=['Sample','Gene','Position','Reference','Mutation','Fungicide','Read']).to_csv("$FINAL_RESULTS", index=False)

# Summary counting all supporting reads
summary_counter = defaultdict(int)
for res in all_results:
    key = (res['Gene'], res['Position'], res['Reference'], res['Mutation'], res['Fungicide'])
    summary_counter[key] += 1

summary_data = []
for (gene,pos,ref,mut,fung), count in summary_counter.items():
    summary_data.append({
        'Gene': gene,
        'Position': pos,
        'Reference': ref,
        'Mutation': mut,
        'Fungicide': fung,
        'Support_Reads': count
    })
pd.DataFrame(summary_data).to_csv("$SUMMARY_RESULTS", index=False)
EOF

# Cleanup intermediate files
if ! $KEEP_FILES; then
    echo "[INFO] Cleaning up intermediate files..." | tee -a "$LOG_FILE"
    if [ "$INPUT_MODE" == "paired" ]; then
        rm -f "$R1_OUT" "$R2_OUT"
    else
        rm -f "$SE_OUT"
    fi
else
    echo "[INFO] Keeping intermediate files (--keep flag used)" | tee -a "$LOG_FILE"
fi

echo "[INFO] Analysis complete for sample $SAMPLE" | tee -a "$LOG_FILE"
echo "[INFO] Final results saved to: $FINAL_RESULTS" | tee -a "$LOG_FILE"
echo "[INFO] Summary results saved to: $SUMMARY_RESULTS" | tee -a "$LOG_FILE"
