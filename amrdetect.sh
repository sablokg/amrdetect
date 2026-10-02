#!/usr/bin/env bash
#
# amr_pipeline.sh — AMR gene detection from paired-end short reads Gaurav Sablok gsablok@proton.me
#
# Usage: ./amr_pipeline.sh -1 R1.fastq.gz -2 R2.fastq.gz -o outdir -s sample_name [-t threads]
#
# Requires (install via conda/mamba):
#   fastp, spades.py, quast.py, amrfinder (+ amrfinder_update), abricate
#
# conda create -n amr-pipeline -c bioconda -c conda-forge \
#   fastp spades quast ncbi-amrfinderplus abricate -y

set -euo pipefail
IFS=$'\n\t'

# ---------------------------
# Defaults
# ---------------------------
THREADS=8
OUTDIR="amr_results"
SAMPLE="sample"
R1=""
R2=""
ASSEMBLY=""     # optional: skip assembly if you already have a FASTA

# ---------------------------
# Argument parsing
# ---------------------------
usage() {
    cat <<EOF
Usage: $0 -1 R1.fastq.gz -2 R2.fastq.gz -s SAMPLE_NAME [-o OUTDIR] [-t THREADS]
       $0 -a assembly.fasta -s SAMPLE_NAME [-o OUTDIR]   # skip assembly step

  -1  Forward reads (fastq.gz)
  -2  Reverse reads (fastq.gz)
  -a  Pre-assembled genome FASTA (skips QC+assembly steps)
  -s  Sample name (required)
  -o  Output directory (default: amr_results)
  -t  Threads (default: 8)
  -h  Show this help
EOF
    exit 1
}

while getopts "1:2:a:s:o:t:h" opt; do
    case $opt in
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        a) ASSEMBLY="$OPTARG" ;;
        s) SAMPLE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -z "$SAMPLE" ]] && { echo "ERROR: sample name (-s) required"; usage; }
if [[ -z "$ASSEMBLY" && ( -z "$R1" || -z "$R2" ) ]]; then
    echo "ERROR: provide either -1/-2 reads or -a assembly"
    usage
fi

mkdir -p "$OUTDIR"/{qc,assembly,quast,amrfinder,abricate,logs,report}
LOG="$OUTDIR/logs/${SAMPLE}.log"
exec > >(tee -a "$LOG") 2>&1

echo "=== AMR Pipeline started: $(date) ==="
echo "Sample: $SAMPLE | Threads: $THREADS | Outdir: $OUTDIR"

# ---------------------------
# Dependency check
# ---------------------------
check_tool() {
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH. Install it first."; exit 1; }
}

if [[ -z "$ASSEMBLY" ]]; then
    for t in fastp spades.py quast.py; do check_tool "$t"; done
fi
check_tool amrfinder
check_tool abricate

# ---------------------------
# 1. Read QC & trimming
# ---------------------------
if [[ -z "$ASSEMBLY" ]]; then
    echo "--- Step 1: Read QC (fastp) ---"
    TRIM_R1="$OUTDIR/qc/${SAMPLE}_R1.trim.fastq.gz"
    TRIM_R2="$OUTDIR/qc/${SAMPLE}_R2.trim.fastq.gz"

    fastp \
        -i "$R1" -I "$R2" \
        -o "$TRIM_R1" -O "$TRIM_R2" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread "$THREADS" \
        --json "$OUTDIR/qc/${SAMPLE}.fastp.json" \
        --html "$OUTDIR/qc/${SAMPLE}.fastp.html"

    # ---------------------------
    # 2. Assembly
    # ---------------------------
    echo "--- Step 2: Genome assembly (SPAdes) ---"
    SPADES_DIR="$OUTDIR/assembly/${SAMPLE}_spades"

    spades.py \
        -1 "$TRIM_R1" -2 "$TRIM_R2" \
        --isolate \
        -o "$SPADES_DIR" \
        -t "$THREADS" \
        -m 32

    ASSEMBLY="$OUTDIR/assembly/${SAMPLE}.fasta"
    # Filter short/low-coverage contigs (<500bp)
    seqkit seq -m 500 "$SPADES_DIR/contigs.fasta" > "$ASSEMBLY" 2>/dev/null \
        || cp "$SPADES_DIR/contigs.fasta" "$ASSEMBLY"

    echo "--- Step 3: Assembly QC (QUAST) ---"
    quast.py "$ASSEMBLY" \
        -o "$OUTDIR/quast/${SAMPLE}" \
        -t "$THREADS" \
        --min-contig 200
else
    echo "--- Skipping QC/assembly: using provided assembly $ASSEMBLY ---"
fi

# ---------------------------
# 4. AMR detection — AMRFinderPlus (NCBI, curated reference DB)
# ---------------------------
echo "--- Step 4: AMR detection (AMRFinderPlus) ---"
amrfinder \
    -n "$ASSEMBLY" \
    --threads "$THREADS" \
    --plus \
    -o "$OUTDIR/amrfinder/${SAMPLE}_amrfinder.tsv" \
    || echo "WARNING: amrfinder failed — check organism/DB setup (run 'amrfinder_update' first)"

# ---------------------------
# 5. AMR detection — ABRicate (cross-validation with multiple DBs)
# ---------------------------
echo "--- Step 5: AMR detection (ABRicate, multi-database) ---"
for DB in card resfinder ncbi argannot; do
    echo "  -> database: $DB"
    abricate --db "$DB" --threads "$THREADS" "$ASSEMBLY" \
        > "$OUTDIR/abricate/${SAMPLE}_${DB}.tsv" \
        || echo "  WARNING: abricate DB '$DB' not installed (abricate --setupdb)"
done

# Summary across abricate DBs
abricate --summary "$OUTDIR"/abricate/${SAMPLE}_*.tsv \
    > "$OUTDIR/abricate/${SAMPLE}_summary.tsv" 2>/dev/null || true

# ---------------------------
# 6. Consolidated report
# ---------------------------
echo "--- Step 6: Building consolidated report ---"
REPORT="$OUTDIR/report/${SAMPLE}_AMR_report.tsv"

{
    echo -e "sample\tsource_tool\tgene\tcoverage\tidentity\tdrug_class\tcontig"
    if [[ -f "$OUTDIR/amrfinder/${SAMPLE}_amrfinder.tsv" ]]; then
        tail -n +2 "$OUTDIR/amrfinder/${SAMPLE}_amrfinder.tsv" | \
        awk -F'\t' -v s="$SAMPLE" 'BEGIN{OFS="\t"} {print s,"AMRFinderPlus",$6,$16,$17,$11,$1}'
    fi
    for f in "$OUTDIR"/abricate/${SAMPLE}_*.tsv; do
        [[ "$f" == *summary* ]] && continue
        db=$(basename "$f" .tsv | sed "s/${SAMPLE}_//")
        tail -n +2 "$f" 2>/dev/null | \
        awk -F'\t' -v s="$SAMPLE" -v db="$db" 'BEGIN{OFS="\t"} {print s,"ABRicate:"db,$6,$10,$11,$14,$2}'
    done
} > "$REPORT"

echo "=== AMR Pipeline completed: $(date) ==="
echo "Final report: $REPORT"
echo "Assembly:     $ASSEMBLY"
echo "Log:          $LOG"