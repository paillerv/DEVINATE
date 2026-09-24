#!/usr/bin/env bash

set -euo pipefail

# -------------------------
# USAGE
# -------------------------
if [ "$#" -lt 4 ]; then
    echo "Usage:"
    echo "  $0 <reads.fasta.gz> <TE.fasta> <genome.fa> <prefix>"
    exit 1
fi

READS=$1
TE=$2
GENOME=$3
PREFIX=$4

# Chrono - Début
START_TIME=$(date +%s)

# -------------------------
# WORKDIR
# -------------------------
WORKDIR="${PREFIX}"
SAMPLE_NAME=$(basename "$PREFIX")
OUTDIR="${WORKDIR}"/"${SAMPLE_NAME}_classification"

mkdir -p "$WORKDIR"
mkdir -p "$OUTDIR" 

# Extract TE size with python Bio for the #11 module
TE_SIZE=$(python3 -c "import sys; from Bio import SeqIO; print(len(next(SeqIO.parse('$TE', 'fasta'))))")
TE_MIN=$((TE_SIZE - 150))
TE_MAX=$((TE_SIZE + 150))

# Remove TE extension
TE_BASENAME=$(basename "$TE")
TE_NAME="${TE_BASENAME%.*}"

# Extract script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Threads 
THREADS=32

echo -e "\n" 

# -------------------------
# FILES output
# -------------------------
READS_FA="${WORKDIR}/${SAMPLE_NAME}.reads.fasta.gz"

TE_BAM="${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.vs_TE.raw.bam"
TE_BAM_FILT="${WORKDIR}/${SAMPLE_NAME}.all_reads_vs_${TE_NAME}.nosupp.q20.bam"

CLASSIF="${OUTDIR}/classification.tsv"

READS_KEEP="${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.KEEP_readsIds.txt"
READS_KEEP_FA="${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.KEEP_reads.fasta"

GENOME_BAM="${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.KEEP_reads_vs_genome.bam"
GENOME_BAM_FILT="${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.KEEP_reads_vs_genome.nosupp.q20.bam"

ALL_GENOME_BAM="${WORKDIR}/${SAMPLE_NAME}.all_reads_vs_genome.bam"
ALL_GENOME_BAM_FILT="${WORKDIR}/${SAMPLE_NAME}.all_reads_vs_genome.nosupp.q20.bam"

FLANKS_FA="${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.flanks.fasta"
FLANKS_BAM="${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.flanks_vs_genome.bam"

INSERTIONS="${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.insertions.bed"
CLUSTERS="${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.clusters.tsv"

# -------------------------
# 0. REMOVE CHIMERAS
# -------------------------
echo "[STEP 0] Remove chimerical reads"

seqkit grep -r -v -p ";" "$READS" \
    -o "$READS_FA"

echo -e "\n" 

# -------------------------
# 1. Mapping reads on TE
# -------------------------
echo "[STEP 1] Mapping reads -> TE"

minimap2 -ax map-ont --secondary=no -t $THREADS "$TE" "$READS_FA" \
    | samtools sort -@ $THREADS -o "$TE_BAM"

echo -e "\n" 

# -------------------------
# 2. Filter alignments
# -------------------------
echo "[STEP 2] Filtering TE alignments"

samtools view -F 4 -F 2048 -q 20 "$TE_BAM" -b > "$TE_BAM_FILT"
samtools index "$TE_BAM_FILT"

rm -f "$TE_BAM"

echo -e "\n" 

# -------------------------
# 2.5.1 AUTOMATIC VARIANT DISCOVERY
# -------------------------
echo "[STEP 2.5.1] Automatic Variant Discovery"
VARIANTS_FILE="${OUTDIR}/auto_discovered_variants.txt"

python3 "$SCRIPT_DIR"/scripts/discover_variants_only.py \
    "$TE_BAM_FILT" \
    "$VARIANTS_FILE"

# Generate a BED variants file for bamdash
awk '
{
    if ($0 ~ /:DEL:/) {
        split($0, a, ":DEL:")
    } else if ($0 ~ /:INS:/) {
        split($0, a, ":INS:")
    } else {
        next
    }

    name = a[1]
    split(a[2], b, ":")
    split(b[1], pos, "-")

    print name "\t" pos[1] "\t" pos[2]
}
' "$VARIANTS_FILE" > "${OUTDIR}/${SAMPLE_NAME}.variants.bed" 

echo -e "\n" 

# ---------------------------
# 2.5.2 PLOT COVERAGE ALONG TE WITH BAMDASH
# ---------------------------
echo "[STEP 2.5.2] Plot coverage along the TE\n"

bamdash --verbose \
   -b "${TE_BAM_FILT}" \
   -q 0 \
   -p "${OUTDIR}/${SAMPLE_NAME}.coverage" \
   -t "${OUTDIR}/${SAMPLE_NAME}.variants.bed" 


echo "The HTML coverage file is generated in ${OUTDIR}/${SAMPLE_NAME}.coverage.html is generated." 

echo -e "\n" 

# -------------------------
# 3. Classification
# -------------------------
echo "[STEP 3] Reads classification"

python3 "$SCRIPT_DIR"/scripts/classify_variants.py \
    "$TE_BAM_FILT" \
    "$VARIANTS_FILE" \
    "$OUTDIR"

echo -e "\n" 

# -------------------------
# 4. Extract KEEP reads
# -------------------------
echo "[STEP 4] Extract KEEP-tagged reads"

awk 'NR>1 && $12=="KEEP"{print $1}' \
    "$CLASSIF" \
    | sort -u > "$READS_KEEP"

N=$(wc -l < "$READS_KEEP")

echo "[INFO] KEEP reads: $N"

if [ "$N" -eq 0 ]; then
    echo "[ERROR] No KEEP reads found"
    exit 1
fi

seqkit grep -f "$READS_KEEP" \
    "$READS_FA" \
    -o "$READS_KEEP_FA"

echo -e "\n" 

# -------------------------
# 5. Mapping reads -> genome
# -------------------------
echo "[STEP 5] Mapping reads -> genome"

minimap2 -ax map-ont -t $THREADS \
    "$GENOME" "$READS_KEEP_FA" \
    | samtools sort -@ $THREADS -o "$GENOME_BAM"

samtools view -F 4 -F 2048 -q 20 \
    "$GENOME_BAM" -b > "$GENOME_BAM_FILT"

samtools index "$GENOME_BAM_FILT"

rm -f "$GENOME_BAM"

echo -e "\n" 

# -------------------------
# 6. Extract flanks
# -------------------------
echo "[STEP 6] Extract flanks sequences"

python3 "$SCRIPT_DIR"/scripts/extract_flanks.py \
    "$TE_BAM_FILT" \
    "$CLASSIF" \
    "$FLANKS_FA" \
    --min-flank 200

echo -e "\n" 

# -------------------------
# 7. Mapping flanks on genome
# -------------------------
echo "[STEP 7] Mapping flanks on dm6"

minimap2 -ax map-ont --secondary=no -t $THREADS \
    "$GENOME" "$FLANKS_FA" \
    | samtools sort -@ $THREADS -o "$FLANKS_BAM"

samtools index "$FLANKS_BAM"

echo -e "\n" 

# -------------------------
# 8. Filter flanks + pairing
# -------------------------
echo "[STEP 8] Filtering flanks + pairing"

python3 "$SCRIPT_DIR"/scripts/filter_flanking_pairs.py \
    "$FLANKS_BAM" > "$INSERTIONS"

echo -e "\n" 

# -------------------------
# 9. Clustering
# -------------------------
echo "[STEP 9] Clustering insertions"

python3 "$SCRIPT_DIR"/scripts/cluster.py \
    "$INSERTIONS" > "$CLUSTERS"

echo -e "\n" 

# -------------------------
# STEP 10. Global mapping of all raw reads on genome
# -------------------------
echo "[STEP 10] Global mapping of all raw reads on genome"

minimap2 -ax map-ont -t $THREADS "$GENOME" "$READS_FA" \
    | samtools sort -@ $THREADS -o "$ALL_GENOME_BAM"

samtools view -F 4 -F 2048 -q 20 "$ALL_GENOME_BAM" -b > "$ALL_GENOME_BAM_FILT"
samtools index "$ALL_GENOME_BAM_FILT"

echo -e "\n" 

# -------------------------
# STEP 11. Evaluate Insertion Sites
# -------------------------
echo "[STEP 11] Evaluating insertion sites status"

FINAL_OUTPUT="${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.empty_sites.VAF_freq.tsv"

python3 "$SCRIPT_DIR"/scripts/evaluate_empty_sites.py \
    "$CLUSTERS" \
    "$INSERTIONS" \
    "$CLASSIF" \
    "$ALL_GENOME_BAM_FILT" \
    "$FINAL_OUTPUT"

echo -e "\n"

# -------------------------
# STEP 12. Mapping concordance flanks vs full sequence on dm6
# -------------------------
echo "[STEP 12] Mapping corcodance - Flanks vs Full sequence"

python3 "$SCRIPT_DIR"/scripts/compare_flanks_vs_fullreads.py \
    "$INSERTIONS" \
    "$ALL_GENOME_BAM_FILT" \
    "$CLASSIF" > "${WORKDIR}/${SAMPLE_NAME}.${TE_NAME}.compare_flanks_vs_fullreads.tsv"

# -------------------------
# 13. Cleanup
# -------------------------
echo "[STEP 13] Cleanup temporary files"

rm -f "$READS_FA" "$ALL_GENOME_BAM"

# -------------------------
# DONE
# -------------------------
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo ""
echo "[DONE] Pipeline finished in ${MINUTES} min and ${SECONDS} sec."
echo "Results:"
echo "  $WORKDIR"
