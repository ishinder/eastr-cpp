#!/bin/bash
set -euo pipefail

# Build hg38 bowtie2 index with pseudoautosomal regions (PARs) masked
#
# The PAR regions on chrY are identical to chrX, causing multi-mapping issues
# in RNA-seq analysis. This script masks them with N's before building the index.
#
# PAR coordinates (GRCh38/hg38):
#   PAR1: chrY:10,001-2,781,479
#   PAR2: chrY:56,887,903-57,217,415
#
# Usage: ./build_hg38_noPARs_index.sh [output_dir] [threads]
#
# Requirements: wget, samtools, bedtools, bowtie2-build

OUTPUT_DIR="${1:-.}"
THREADS="${2:-8}"

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "=== Building hg38 bowtie2 index with PAR regions masked ==="
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo ""

# Step 1: Download hg38 reference if not present
if [[ ! -f hg38.fa && ! -f hg38.fa.gz ]]; then
    echo "Downloading hg38 reference from UCSC..."
    wget -q --show-progress https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
fi

if [[ -f hg38.fa.gz && ! -f hg38.fa ]]; then
    echo "Decompressing hg38.fa.gz..."
    gunzip -k hg38.fa.gz
fi

# Step 2: Index the reference
if [[ ! -f hg38.fa.fai ]]; then
    echo "Indexing reference..."
    samtools faidx hg38.fa
fi

# Step 3: Create BED file with PAR regions to mask
echo "Creating PAR mask BED file..."
cat > par_regions.bed << 'EOF'
chrY	10000	2781479	PAR1
chrY	56887902	57217415	PAR2
EOF

# Step 4: Mask PAR regions with N's
echo "Masking PAR regions on chrY..."
bedtools maskfasta -fi hg38.fa -bed par_regions.bed -fo hg38_noPARs.fa

# Step 5: Index the masked reference
echo "Indexing masked reference..."
samtools faidx hg38_noPARs.fa

# Step 6: Build bowtie2 index
echo "Building bowtie2 index (this may take 1-2 hours)..."
bowtie2-build --threads "$THREADS" hg38_noPARs.fa hg38_noPARs

# Cleanup intermediate files (optional - comment out to keep)
# rm -f hg38.fa hg38.fa.gz hg38.fa.fai par_regions.bed

echo ""
echo "=== Done! ==="
echo "Bowtie2 index created: $OUTPUT_DIR/hg38_noPARs"
echo ""
echo "Use with EASTR:"
echo "  eastr --bam input.bam -r $OUTPUT_DIR/hg38_noPARs.fa -i $OUTPUT_DIR/hg38_noPARs ..."
