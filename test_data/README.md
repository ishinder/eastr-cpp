# Test Data

This directory contains test data for EASTR C++. The test dataset is a subset of Arabidopsis thaliana chromosome 4 RNA-seq alignments.

## Contents

### `chr4_subset/`
Three sample BAM files containing RNA-seq alignments for chromosome 4:
- `SRR14056780_chr4.bam` - Sample 1 (~113MB)
- `SRR14056781_chr4.bam` - Sample 2 (~154MB)
- `SRR14056782_chr4.bam` - Sample 3 (~114MB)

Each BAM file includes its corresponding `.bam.bai` index.

### `reference/`
Arabidopsis thaliana TAIR10 reference genome (compressed):
- `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz` - Reference FASTA (~35MB)
- `Arabidopsis_thaliana.TAIR10.58.gtf.gz` - Gene annotations (~11MB)

### `indices/bowtie2/`
Pre-built Bowtie2 index for the reference genome (~170MB total):
- `arabidopsis_ref.{1,2,3,4}.bt2`
- `arabidopsis_ref.rev.{1,2}.bt2`

### `bam_list.txt`
A ready-to-use BAM list file with relative paths for testing EASTR.

## Usage

Run from the repository root directory:

```bash
# Using the BAM list
./build/eastr \
    --bam test_data/bam_list.txt \
    --reference test_data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
    --bowtie2_index test_data/indices/bowtie2/arabidopsis_ref \
    --out_filtered_bam test_output/ \
    --out_removed_junctions test_output/spurious.bed \
    --verbose

# Using the GTF
./build/eastr \
    --gtf test_data/reference/Arabidopsis_thaliana.TAIR10.58.gtf.gz \
    --reference test_data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
    --bowtie2_index test_data/indices/bowtie2/arabidopsis_ref \
    --out_removed_junctions test_output/spurious.bed \
    --verbose
```

## Data Source

- **Reference genome**: Ensembl Plants, Arabidopsis thaliana TAIR10
- **BAM files**: SRA accession SRR14056780, SRR14056781, SRR14056782 (chr4 subset)

## Storage Notes

Large files in this directory are tracked using Git LFS. After cloning, run:
```bash
git lfs pull
```
to download the actual file contents.
