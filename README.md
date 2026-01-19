# EASTR

<p align="center">
(\(\<br>
(-.-)<br>
o('')('')<br>
</p>

<p align="center">
<em>Emending Alignments of Spliced Transcript Reads</em>
</p>

<p align="center">
<a href="https://github.com/ishinder/eastr-cpp/blob/main/LICENSE" target="_blank">
    <img src="https://img.shields.io/github/license/ishinder/eastr-cpp" alt="License">
</a>
<a href="https://github.com/ishinder/eastr-cpp/releases" target="_blank">
    <img src="https://img.shields.io/github/v/release/ishinder/eastr-cpp" alt="Release">
</a>
</p>

---

EASTR is a tool for detecting and eliminating spuriously spliced alignments in RNA-seq datasets. It improves the accuracy of transcriptome assembly by identifying and removing misaligned spliced alignments. The tool can process GTF, BED, and BAM files as input and can be applied to any RNA-seq dataset regardless of the alignment software used.

## Installation

### Conda (recommended)

```bash
conda install -c bioconda eastr-cpp
```

<details>
<summary><b>Build from source</b></summary>

**Dependencies:**
- C++17 compiler (GCC >=7 or Clang >=5)
- CMake >=3.16
- htslib >=1.15
- bowtie2 >=2.5.2 (must be in PATH)

```bash
# Install dependencies
# macOS
brew install htslib bowtie2

# Ubuntu/Debian
sudo apt-get install libhts-dev bowtie2

# Conda environment
conda install -c bioconda htslib bowtie2

# Build EASTR
git clone --recurse-submodules https://github.com/ishinder/eastr-cpp.git
cd eastr-cpp
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Optional: install system-wide
sudo make install
```

</details>

## Quick Start

```bash
eastr --bam input.bam -r genome.fa --out_filtered_bam filtered.bam
```

EASTR automatically builds a bowtie2 index if one isn't provided. To use an existing index, add `-i /path/to/bt2_index`.

## Usage Examples

### Filter multiple BAM files in parallel

```bash
# Create a file listing your BAMs (one path per line)
ls /path/to/*.bam > samples.txt

# Filter all samples using 8 threads
eastr --bam samples.txt \
      -r genome.fa \
      --out_filtered_bam filtered_bams/ \
      -p 8 \
      --verbose
```

### Identify spurious junctions in a GTF annotation

```bash
eastr --gtf annotation.gtf \
      -r genome.fa \
      -i /path/to/bt2_index \
      --out_removed_junctions spurious_junctions.bed
```

> **Note:** The `-i` (bowtie2 index) flag is optional. If omitted, EASTR will automatically build an index from the reference FASTA, which adds a few minutes for large genomes. Providing a pre-built index speeds up repeated runs.

### Save removed alignments for inspection

```bash
eastr --bam input.bam \
      -r genome.fa \
      --out_filtered_bam filtered.bam \
      --removed_alignments_bam

# Output:
#   filtered.bam                        - alignments with spurious junctions removed
#   filtered_removed_alignments.bam     - the removed alignments
```

### Use trusted junctions that should never be removed

```bash
eastr --bam input.bam \
      -r genome.fa \
      --out_filtered_bam filtered.bam \
      --trusted_bed validated_junctions.bed
```

## Options

### Required (one input type)

| Option | Description |
|--------|-------------|
| `--bam FILE` | BAM file, or text file listing BAM paths (one per line) |
| `--gtf FILE` | GTF annotation file |
| `--bed FILE` | BED file with junction coordinates |
| `-r, --reference FILE` | Reference genome FASTA (required for all input types) |

### Common options

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --bowtie2_index PATH` | auto-built | Bowtie2 index prefix (built automatically if not provided) |
| `-p INT` | 1 | Number of threads |
| `--out_filtered_bam PATH` | — | Output filtered BAM file or directory |
| `--out_removed_junctions PATH` | stdout | Output spurious junctions (BED format) |
| `--verbose` | off | Show progress information |

<details>
<summary><b>Advanced options</b></summary>

### Algorithm parameters

| Option | Default | Description |
|--------|---------|-------------|
| `--bt2_k INT` | 10 | Minimum distinct bowtie2 alignments for spurious classification |
| `-o INT` | 50 | Flanking sequence length on each side of junction |
| `-a INT` | 7 | Minimum anchor length in each exon |
| `--min_duplicate_exon_length INT` | 27 | Minimum length for duplicated exon detection |
| `--min_junc_score INT` | 1 | Minimum supporting reads per junction |
| `--trusted_bed FILE` | — | BED file of junctions to never remove |

### Alignment scoring (minimap2 parameters)

| Option | Default | Description |
|--------|---------|-------------|
| `-A INT` | 3 | Match score |
| `-B INT` | 4 | Mismatch penalty |
| `-k INT` | 3 | K-mer length |
| `-w INT` | 2 | Minimizer window size |
| `-m INT` | 25 | Minimum chain score |

### Additional output options

| Option | Description |
|--------|-------------|
| `--out_original_junctions PATH` | Write all junctions before filtering |
| `--out_kept_junctions PATH` | Write non-spurious junctions |
| `--removed_alignments_bam` | Write removed alignments to separate BAM |
| `--filtered_bam_suffix STR` | Suffix for output BAMs (default: `_EASTR_filtered`) |

</details>

## Input requirements

- **BAM files** must be coordinate-sorted
- **Reference FASTA** can be uncompressed or bgzip-compressed (`.fa`, `.fa.gz`, `.fasta`, `.fasta.gz`)
- **GTF/BED files** can be plain text or gzipped

## Citation

If you use EASTR in your research, please cite:

> Shinder I, Hu R, Ji HJ, Chao KH, Pertea M. EASTR: Identifying and eliminating systematic alignment errors in multi-exon genes. *Nat Commun.* 2023;14:7223. doi: [10.1038/s41467-023-43017-4](https://doi.org/10.1038/s41467-023-43017-4)

## License

[MIT License](LICENSE)

---

<details>
<summary><b>Performance notes</b></summary>

This C++ implementation provides significant speedups over the [original Python version](https://github.com/ishinder/EASTR) through parallelization at every pipeline stage:

| Stage | Parallelization strategy |
|-------|--------------------------|
| BAM junction extraction | Inter-file + intra-file tiled parallelism |
| FASTA sequence retrieval | Multi-threaded with coordinate-sorted access |
| Self-alignment (minimap2) | Thread pool for parallel alignment |
| Bowtie2 validation | Native multi-threading |
| BAM filtering | Parallel files + multi-threaded BGZF compression |

### Benchmark

Test dataset: 3 Arabidopsis chr4 BAM files (~380MB total), 6 threads

| Version | Time | Speedup |
|---------|------|---------|
| C++ | 5.8s | **10x** |
| Python | 60s | — |

### Threading recommendations

The `-p` option controls parallelism across all stages. For best performance:

- Use `-p` equal to or greater than the number of input BAM files
- When `-p` exceeds the number of BAM files, EASTR uses tiled intra-file parallelism (requires BAM index files)
- For large datasets, `-p 8` to `-p 16` typically provides good throughput

</details>

<details>
<summary><b>Development</b></summary>

### Building with tests

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DEASTR_BUILD_TESTS=ON
make -j$(nproc)
ctest --output-on-failure
```

### Test data

Test data requires Git LFS:

```bash
git lfs install
git lfs pull
```

### Project structure

```
eastr-cpp/
├── include/eastr/     # Header files
├── src/               # Implementation
├── tests/             # Unit tests
├── external/          # Dependencies (minimap2)
├── conda/             # Conda recipe
└── test_data/         # Test datasets (Git LFS)
```

</details>
