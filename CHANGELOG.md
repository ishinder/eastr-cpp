# Changelog

All notable changes to this project are documented here. This project adheres to
[Semantic Versioning](https://semver.org/).

## [2.1.1] - unreleased

### Fixed
- `--out_filtered_bam <file>` with a single BAM and a bare filename (no directory
  prefix, e.g. `--out_filtered_bam filtered.bam`) created a **directory** named after
  the output (`filtered.bam/`) and wrote the filtered BAM *inside* it, instead of
  writing the file. 2.1.0 corrected output-path generation but missed two inline
  sites that run first and create the directory: the auto-built bowtie2 index
  directory and the temporary spurious-junction BED directory both treated a bare
  output filename as a directory. Creating that directory also flipped the
  file/directory heuristic (which checks `fs::is_directory`), defeating the 2.1.0
  change. Now a bare single-file `--out_filtered_bam` writes exactly that file (and
  `--removed_alignments_bam` writes `<name>_removed_alignments.bam` next to it), with
  no stray directory. For a bare filename the auto-built bowtie2 index is written to
  its own `<reference>_bt2_idx/` directory next to the reference FASTA (matching the
  default when no output path is given).
- The temporary spurious-junction BED now uses a process-unique name, so concurrent
  EASTR runs sharing an output directory no longer clobber each other's temp file.

## [2.1.0] - 2026-06-08

### Added
- `--version` flag that prints the version and exits 0. The version is generated
  from the CMake `project(... VERSION ...)` (a single source of truth) via a
  configured `eastr/version.hpp`, so the git tag, CMake version, `--version`
  output, and conda package can no longer drift. A CI check fails on any drift
  between `CMakeLists.txt` and `conda/meta.yaml`.
- `--bed_list` / `--bam_list` flags to explicitly declare that the `--bed` / `--bam`
  argument is a text file listing input paths (one per line).
- GitHub Actions CI building and testing on `ubuntu-latest` (x86_64),
  `ubuntu-24.04-arm` (linux-aarch64), and `macos-latest` (osx-arm64).
- Unit tests for output path generation and content-based input detection,
  including a Galaxy-style `*.dat`-named input case and a single-file
  `--out_original_junctions` case.

### Changed
- Input-type detection no longer depends on the file extension. `--bam` is detected
  by magic bytes (BGZF/BAM, CRAM, SAM); `--bed` is detected by content (a BED record
  vs. a list of existing paths). This fixes inputs named `*.dat` (e.g. every Galaxy
  dataset). Auto-detecting a path list still works but now emits a deprecation
  warning recommending `--bed_list` / `--bam_list`. `--gtf` already accepted any
  filename. (Issue #2)
- `--out_original_junctions` now accepts a **file path** for a single input (BAM,
  GTF, or BED) and writes one file, consistent with `--out_removed_junctions` /
  `--out_kept_junctions`. It previously behaved as a directory for BAM input and was
  a silent no-op for GTF/BED. Directory mode is still used for multi-input runs.
  (Issues #1, #3)
- `--out_filtered_bam` accepts a **file path** for a single BAM input (writes the
  filtered BAM directly), keeping directory mode for BAM lists. (Issue #4) **Note:**
  this was only partially fixed in 2.1.0 — a bare filename still created a directory;
  fully fixed in 2.1.1.

### Fixed
- Bare output filenames without a directory prefix (e.g. `out.bed`, `out.bam`) were
  misinterpreted as directories, causing "Cannot open output file" errors. Both bare
  filenames and `dir/out.bed` now work. (Issues #1, #3, #4)
- `--out_filtered_bam` with a file-looking path and multiple BAM inputs previously
  produced paths like `out.bam/<sample>_EASTR_filtered.bam`; it now falls back to the
  parent directory.

### Notes
- `conda/meta.yaml` `sha256` for the source tarball must be updated when the `v2.1.0`
  tag is cut.
- Adding `additional_platforms: [linux-aarch64]` to the Bioconda recipe is a separate
  PR in `bioconda-recipes`.
