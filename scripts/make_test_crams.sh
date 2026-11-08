#!/usr/bin/env bash
# Build CRAM copies of the chr4 test BAMs (test_data/chr4_subset/*.bam) in
# test_data/chr4_subset_cram/, for testing CRAM input end to end.
#
# The CRAMs are encoded against a temporary uncompressed copy of the test
# reference, which is deleted afterwards. Their @SQ UR: tags therefore point
# to a path that no longer exists, like a CRAM made on another machine, so
# EASTR can only decode them with the reference passed via -r.
#
# Requires samtools. Run from the repository root:
#   scripts/make_test_crams.sh
set -euo pipefail

src_dir=test_data/chr4_subset
out_dir=test_data/chr4_subset_cram
ref_gz=test_data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

command -v samtools >/dev/null || { echo "samtools is required" >&2; exit 1; }
[ -f "$ref_gz" ] || { echo "missing $ref_gz (run 'git lfs pull')" >&2; exit 1; }

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
gunzip -c "$ref_gz" > "$tmp/reference.fa"
samtools faidx "$tmp/reference.fa"

mkdir -p "$out_dir"
for bam in "$src_dir"/*.bam; do
    cram="$out_dir/$(basename "${bam%.bam}").cram"
    echo "Writing $cram"
    samtools view -C -T "$tmp/reference.fa" -o "$cram" "$bam"
    samtools index "$cram"
done
echo "Done. The reference used for encoding has been removed."
