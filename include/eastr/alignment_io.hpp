#pragma once

#include <htslib/sam.h>

#include <string>

namespace eastr {

// Set the reference FASTA used to decode CRAM inputs. Called once from main()
// with -r before any alignment file is opened. Without it, htslib can only
// decode CRAM via REF_PATH/REF_CACHE or by downloading the reference.
void set_alignment_reference(const std::string& fasta_path);

// sam_open(path, "r") for reading BAM/CRAM/SAM input. Applies the reference set
// by set_alignment_reference() so CRAM can be decoded (BAM/SAM ignore it), and
// attaches `threads` htslib decompression threads when threads > 1.
//
// CRAM always gets an htslib thread pool (at least 1 thread), so that CRAM
// slices are decoded on htslib's own worker threads, which have >= 3 MiB
// stacks. htscodecs' NEON rANS-Nx16 order-1 decoder
// (rans_uncompress_O1_32x16_neon) keeps ~1.25 MiB of tables on the stack,
// which overflows the 512 KiB default stack of secondary threads on macOS
// (std::thread / std::async), corrupting adjacent thread memory. Do not call
// hts_set_threads() again on the returned handle: htslib ignores a second
// pool for CRAM.
// Returns nullptr on failure, like sam_open.
samFile* open_alignment(const std::string& path, int threads = 1);

} // namespace eastr
