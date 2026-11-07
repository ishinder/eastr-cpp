#pragma once

#include <htslib/sam.h>

#include <string>

namespace eastr {

// Set the reference FASTA used to decode CRAM inputs. Called once from main()
// with -r before any alignment file is opened. Without it, htslib can only
// decode CRAM via REF_PATH/REF_CACHE or by downloading the reference.
void set_alignment_reference(const std::string& fasta_path);

// sam_open(path, "r") for reading BAM/CRAM/SAM input. Applies the reference set
// by set_alignment_reference() so CRAM can be decoded (BAM/SAM ignore it).
// Returns nullptr on failure, like sam_open.
samFile* open_alignment(const std::string& path);

} // namespace eastr
