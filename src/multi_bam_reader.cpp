#include "eastr/multi_bam_reader.hpp"

#include <cstring>
#include <stdexcept>
#include <filesystem>

namespace fs = std::filesystem;

namespace eastr {

MultiBamReader::MultiBamReader(const std::vector<std::string>& bam_paths) {
    if (bam_paths.empty()) {
        throw std::runtime_error("No BAM files provided");
    }

    // Allocate owned copy for sync mode
    sync_current_aln_ = bam_init1();

    files_.reserve(bam_paths.size());
    records_.reserve(bam_paths.size());
    sample_names_.reserve(bam_paths.size());
    exhausted_.resize(bam_paths.size(), false);

    for (size_t i = 0; i < bam_paths.size(); ++i) {
        const auto& path = bam_paths[i];

        // Open BAM file
        samFile* fp = sam_open(path.c_str(), "r");
        if (!fp) {
            throw std::runtime_error("Cannot open BAM file: " + path);
        }

        // Read header
        sam_hdr_t* hdr = sam_hdr_read(fp);
        if (!hdr) {
            sam_close(fp);
            throw std::runtime_error("Cannot read BAM header: " + path);
        }

        // Validate sorted
        validate_sorted(fp, hdr, path);

        // Use first file's header as the shared header
        if (i == 0) {
            header_ = hdr;
        } else {
            sam_hdr_destroy(hdr);
        }

        files_.push_back(fp);

        // Allocate record buffer
        bam1_t* rec = bam_init1();
        records_.push_back(rec);

        // Extract sample name from filename
        fs::path p(path);
        sample_names_.push_back(p.stem().string());

        // Read first record and add to queue
        if (read_next(static_cast<int>(i))) {
            queue_.push({static_cast<int>(i), records_[i]->core.tid,
                        records_[i]->core.pos});
        }
    }
}

MultiBamReader::~MultiBamReader() {
    // Stop prefetch threads if running
    disable_prefetching();

    for (auto* rec : records_) {
        if (rec) bam_destroy1(rec);
    }
    for (auto* fp : files_) {
        if (fp) sam_close(fp);
    }
    if (header_) {
        sam_hdr_destroy(header_);
    }
    if (sync_current_aln_) {
        bam_destroy1(sync_current_aln_);
    }
    if (prefetch_current_aln_) {
        bam_destroy1(prefetch_current_aln_);
    }
}

void MultiBamReader::validate_sorted(samFile* fp, sam_hdr_t* hdr,
                                      const std::string& path) {
    kstring_t val = {0, 0, nullptr};

    int ret = sam_hdr_find_tag_hd(hdr, "SO", &val);
    if (ret < 0 || !val.s || strcmp(val.s, "coordinate") != 0) {
        ks_free(&val);
        throw std::runtime_error(
            "BAM file must be coordinate-sorted: " + path +
            "\nRun: samtools sort -o sorted.bam input.bam");
    }

    ks_free(&val);
}

bool MultiBamReader::read_next(int file_idx) {
    if (exhausted_[file_idx]) {
        return false;
    }

    int ret = sam_read1(files_[file_idx], header_, records_[file_idx]);
    if (ret < 0) {
        exhausted_[file_idx] = true;
        return false;
    }

    return true;
}

MultiBamReader::AlignmentRecord* MultiBamReader::next_sync() {
    if (queue_.empty()) {
        return nullptr;
    }

    // Get the file with the smallest coordinate
    QueueEntry entry = queue_.top();
    queue_.pop();

    int idx = entry.file_idx;

    // Deep copy the record to our owned buffer BEFORE reading next
    // (read_next overwrites records_[idx])
    bam_copy1(sync_current_aln_, records_[idx]);

    // Set up return value with our owned copy
    current_record_.aln = sync_current_aln_;
    current_record_.file_index = idx;
    current_record_.sample_name = sample_names_[idx];
    current_record_.tid = entry.tid;
    current_record_.pos = entry.pos;

    // Read next record from this file and re-add to queue
    if (read_next(idx)) {
        queue_.push({idx, records_[idx]->core.tid, records_[idx]->core.pos});
    }

    return &current_record_;
}

MultiBamReader::AlignmentRecord* MultiBamReader::next() {
    if (!prefetching_enabled_) {
        return next_sync();
    }
    return next_prefetch();
}

void MultiBamReader::enable_prefetching(size_t buffer_size) {
    if (prefetching_enabled_) {
        return;  // Already enabled
    }

    size_t num_files = files_.size();
    per_file_buffer_size_ = std::max(size_t(100), buffer_size / num_files);

    // Allocate the current record alignment for prefetch mode
    prefetch_current_aln_ = bam_init1();

    // Initialize per-file states
    per_file_states_.resize(num_files);
    current_from_file_.resize(num_files);

    for (size_t i = 0; i < num_files; ++i) {
        per_file_states_[i] = std::make_unique<PerFileState>();
        auto& state = *per_file_states_[i];
        state.buffer_size = per_file_buffer_size_;
        state.buffer.resize(per_file_buffer_size_);

        // Initialize buffer records
        for (auto& rec : state.buffer) {
            rec.aln = bam_init1();
            rec.valid = false;
        }

        // Initialize current record from this file
        current_from_file_[i].aln = bam_init1();
        current_from_file_[i].valid = false;
    }

    prefetching_enabled_ = true;

    // Clear the sync-mode queue (we'll rebuild using prefetch buffers)
    while (!queue_.empty()) queue_.pop();

    // Start per-file prefetch threads
    for (size_t i = 0; i < num_files; ++i) {
        per_file_states_[i]->thread = std::thread(
            &MultiBamReader::per_file_prefetch_worker, this, static_cast<int>(i));
    }

    // Prime the queue: pop first record from each file's buffer
    for (size_t i = 0; i < num_files; ++i) {
        if (refill_queue_entry(static_cast<int>(i))) {
            queue_.push({static_cast<int>(i),
                        current_from_file_[i].tid,
                        current_from_file_[i].pos});
        }
    }
}

void MultiBamReader::disable_prefetching() {
    if (!prefetching_enabled_) {
        return;
    }

    // Signal all threads to stop
    for (auto& state : per_file_states_) {
        if (state) {
            state->stop = true;
            state->not_full.notify_all();
        }
    }

    // Wait for all threads to finish
    for (auto& state : per_file_states_) {
        if (state && state->thread.joinable()) {
            state->thread.join();
        }
    }

    // Clean up buffers
    for (auto& state : per_file_states_) {
        if (state) {
            for (auto& rec : state->buffer) {
                if (rec.aln) {
                    bam_destroy1(rec.aln);
                    rec.aln = nullptr;
                }
            }
        }
    }

    for (auto& rec : current_from_file_) {
        if (rec.aln) {
            bam_destroy1(rec.aln);
            rec.aln = nullptr;
        }
    }

    per_file_states_.clear();
    current_from_file_.clear();

    if (prefetch_current_aln_) {
        bam_destroy1(prefetch_current_aln_);
        prefetch_current_aln_ = nullptr;
    }

    prefetching_enabled_ = false;
}

void MultiBamReader::per_file_prefetch_worker(int file_idx) {
    auto& state = *per_file_states_[file_idx];
    samFile* fp = files_[file_idx];
    bam1_t* temp_rec = bam_init1();

    while (!state.stop) {
        // Wait for space in buffer
        size_t next_write = (state.write_pos.load() + 1) % state.buffer_size;

        std::unique_lock<std::mutex> lock(state.mutex);
        state.not_full.wait(lock, [&] {
            return next_write != state.read_pos.load() || state.stop.load();
        });

        if (state.stop) {
            break;
        }
        lock.unlock();

        // Read next record from file
        int ret = sam_read1(fp, header_, temp_rec);
        if (ret < 0) {
            // EOF or error
            state.eof = true;
            state.not_empty.notify_all();
            break;
        }

        // Copy to buffer
        auto& buf_rec = state.buffer[state.write_pos.load()];
        bam_copy1(buf_rec.aln, temp_rec);
        buf_rec.file_index = file_idx;
        buf_rec.sample_name = sample_names_[file_idx];
        buf_rec.tid = temp_rec->core.tid;
        buf_rec.pos = temp_rec->core.pos;
        buf_rec.valid = true;

        // Advance write position
        state.write_pos = next_write;
        state.not_empty.notify_one();
    }

    bam_destroy1(temp_rec);
}

bool MultiBamReader::pop_from_file_buffer(int file_idx, BufferedRecord& out) {
    auto& state = *per_file_states_[file_idx];

    std::unique_lock<std::mutex> lock(state.mutex);
    state.not_empty.wait(lock, [&] {
        return state.read_pos.load() != state.write_pos.load() || state.eof.load();
    });

    // Check for EOF with empty buffer
    if (state.read_pos.load() == state.write_pos.load() && state.eof.load()) {
        return false;
    }

    // Copy from buffer
    auto& buf_rec = state.buffer[state.read_pos.load()];
    bam_copy1(out.aln, buf_rec.aln);
    out.file_index = buf_rec.file_index;
    out.sample_name = buf_rec.sample_name;
    out.tid = buf_rec.tid;
    out.pos = buf_rec.pos;
    out.valid = true;

    // Advance read position
    state.read_pos = (state.read_pos.load() + 1) % state.buffer_size;
    lock.unlock();
    state.not_full.notify_one();

    return true;
}

bool MultiBamReader::refill_queue_entry(int file_idx) {
    return pop_from_file_buffer(file_idx, current_from_file_[file_idx]);
}

MultiBamReader::AlignmentRecord* MultiBamReader::next_prefetch() {
    if (queue_.empty()) {
        return nullptr;
    }

    // Get the file with the smallest coordinate
    QueueEntry entry = queue_.top();
    queue_.pop();

    int idx = entry.file_idx;

    // Copy from current_from_file_ to our return value
    auto& src = current_from_file_[idx];
    bam_copy1(prefetch_current_aln_, src.aln);

    current_record_.aln = prefetch_current_aln_;
    current_record_.file_index = src.file_index;
    current_record_.sample_name = src.sample_name;
    current_record_.tid = src.tid;
    current_record_.pos = src.pos;

    // Refill from this file's buffer and re-add to queue if successful
    if (refill_queue_entry(idx)) {
        queue_.push({idx, current_from_file_[idx].tid, current_from_file_[idx].pos});
    }

    return &current_record_;
}

size_t MultiBamReader::buffer_occupancy() const {
    if (!prefetching_enabled_) {
        return 0;
    }

    size_t total = 0;
    for (const auto& state : per_file_states_) {
        if (state) {
            size_t wp = state->write_pos.load();
            size_t rp = state->read_pos.load();
            if (wp >= rp) {
                total += wp - rp;
            } else {
                total += state->buffer_size - rp + wp;
            }
        }
    }
    return total;
}

} // namespace eastr
