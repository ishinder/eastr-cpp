#include "eastr/path_utils.hpp"

#include <filesystem>

namespace fs = std::filesystem;

namespace eastr {

bool is_file_path(const std::string& path) {
    if (path.empty()) {
        return false;
    }

    // An existing directory is unambiguously a directory.
    std::error_code ec;
    if (fs::is_directory(path, ec)) {
        return false;
    }

    // A trailing path separator means a directory was intended.
    char last = path.back();
    if (last == '/' || last == '\\') {
        return false;
    }

    // Otherwise: treat as a file if there is an extension after the last
    // separator. Guarding on last_slash avoids the npos pitfall where a bare
    // "out.bed" (no separator) would otherwise be misread as a directory.
    size_t last_dot = path.find_last_of('.');
    size_t last_slash = path.find_last_of("/\\");
    return last_dot != std::string::npos &&
           (last_slash == std::string::npos || last_dot > last_slash);
}

std::string path_basename(const std::string& path) {
    size_t last_slash = path.find_last_of("/\\");
    size_t last_dot = path.find_last_of('.');
    if (last_slash == std::string::npos) {
        last_slash = 0;
    } else {
        last_slash++;
    }
    if (last_dot == std::string::npos || last_dot < last_slash) {
        return path.substr(last_slash);
    }
    return path.substr(last_slash, last_dot - last_slash);
}

std::string path_parent_dir(const std::string& path) {
    size_t last_slash = path.find_last_of("/\\");
    if (last_slash == std::string::npos) {
        return std::string(".");
    }
    return path.substr(0, last_slash);
}

std::string index_output_dir(const std::string& out_path) {
    if (out_path.empty()) {
        return std::string();
    }
    if (is_file_path(out_path)) {
        // Place artifacts beside the output file; "" for a bare filename so the
        // caller falls back to its default location (next to the reference).
        size_t last_slash = out_path.find_last_of("/\\");
        return (last_slash == std::string::npos) ? std::string()
                                                 : out_path.substr(0, last_slash);
    }
    // A directory output: use it directly.
    return out_path;
}

std::string sidecar_dir(const std::string& out_path) {
    if (out_path.empty()) {
        return std::string(".");
    }
    if (is_file_path(out_path)) {
        // Parent of the output file, or "." for a bare filename — never the
        // file path itself (which would create a directory named after it).
        size_t last_slash = out_path.find_last_of("/\\");
        return (last_slash == std::string::npos) ? std::string(".")
                                                 : out_path.substr(0, last_slash);
    }
    return out_path;
}

} // namespace eastr
