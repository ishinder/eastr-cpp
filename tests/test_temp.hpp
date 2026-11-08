#pragma once

#include <atomic>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unistd.h>

namespace eastr_test {

// A unique scratch directory per fixture, including when CTest runs in parallel.
class TempDir {
public:
    TempDir() {
        for (int attempt = 0; attempt < 100; ++attempt) {
            path_ = std::filesystem::temp_directory_path() /
                    ("eastr_test_" + std::to_string(::getpid()) + "_" +
                     std::to_string(next_id_++));
            if (std::filesystem::create_directory(path_)) return;
        }
        throw std::runtime_error("Cannot create a unique test directory");
    }

    ~TempDir() {
        std::error_code ignored;
        std::filesystem::remove_all(path_, ignored);
    }

    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;

    std::string str() const { return path_.string(); }
    std::string file(const std::string& name) const { return (path_ / name).string(); }

private:
    inline static std::atomic<unsigned> next_id_{0};
    std::filesystem::path path_;
};

} // namespace eastr_test
