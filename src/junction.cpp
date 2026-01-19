#include "eastr/junction.hpp"

#include <sstream>

namespace eastr {

std::string JunctionKey::to_string() const {
    std::ostringstream oss;
    oss << chrom << ":" << start << "-" << end << "(" << to_char(strand) << ")";
    return oss.str();
}

bool parse_region_string(const std::string& region, std::string& chrom,
                         int64_t& start, int64_t& end) {
    // Parse format: "chr:start-end"
    auto colon_pos = region.find(':');
    if (colon_pos == std::string::npos) {
        return false;
    }

    auto dash_pos = region.find('-', colon_pos);
    if (dash_pos == std::string::npos) {
        return false;
    }

    chrom = region.substr(0, colon_pos);

    try {
        start = std::stoll(region.substr(colon_pos + 1, dash_pos - colon_pos - 1));
        end = std::stoll(region.substr(dash_pos + 1));
    } catch (const std::exception&) {
        return false;
    }

    return true;
}

} // namespace eastr
