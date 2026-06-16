#pragma once
#include "dbg.hpp"
#include <vector>
#include <unordered_map>
#include <cstdint>

struct VecHash {
    size_t operator()(const std::vector<uint32_t>& v) const {
        size_t h = 0;
        for (uint32_t x : v) {
            h ^= x + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

using ECCounts = std::unordered_map<std::vector<uint32_t>, uint64_t, VecHash>;

struct EMResult {
    std::vector<double> counts;  // estimated read count per transcript
    std::vector<double> tpm;     // transcripts per million
};

EMResult run_em(const ECCounts& ec_counts, const DBG& g,
                double mean_fl = 200.0, double sd_fl = 30.0,
                int max_iter = 500, double tol = 1e-4);
