#pragma once
#include "hash_map.hpp"
#include <vector>
#include <string>
#include <cstdint>

struct Unitig {
    uint32_t              id;
    std::vector<uint32_t> tids;
};

struct DBG {
    int k = 0;
    std::vector<Unitig>      unitigs;
    KmerMap                  kmer_index;
    std::vector<std::string> transcript_names;
    std::vector<uint32_t>    transcript_lengths;

    using Transcript = std::pair<std::string, std::string>;

    void build(const std::vector<Transcript>& transcripts, int k);
    void save(const std::string& path) const;
    void load(const std::string& path);
    void print_stats() const;
};
