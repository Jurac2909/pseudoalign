#pragma once
#include "dbg.hpp"
#include <vector>
#include <string>
#include <atomic>

extern std::atomic<uint64_t> g_skip_count;
extern std::atomic<uint64_t> g_kmer_count;

// min_len: transcripts shorter than this are removed from the result.
// Set to read.size() to match kallisto 0.50+'s behaviour (a read can't
// come from a transcript shorter than itself).
std::vector<uint32_t> pseudoalign(const std::string& read, const DBG& g,
                                  uint32_t min_len = 0);
