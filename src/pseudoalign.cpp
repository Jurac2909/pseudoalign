#include "pseudoalign.hpp"
#include <algorithm>
#include <atomic>
#include <unordered_map>

std::atomic<uint64_t> g_skip_count{0};
std::atomic<uint64_t> g_kmer_count{0};

// Pseudoalignment by k-compatibility-class intersection (kallisto-style).
//
// We walk the read's k-mers, look up each one's equivalence class (the set of
// transcripts containing it), and report the transcripts consistent with the
// whole read. Instead of intersecting greedily in sequence order — which
// anchors the answer to whichever EC the first k-mer happened to hit and can
// only ever stay the same or grow — we count, for every transcript, how many of
// the read's k-mer positions fall in an EC that contains it. The
// pseudoalignment is then:
//
//   * the transcripts present in *every* hit k-mer (support == n_hit): a true
//     intersection of all k-compatibility classes, as in kallisto; or, when
//     that intersection is empty (a few k-mers land on a sequencing error or a
//     paralogous repeat),
//   * the transcripts with maximal support (plurality), so the read survives.
//
// Both cases reduce to one rule: keep the transcripts whose support equals the
// best support observed.
std::vector<uint32_t> pseudoalign(const std::string& read, const DBG& g,
                                  uint32_t min_len)
{
    if ((int)read.size() < g.k) return {};

    // ec_support[ec] = number of read k-mer positions whose EC is `ec`.
    // thread_local so the buffer is reused across calls rather than reallocated
    // per read (and ready for the planned per-thread alignment).
    thread_local std::unordered_map<uint32_t, uint32_t> ec_support;
    ec_support.clear();
    uint64_t n_hit = 0;  // read k-mer positions found in the index

    size_t i = 0;
    while (i + (size_t)g.k <= read.size()) {
        Kmer km;
        if (!encode_kmer(read.c_str() + i, g.k, km)) { ++i; continue; }
        Kmer ck = canonical(km, g.k);

        ++g_kmer_count;
        const KmerEntry* e = g.kmer_index.find(ck);
        if (!e) { ++i; continue; }

        // Contig-skip: consecutive k-mers inside one contig share an EC, so if
        // the k-mer `remaining` positions ahead is confirmed to sit in the same
        // contig, credit the whole run at once instead of looking each one up.
        uint32_t run = 1;
        uint32_t remaining = e->contig_len - e->pos - 1;
        if (remaining > 0) {
            size_t check = i + remaining;
            if (check + (size_t)g.k > read.size())
                check = read.size() - (size_t)g.k;
            if (check > i) {
                Kmer km2;
                if (encode_kmer(read.c_str() + check, g.k, km2)) {
                    const KmerEntry* e2 = g.kmer_index.find(canonical(km2, g.k));
                    if (e2 && e2->contig_id == e->contig_id) {
                        uint64_t skipped = check - i;
                        g_skip_count += skipped;
                        run += static_cast<uint32_t>(skipped);
                        i = check;  // ++i below advances past the confirmed k-mer
                    }
                }
            }
        }

        ec_support[e->ec_id] += run;
        n_hit += run;
        ++i;
    }

    if (n_hit == 0) return {};

    // Tally per-transcript support across the hit ECs.
    thread_local std::unordered_map<uint32_t, uint64_t> tx_support;
    tx_support.clear();
    for (const auto& [ec, cnt] : ec_support)
        for (uint32_t t : g.unitigs[ec].tids)
            tx_support[t] += cnt;

    uint64_t best = 0;
    for (const auto& [t, s] : tx_support) best = std::max(best, s);

    // Reject reads where even the best-supported transcript is backed by fewer
    // than half the hit k-mers — likely chimeric or dominated by errors.
    if (best * 2 < n_hit) return {};

    std::vector<uint32_t> candidates;
    for (const auto& [t, s] : tx_support)
        if (s == best) candidates.push_back(t);
    std::sort(candidates.begin(), candidates.end());

    if (min_len > 0) {
        candidates.erase(
            std::remove_if(candidates.begin(), candidates.end(),
                [&](uint32_t t){ return g.transcript_lengths[t] < min_len; }),
            candidates.end());
    }
    return candidates;
}
