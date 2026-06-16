#include "dbg.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

// Strip trailing poly-A if the run is longer than min_run bases.
// Matches kallisto's preprocessing: clips >10 consecutive 3'-terminal A's.
static std::string clip_polya(const std::string& seq, size_t min_run = 10)
{
    size_t end = seq.size();
    while (end > 0 && (seq[end - 1] == 'A' || seq[end - 1] == 'a'))
        --end;
    if (seq.size() - end > min_run)
        return seq.substr(0, end);
    return seq;
}

struct TidSetHash {
    size_t operator()(const std::vector<uint32_t>& v) const {
        size_t h = 0;
        for (uint32_t x : v) {
            h ^= x + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

// Kallisto-style streaming build: process one transcript at a time,
// updating each k-mer's EC in place. Never stores all (kmer, tid) pairs
// simultaneously, so memory scales with distinct k-mers + distinct EC sets
// rather than total k-mer positions.
//
// Key property: transcripts are processed in order t=0,1,...,N-1, so each
// new tid is always the largest seen for a k-mer. Extending an EC = append.
void DBG::build(const std::vector<Transcript>& transcripts, int ksize)
{
    if (ksize < 1 || ksize > 31) throw std::invalid_argument("k must be 1-31");
    k = ksize;
    uint32_t total_t = static_cast<uint32_t>(transcripts.size());

    // Estimate number of distinct k-mers to pre-size the table.
    // Starts small and grows automatically if the estimate is low.
    size_t total_bases = 0;
    for (const auto& tr : transcripts) total_bases += tr.second.size();
    KmerMap32 kmer_to_ecid(total_bases / 20);

    std::vector<std::vector<uint32_t>> ec_tids;
    std::unordered_map<std::vector<uint32_t>, uint32_t, TidSetHash> set_to_uid;

    auto get_or_create_ec = [&](std::vector<uint32_t> tids) -> uint32_t {
        auto it = set_to_uid.find(tids);
        if (it != set_to_uid.end()) return it->second;
        uint32_t id = static_cast<uint32_t>(ec_tids.size());
        ec_tids.push_back(tids);
        set_to_uid.emplace(std::move(tids), id);
        return id;
    };

    for (uint32_t tid = 0; tid < total_t; tid++) {
        const std::string& name = transcripts[tid].first;
        const std::string seq   = clip_polya(transcripts[tid].second);
        transcript_names.push_back(name);
        transcript_lengths.push_back(static_cast<uint32_t>(seq.size()));

        if ((int)seq.size() < k) continue;

        // Distinct k-mers for this transcript (per-transcript sort is cheap:
        // ~3000 elements fit in L2 cache).
        std::vector<Kmer> local;
        local.reserve(seq.size() - k + 1);
        Kmer km;
        for (size_t i = 0; i + (size_t)k <= seq.size(); i++) {
            if (i == 0) { if (!encode_kmer(seq.c_str(), k, km)) continue; }
            else        { if (!roll_kmer(km, seq[i+k-1], k, km)) continue; }
            local.push_back(canonical(km, k));
        }
        std::sort(local.begin(), local.end());
        local.erase(std::unique(local.begin(), local.end()), local.end());

        for (Kmer ck : local) {
            uint32_t* ec_ptr = kmer_to_ecid.find(ck);

            if (!ec_ptr) {
                // New k-mer: assign singleton EC.
                // insert() may grow the table, but ec_ptr is null so no alias issue.
                uint32_t ec_id = get_or_create_ec({tid});
                kmer_to_ecid.insert(ck, ec_id);
            } else {
                // Known k-mer: extend its EC with this transcript.
                // We copy the old tid list and append tid (always the largest).
                // get_or_create_ec may push_back to ec_tids but does not touch
                // kmer_to_ecid, so ec_ptr stays valid.
                std::vector<uint32_t> new_tids = ec_tids[*ec_ptr];
                new_tids.push_back(tid);
                *ec_ptr = get_or_create_ec(std::move(new_tids));
            }
        }
    }
    std::cerr << "\rDone — " << kmer_to_ecid.count << " k-mers, "
              << ec_tids.size() << " ECs\n";

    // Build unitigs from ec_tids, then free the intermediates.
    unitigs.resize(ec_tids.size());
    for (size_t i = 0; i < ec_tids.size(); i++)
        unitigs[i] = {static_cast<uint32_t>(i), ec_tids[i]};
    { decltype(ec_tids) tmp; ec_tids.swap(tmp); }
    { decltype(set_to_uid) tmp; set_to_uid.swap(tmp); }

    // Copy kmer→ec pairs out, free kmer_to_ecid, then build the final KmerMap.
    std::vector<std::pair<Kmer, uint32_t>> pairs;
    pairs.reserve(kmer_to_ecid.count);
    for (const auto& slot : kmer_to_ecid.table)
        if (slot.key != KmerMap32::EMPTY)
            pairs.push_back({slot.key, slot.value});
    { KmerMap32 tmp; kmer_to_ecid = std::move(tmp); }

    kmer_index = KmerMap(pairs.size());
    for (const auto& [ck, ec] : pairs)
        kmer_index.insert(ck, {ec, 0, 0, 1});

    // Second pass: assign contig IDs, positions, and lengths.
    //
    // A contig is a maximal run of consecutive k-mers in one transcript that
    // all share the same EC. During pseudoalignment the skip optimisation uses
    // this: if k-mer A is at contig position p of length n, and the k-mer
    // n-1-p positions later in the read is confirmed to be in the same contig,
    // every intermediate k-mer is guaranteed to have the same EC and can be
    // skipped without individual lookup.
    //
    // A k-mer is assigned on first encounter; later transcripts that contain
    // the same k-mer see contig_id != 0 and do not re-assign it.
    {
        uint32_t next_cid = 1;
        for (uint32_t tid = 0; tid < total_t; tid++) {
            const std::string seq = clip_polya(transcripts[tid].second);
            if ((int)seq.size() < k) continue;

            // Ordered list of KmerEntry pointers for this transcript.
            std::vector<KmerEntry*> ev;
            ev.reserve(seq.size() - k + 1);
            Kmer km;
            for (size_t i = 0; i + (size_t)k <= seq.size(); i++) {
                bool ok = (i == 0) ? encode_kmer(seq.c_str(), k, km)
                                   : roll_kmer(km, seq[i + k - 1], k, km);
                ev.push_back(ok ? kmer_index.find(canonical(km, k)) : nullptr);
            }

            // Find maximal runs of unassigned same-EC k-mers.
            for (size_t j = 0; j < ev.size(); ) {
                KmerEntry* e0 = ev[j];
                if (!e0 || e0->contig_id != 0) { ++j; continue; }

                uint32_t ec = e0->ec_id;
                size_t end = j;
                while (end + 1 < ev.size() &&
                       ev[end + 1] &&
                       ev[end + 1]->ec_id == ec &&
                       ev[end + 1]->contig_id == 0)
                    ++end;

                uint32_t len = static_cast<uint32_t>(end - j + 1);
                if (len > 1) {
                    uint32_t cid = next_cid++;
                    for (size_t p = j; p <= end; p++) {
                        ev[p]->contig_id  = cid;
                        ev[p]->pos        = static_cast<uint32_t>(p - j);
                        ev[p]->contig_len = len;
                    }
                }
                j = end + 1;
            }
        }
    }
}

void DBG::save(const std::string& path) const
{
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open index file for writing: " + path);

    uint32_t magic = 0x58494B50;
    f.write(reinterpret_cast<const char*>(&magic), 4);
    int32_t kk = k;
    f.write(reinterpret_cast<const char*>(&kk), 4);

    uint32_t nt = static_cast<uint32_t>(transcript_names.size());
    f.write(reinterpret_cast<const char*>(&nt), 4);
    for (uint32_t i = 0; i < nt; i++) {
        uint32_t len = static_cast<uint32_t>(transcript_names[i].size());
        f.write(reinterpret_cast<const char*>(&len), 4);
        f.write(transcript_names[i].data(), len);
        f.write(reinterpret_cast<const char*>(&transcript_lengths[i]), 4);
    }

    uint32_t nec = static_cast<uint32_t>(unitigs.size());
    f.write(reinterpret_cast<const char*>(&nec), 4);
    for (const auto& u : unitigs) {
        uint32_t ntids = static_cast<uint32_t>(u.tids.size());
        f.write(reinterpret_cast<const char*>(&ntids), 4);
        f.write(reinterpret_cast<const char*>(u.tids.data()), ntids * 4);
    }

    uint64_t cap = kmer_index.table.size();
    f.write(reinterpret_cast<const char*>(&cap), 8);
    f.write(reinterpret_cast<const char*>(kmer_index.table.data()),
            cap * sizeof(KmerMap::Slot));
}

void DBG::load(const std::string& path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open index file: " + path);

    uint32_t magic;
    f.read(reinterpret_cast<char*>(&magic), 4);
    if (magic != 0x58494B50) throw std::runtime_error("bad index file magic");

    int32_t kk;
    f.read(reinterpret_cast<char*>(&kk), 4);
    k = kk;

    uint32_t nt;
    f.read(reinterpret_cast<char*>(&nt), 4);
    transcript_names.resize(nt);
    transcript_lengths.resize(nt);
    for (uint32_t i = 0; i < nt; i++) {
        uint32_t len;
        f.read(reinterpret_cast<char*>(&len), 4);
        transcript_names[i].resize(len);
        f.read(transcript_names[i].data(), len);
        f.read(reinterpret_cast<char*>(&transcript_lengths[i]), 4);
    }

    uint32_t nec;
    f.read(reinterpret_cast<char*>(&nec), 4);
    unitigs.resize(nec);
    for (uint32_t i = 0; i < nec; i++) {
        unitigs[i].id = i;
        uint32_t ntids;
        f.read(reinterpret_cast<char*>(&ntids), 4);
        unitigs[i].tids.resize(ntids);
        f.read(reinterpret_cast<char*>(unitigs[i].tids.data()), ntids * 4);
    }

    uint64_t cap;
    f.read(reinterpret_cast<char*>(&cap), 8);
    kmer_index.table.resize(cap);
    f.read(reinterpret_cast<char*>(kmer_index.table.data()),
           cap * sizeof(KmerMap::Slot));

    kmer_index.count = 0;
    for (const auto& s : kmer_index.table)
        if (s.key != KmerMap::EMPTY) ++kmer_index.count;
}

void DBG::print_stats() const
{
    std::cout << "k                  : " << k << "\n";
    std::cout << "Transcripts        : " << transcript_names.size() << "\n";
    std::cout << "Distinct k-mers    : " << kmer_index.size() << "\n";
    std::cout << "Equivalence classes: " << unitigs.size() << "\n";
}
