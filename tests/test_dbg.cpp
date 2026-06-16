#include "kmer.hpp"
#include "fasta.hpp"
#include "dbg.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static int g_failures = 0;

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::cerr << "  FAIL [" << __FILE__ << ":" << __LINE__ << "]  " #cond "\n"; \
            ++g_failures; \
        } \
    } while (0)

#define CHECK_EQ(a, b) \
    do { \
        if ((a) != (b)) { \
            std::cerr << "  FAIL [" << __FILE__ << ":" << __LINE__ << "]  " \
                      << #a << " == " #b "  (got " << (a) << " vs " << (b) << ")\n"; \
            ++g_failures; \
        } \
    } while (0)

static void section(const char* name) {
    std::cout << "\n--- " << name << " ---\n";
}

static std::vector<uint32_t> lookup(const DBG& g, const std::string& s) {
    if (static_cast<int>(s.size()) != g.k) return {};
    Kmer km;
    if (!encode_kmer(s.c_str(), g.k, km)) return {};
    Kmer ck = canonical(km, g.k);
    auto it = g.kmer_index.find(ck);
    if (it == g.kmer_index.end()) return {};
    return g.unitigs[it->second.ec_id].tids;
}

static bool tids_eq(const std::vector<uint32_t>& got,
                    std::initializer_list<uint32_t> want_init)
{
    std::vector<uint32_t> want(want_init);
    return got == want;
}

static void test_kmer_encode_decode() {
    section("kmer: encode / decode round-trip");
    for (const char* s : {"ACGT", "TTTT", "AAAA", "CGCG", "TGCA", "ATCGATC"}) {
        Kmer km;
        int len = static_cast<int>(strlen(s));
        CHECK(encode_kmer(s, len, km));
        CHECK_EQ(decode_kmer(km, len), std::string(s));
    }
    { Kmer km; CHECK(!encode_kmer("ACGN", 4, km)); }
    { Kmer km; CHECK(!encode_kmer("N", 1, km));    }
}

static void test_kmer_reverse_complement() {
    section("kmer: reverse complement");
    { Kmer km; encode_kmer("ACGT", 4, km); CHECK_EQ(decode_kmer(reverse_complement(km, 4), 4), "ACGT"); }
    { Kmer km; encode_kmer("AAAC", 4, km); CHECK_EQ(decode_kmer(reverse_complement(km, 4), 4), "GTTT"); }
    { Kmer km; encode_kmer("GTTT", 4, km); CHECK_EQ(decode_kmer(reverse_complement(km, 4), 4), "AAAC"); }
    { Kmer km; encode_kmer("ATGATGA", 7, km); CHECK_EQ(reverse_complement(reverse_complement(km, 7), 7), km); }
    { Kmer km; encode_kmer("ATGATGA", 7, km); CHECK_EQ(decode_kmer(reverse_complement(km, 7), 7), "TCATCAT"); }
}

static void test_kmer_canonical() {
    section("kmer: canonical form");
    for (const char* s : {"ATGATGA", "CGATCGA", "GCTAGCT", "TTACGTT"}) {
        Kmer km, rc;
        encode_kmer(s, 7, km);
        rc = reverse_complement(km, 7);
        CHECK_EQ(canonical(km, 7), canonical(rc, 7));
    }
    for (const char* s : {"ATGATGA", "CGATCGA", "GCTAGCT"}) {
        Kmer km; encode_kmer(s, 7, km);
        Kmer c = canonical(km, 7);
        Kmer rc = reverse_complement(c, 7);
        CHECK(c <= rc);
    }
}

static void test_kmer_rolling() {
    section("kmer: rolling hash matches fresh encode");
    const std::string seq = "ATGATGATGATGATG";
    const int k = 7;
    Kmer km;
    CHECK(encode_kmer(seq.c_str(), k, km));
    for (size_t i = 1; i + k <= seq.size(); ++i) {
        Kmer rolled, fresh;
        CHECK(roll_kmer(km, seq[i + k - 1], k, rolled));
        CHECK(encode_kmer(seq.c_str() + i, k, fresh));
        CHECK_EQ(rolled, fresh);
        km = rolled;
    }
}

static void test_fasta_basic(const std::string& fa) {
    section("fasta: basic parse");
    std::vector<std::pair<std::string,std::string>> recs;
    parse_fasta(fa, [&](const std::string& id, const std::string& seq) {
        recs.push_back({id, seq});
    });
    CHECK_EQ(recs.size(), 4u);
    CHECK_EQ(recs[0].first, "TX1");
    CHECK_EQ(recs[1].first, "TX2");
    CHECK_EQ(recs[2].first, "TX3");
    CHECK_EQ(recs[3].first, "TX4_too_short");
    for (auto& r : recs) CHECK(!r.second.empty());
}

static void test_fasta_multiline() {
    section("fasta: multi-line sequences are concatenated");
    const std::string path = "/tmp/pk_multiline.fa";
    {
        std::ofstream f(path);
        f << ">SEQ description is ignored\n"
          << "ACGT\n"
          << "ACGT\n"
          << ">SEQ2\n"
          << "TTTT\n";
    }
    std::vector<std::pair<std::string,std::string>> recs;
    parse_fasta(path, [&](const std::string& id, const std::string& seq) {
        recs.push_back({id, seq});
    });
    CHECK_EQ(recs.size(), 2u);
    CHECK_EQ(recs[0].first, "SEQ");
    CHECK_EQ(recs[0].second, "ACGTACGT");
    CHECK_EQ(recs[1].second, "TTTT");
}

static DBG build_isoform_dbg(const std::string& fa) {
    std::vector<DBG::Transcript> transcripts;
    parse_fasta(fa, [&](const std::string& id, const std::string& seq) {
        transcripts.push_back({id, seq});
    });
    DBG g;
    g.build(transcripts, 7);
    return g;
}

static void test_dbg_transcript_count(const DBG& g) {
    section("dbg: transcript metadata");
    CHECK_EQ(g.transcript_names.size(), 4u);
    CHECK_EQ(g.transcript_names[0], "TX1");
    CHECK_EQ(g.transcript_names[1], "TX2");
    CHECK_EQ(g.transcript_names[2], "TX3");
    CHECK_EQ(g.transcript_names[3], "TX4_too_short");
    CHECK(g.transcript_lengths[0] > (uint32_t)g.k);
    CHECK(g.transcript_lengths[3] < (uint32_t)g.k);
}

static void test_dbg_shared_exon_a(const DBG& g) {
    section("dbg: EXON_A k-mers map to {TX1, TX2}");
    auto tids = lookup(g, "ATGATGA");
    CHECK(!tids.empty());
    CHECK(tids_eq(tids, {0, 1}));
    auto tids2 = lookup(g, "GATGATG");
    CHECK(!tids2.empty());
    CHECK(tids_eq(tids2, {0, 1}));
}

static void test_dbg_unique_exon_b1(const DBG& g) {
    section("dbg: EXON_B1 k-mers map to {TX1} only");
    auto tids = lookup(g, "CGATCGA");
    CHECK(!tids.empty());
    CHECK(tids_eq(tids, {0}));
    auto tids2 = lookup(g, "ATCGATC");
    CHECK(!tids2.empty());
    CHECK(tids_eq(tids2, {0}));
}

static void test_dbg_unique_exon_b2(const DBG& g) {
    section("dbg: EXON_B2 k-mers map to {TX2} only");
    auto tids = lookup(g, "GCTAGCT");
    CHECK(!tids.empty());
    CHECK(tids_eq(tids, {1}));
    auto tids2 = lookup(g, "CTAGCTA");
    CHECK(!tids2.empty());
    CHECK(tids_eq(tids2, {1}));
}

static void test_dbg_shared_exon_c(const DBG& g) {
    section("dbg: EXON_C k-mers map to {TX1, TX2, TX3}");
    for (const char* s : {"TTACGTT", "TACGTTA", "ACGTTAC", "CGTTACG"}) {
        auto tids = lookup(g, s);
        CHECK(!tids.empty());
        CHECK(tids_eq(tids, {0, 1, 2}));
    }
}

static void test_dbg_rc_same_unitig(const DBG& g) {
    section("dbg: a k-mer and its RC map to the same unitig");
    struct Case { std::string fwd; std::string rc; };
    std::vector<Case> cases = {
        {"ATGATGA", "TCATCAT"},
        {"CGATCGA", "TCGATCG"},
        {"GCTAGCT", "AGCTAGC"},
        {"TTACGTT", "AACGTAA"},
    };
    for (auto& c : cases) {
        Kmer fwd_km, rc_km;
        CHECK(encode_kmer(c.fwd.c_str(), g.k, fwd_km));
        CHECK(encode_kmer(c.rc.c_str(),  g.k, rc_km));
        CHECK_EQ(decode_kmer(reverse_complement(fwd_km, g.k), g.k), c.rc);
        Kmer c1 = canonical(fwd_km, g.k);
        Kmer c2 = canonical(rc_km,  g.k);
        CHECK_EQ(c1, c2);
        auto it1 = g.kmer_index.find(c1);
        auto it2 = g.kmer_index.find(c2);
        if (it1 != g.kmer_index.end()) {
            CHECK(it2 != g.kmer_index.end());
            CHECK_EQ(it1->second.ec_id, it2->second.ec_id);
        }
    }
}

static void test_dbg_short_transcript(const DBG& g) {
    section("dbg: transcript shorter than k contributes no k-mers");
    for (const auto& u : g.unitigs)
        CHECK(std::find(u.tids.begin(), u.tids.end(), 3u) == u.tids.end());
}

static void test_dbg_repeated_kmer(const DBG& g) {
    section("dbg: a k-mer repeated within one transcript is counted once");
    auto tids = lookup(g, "TTACGTT");
    CHECK(!tids.empty());
    int count = std::count(tids.begin(), tids.end(), 2u);
    CHECK_EQ(count, 1);
}

static void test_dbg_absent_kmer(const DBG& g) {
    section("dbg: a k-mer not in the reference is not found");
    auto tids = lookup(g, "AAAAAAA");
    CHECK(tids.empty());
}

static void test_dbg_unitig_invariants(const DBG& g) {
    section("dbg: structural invariants on every unitig");
    CHECK(!g.unitigs.empty());
    CHECK(!g.kmer_index.empty());
    for (const auto& u : g.unitigs) {
        CHECK(!u.tids.empty());
        CHECK(std::is_sorted(u.tids.begin(), u.tids.end()));
        for (uint32_t t : u.tids)
            CHECK(t < g.transcript_names.size());
    }
    for (const auto& [km, e] : g.kmer_index)
        CHECK(e.ec_id < g.unitigs.size());
}

static void test_dbg_ec_count(const DBG& g) {
    section("dbg: expected equivalence classes are present");
    auto has_ec = [&](std::initializer_list<uint32_t> want) {
        std::vector<uint32_t> w(want);
        for (const auto& u : g.unitigs)
            if (u.tids == w) return true;
        return false;
    };
    CHECK(has_ec({0, 1}));
    CHECK(has_ec({0}));
    CHECK(has_ec({1}));
    CHECK(has_ec({0, 1, 2}));
}

int main(int argc, char* argv[]) {
    std::string fa = (argc > 1) ? argv[1] : "tests/data/isoform.fa";

    test_kmer_encode_decode();
    test_kmer_reverse_complement();
    test_kmer_canonical();
    test_kmer_rolling();

    test_fasta_basic(fa);
    test_fasta_multiline();

    DBG g = build_isoform_dbg(fa);
    g.print_stats();

    test_dbg_transcript_count(g);
    test_dbg_unitig_invariants(g);
    test_dbg_ec_count(g);
    test_dbg_shared_exon_a(g);
    test_dbg_unique_exon_b1(g);
    test_dbg_unique_exon_b2(g);
    test_dbg_shared_exon_c(g);
    test_dbg_rc_same_unitig(g);
    test_dbg_short_transcript(g);
    test_dbg_repeated_kmer(g);
    test_dbg_absent_kmer(g);

    std::cout << "\n=== " << (g_failures == 0 ? "ALL PASSED" : "SOME FAILED")
              << " (" << g_failures << " failure(s)) ===\n";
    return g_failures > 0 ? 1 : 0;
}
