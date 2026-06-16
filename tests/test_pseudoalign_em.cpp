#include "dbg.hpp"
#include "pseudoalign.hpp"
#include "em.hpp"
#include "fasta.hpp"
#include <chrono>
#include <iomanip>
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

static std::string result_str(const DBG& g, const std::vector<uint32_t>& r) {
    if (r.empty()) return "(unmapped)";
    std::string s = "{ ";
    for (uint32_t t : r) s += g.transcript_names[t] + " ";
    return s + "}";
}

static void run(const DBG& g, const std::string& read,
                const std::string& label, std::vector<uint32_t> expect)
{
    auto result = pseudoalign(read, g);
    bool pass = (result == expect);
    if (!pass) ++g_failures;
    std::cout << (pass ? "  PASS" : "  FAIL")
              << "  [" << label << "]\n"
              << "        read   : " << read << "\n"
              << "        result : " << result_str(g, result) << "\n\n";
}

int main(int argc, char* argv[]) {
    std::string fa = (argc > 1) ? argv[1] : "tests/data/medium.fa";

    std::vector<DBG::Transcript> transcripts;
    parse_fasta(fa, [&](const std::string& id, const std::string& seq) {
        transcripts.push_back({id, seq});
    });

    auto t0 = std::chrono::steady_clock::now();

    DBG g;
    g.build(transcripts, 7);
    g.print_stats();
    std::cout << "\n";

    run(g, transcripts[0].second.substr(30, 25), "TRANS1 pos 30", {0});
    run(g, transcripts[1].second.substr(30, 25), "TRANS2 pos 30", {1});
    run(g, transcripts[2].second.substr(30, 25), "TRANS3 pos 30", {2});
    run(g, transcripts[0].second.substr(60, 25), "TRANS1 pos 60", {0});
    run(g, transcripts[1].second.substr(60, 25), "TRANS2 pos 60", {1});
    run(g, transcripts[2].second.substr(60, 25), "TRANS3 pos 60", {2});
    run(g, "ATTTAAGCGCATCGCGCGCGATAGA",            "TRANS3 with read error", {});

    std::cout << "--- EM quantification ---\n";

    ECCounts ec_counts;
    uint64_t total = 0, mapped = 0;
    const int read_len = 25;

    for (const auto& [name, seq] : transcripts) {
        for (size_t i = 0; i + read_len <= seq.size(); ++i) {
            std::string read = seq.substr(i, read_len);
            auto result = pseudoalign(read, g);
            ++total;
            if (!result.empty()) { ec_counts[result]++; ++mapped; }
        }
    }

    std::cout << "  Total reads   : " << total << "\n";
    std::cout << "  Mapped reads  : " << mapped << "\n";
    std::cout << "  Mapping rate  : "
              << std::fixed << std::setprecision(1)
              << 100.0 * mapped / total << "%\n\n";

    EMResult em = run_em(ec_counts, g);

    std::cout << "  Transcript abundances (TPM):\n";
    for (size_t t = 0; t < g.transcript_names.size(); ++t) {
        std::cout << "    " << std::left << std::setw(20) << g.transcript_names[t]
                  << "  est. counts: " << std::setw(8) << std::fixed << std::setprecision(1)
                  << em.counts[t]
                  << "  TPM: " << std::setprecision(0) << em.tpm[t] << "\n";
    }

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "\n=== " << (g_failures == 0 ? "ALL PASSED" : "SOME FAILED")
              << " (" << g_failures << " failure(s))  time: " << ms << " ms ===\n";
    return g_failures > 0 ? 1 : 0;
}
