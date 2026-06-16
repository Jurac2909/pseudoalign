#include "dbg.hpp"
#include "pseudoalign.hpp"
#include "fasta.hpp"
#include <iostream>
#include <string>
#include <vector>

static int failures = 0;

static void check(const DBG& g, const std::string& read,
                  const std::string& label, std::vector<uint32_t> expected)
{
    std::vector<uint32_t> result = pseudoalign(read, g);

    bool pass = (result == expected);
    if (!pass) failures++;

    std::cout << (pass ? "PASS" : "FAIL") << "  " << label << "\n";
    std::cout << "  read:   " << read << "\n";
    std::cout << "  result: ";
    if (result.empty()) {
        std::cout << "(unmapped)";
    } else {
        for (int i = 0; i < result.size(); i++) {
            if (i > 0) std::cout << ", ";
            std::cout << g.transcript_names[result[i]];
        }
    }
    std::cout << "\n\n";
}

int main(int argc, char* argv[]) {
    std::string fa = (argc > 1) ? argv[1] : "tests/data/exon.fa";

    std::vector<DBG::Transcript> transcripts;
    parse_fasta(fa, [&](const std::string& id, const std::string& seq) {
        transcripts.push_back({id, seq});
    });

    DBG g;
    g.build(transcripts, 7);
    g.print_stats();
    std::cout << "\n";

    check(g, "ATGATGATGATGATGCGATCGA", "TX1 unique region",       {0});
    check(g, "ATGATGATGATGATGGCTAGCT", "TX2 unique region",       {1});
    check(g, "TTACGTTACGTTACGTTACGTT", "shared TX1+TX2+TX3",      {0, 1, 2});
    check(g, "ATGATGATGATGATGATGATGA", "shared TX1+TX2",          {0, 1});
    check(g, "TTACGTTACGTTACGT",       "EXON_C only",             {0, 1, 2});
    check(g, "AAACCCGGGTTTAAACCCGGGT", "not in reference",        {});

    std::cout << (failures == 0 ? "ALL PASSED" : "SOME FAILED")
              << " (" << failures << " failure(s))\n";
    return failures > 0 ? 1 : 0;
}
