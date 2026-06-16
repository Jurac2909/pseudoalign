// Diagnostic: for each k-mer in a read, print its EC from the index
#include "dbg.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: dbg_kmer <index> <read_seq>\n";
        return 1;
    }

    DBG g;
    g.load(argv[1]);

    std::string read = argv[2];
    int k = g.k;

    std::cout << "k=" << k << "  transcripts=" << g.transcript_names.size() << "\n\n";

    for (int i = 0; i + k <= (int)read.size(); ++i) {
        Kmer km;
        if (!encode_kmer(read.c_str() + i, k, km)) {
            std::cout << "pos " << i << ": invalid k-mer\n";
            continue;
        }
        Kmer ck = canonical(km, k);
        const KmerEntry* e = g.kmer_index.find(ck);
        if (!e) {
            std::cout << "pos " << i << ": not in index\n";
        } else {
            const auto& tids = g.unitigs[e->ec_id].tids;
            std::cout << "pos " << i << ": contig=" << e->contig_id
                      << "  ec_id=" << e->ec_id
                      << "  EC_size=" << tids.size() << " [";
            for (size_t j = 0; j < std::min((size_t)5, tids.size()); ++j) {
                if (j) std::cout << ",";
                std::cout << g.transcript_names[tids[j]];
            }
            if (tids.size() > 5) std::cout << ",...";
            std::cout << "]\n";
        }
    }
    return 0;
}
