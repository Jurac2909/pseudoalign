#include "dbg.hpp"
#include "fasta.hpp"
#include "fastq.hpp"
#include "pseudoalign.hpp"
#include "em.hpp"
#include <atomic>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

static void cmd_build(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: pseudoalign build <transcripts.fa> <index> [k]\n";
        std::exit(1);
    }

    std::string fa         = argv[2];
    std::string index_path = argv[3];
    int k                  = (argc >= 5) ? std::stoi(argv[4]) : 31;

    auto t0 = std::chrono::steady_clock::now();

    std::vector<DBG::Transcript> transcripts;
    parse_fasta(fa, [&](const std::string& id, const std::string& seq) {
        transcripts.push_back({id, seq});
    });
    std::cout << "Loaded " << transcripts.size() << " transcripts\n";

    DBG g;
    g.build(transcripts, k);
    g.print_stats();

    g.save(index_path);
    std::cout << "Index written to: " << index_path << "\n";

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "Done in " << ms << " ms\n";
}

static void cmd_align(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: pseudoalign align <index> <reads.fastq> [-l mean_frag] [-s sd_frag] [--pseudobam <out.sam>]\n";
        std::exit(1);
    }

    std::string index_path = argv[2];
    std::string fastq      = argv[3];
    double mean_fl = 200.0, sd_fl = 30.0;
    std::string sam_path;

    for (int i = 4; i < argc; ++i) {
        std::string flag = argv[i];
        if ((flag == "-l" || flag == "-s") && i + 1 < argc) {
            if (flag == "-l") mean_fl = std::stod(argv[++i]);
            else              sd_fl   = std::stod(argv[++i]);
        } else if (flag == "--pseudobam" && i + 1 < argc) {
            sam_path = argv[++i];
        }
    }

    auto t0 = std::chrono::steady_clock::now();

    DBG g;
    g.load(index_path);
    std::cout << "Loaded index: k=" << g.k
              << "  transcripts=" << g.transcript_names.size()
              << "  k-mers=" << g.kmer_index.size()
              << "  ECs=" << g.unitigs.size() << "\n\n";

    std::ofstream sam_out;
    if (!sam_path.empty()) {
        sam_out.open(sam_path);
        if (!sam_out.is_open()) {
            std::cerr << "Cannot open SAM output: " << sam_path << "\n";
            std::exit(1);
        }
        // SAM header
        sam_out << "@HD\tVN:1.6\tSO:unsorted\n";
        for (size_t t = 0; t < g.transcript_names.size(); ++t)
            sam_out << "@SQ\tSN:" << g.transcript_names[t]
                    << "\tLN:" << g.transcript_lengths[t] << "\n";
        sam_out << "@PG\tID:pseudokall\tPN:pseudokall\n";
    }

    ECCounts ec_counts;
    uint64_t total = 0, mapped = 0;

    parse_fastq(fastq, [&](const std::string& name, const std::string& seq) {
        ++total;
        std::vector<uint32_t> result = pseudoalign(seq, g,
                                           static_cast<uint32_t>(seq.size()));

        if (!result.empty()) {
            ec_counts[result]++;
            ++mapped;
        }
        if (sam_out.is_open()) {
            if (result.empty()) {
                // unmapped
                sam_out << name << "\t4\t*\t0\t0\t*\t*\t0\t0\t" << seq << "\t*\n";
            } else {
                // primary alignment to first transcript in EC; secondary records for the rest
                for (size_t j = 0; j < result.size(); ++j) {
                    uint32_t t = result[j];
                    int flag = (j == 0) ? 0 : 256; // 256 = secondary
                    sam_out << name << "\t" << flag << "\t"
                            << g.transcript_names[t] << "\t1\t0\t"
                            << seq.size() << "M\t*\t0\t0\t"
                            << seq << "\t*\n";
                }
            }
        }
    });

    // Close SAM before printing to stdout so a SIGPIPE from a downstream pipe
    // (e.g. | head) doesn't kill the process before the file buffer is flushed.
    if (sam_out.is_open()) sam_out.close();

    std::cout << "Reads processed : " << total << "\n";
    std::cout << "Mapped          : " << mapped << "\n";
    std::cout << "Mapping rate    : "
              << std::fixed << std::setprecision(1)
              << 100.0 * mapped / total << "%\n";
    {
        uint64_t looked = g_kmer_count.load();
        uint64_t skipped = g_skip_count.load();
        uint64_t total_pos = looked + skipped;
        std::cout << "K-mers skipped  : " << skipped << " / " << total_pos
                  << "  (" << std::setprecision(1)
                  << 100.0 * skipped / total_pos << "%)\n";
    }
    std::cout << "\n";

    std::cout << "Fragment length  : mean=" << mean_fl << "  sd=" << sd_fl << "\n\n";
    EMResult em = run_em(ec_counts, g, mean_fl, sd_fl);

    std::string tsv_path = "results.tsv";
    std::ofstream tsv(tsv_path);
    tsv << "transcript\test_counts\ttpm\n";
    for (size_t t = 0; t < g.transcript_names.size(); t++) {
        tsv << g.transcript_names[t] << "\t"
            << std::fixed << std::setprecision(1) << em.counts[t] << "\t"
            << std::setprecision(0) << em.tpm[t] << "\n";
    }
    std::cout << "TSV written to: " << tsv_path << "\n";

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "Done in " << ms << " ms\n";
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage:\n"
                  << "  pseudoalign build <transcripts.fa> <index> [k]\n"
                  << "  pseudoalign align  <index> <reads.fastq> [-l mean_frag] [-s sd_frag]\n";
        return 1;
    }

    std::string cmd = argv[1];
    if (cmd == "build") cmd_build(argc, argv);
    else if (cmd == "align") cmd_align(argc, argv);
    else {
        std::cerr << "Unknown command: " << cmd << "\n";
        return 1;
    }

    return 0;
}
