#include "dbg.hpp"
#include "fasta.hpp"
#include "fastq.hpp"
#include "pseudoalign.hpp"
#include "em.hpp"
#include <atomic>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_set>
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

// Return the idx-th tab-separated field of a SAM record.
static std::string sam_field(const std::string& s, int idx)
{
    int f = 0; size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == '\t') {
            if (f == idx) return s.substr(start, i - start);
            ++f; start = i + 1;
        }
    }
    return "";
}

// EM-aware pseudobam: rewrite the SAM keeping, for each read, only alignments to
// transcripts that survive the EM ("active"). This mirrors kallisto, whose
// --pseudobam reports reads only against expressed transcripts rather than the
// full raw equivalence class. Records arrive grouped by read (written together),
// so we filter one read's group at a time -- O(1) memory. If every alignment of
// a read is filtered out, it is re-emitted as unmapped.
static void filter_sam_emaware(const std::string& path,
                               const std::unordered_set<std::string>& active)
{
    std::ifstream in(path);
    std::string tmp = path + ".tmp";
    std::ofstream out(tmp);

    std::string line, cur_q;
    std::vector<std::string> group;

    auto flush = [&]() {
        if (group.empty()) return;
        if (sam_field(group[0], 2) == "*") {   // unmapped read, pass through
            out << group[0] << "\n";
            return;
        }
        std::vector<const std::string*> keep;
        for (const auto& rec : group)
            if (active.count(sam_field(rec, 2))) keep.push_back(&rec);
        if (keep.empty()) {                     // all assignments filtered -> unmapped
            out << cur_q << "\t4\t*\t0\t0\t*\t*\t0\t0\t"
                << sam_field(group[0], 9) << "\t*\n";
            return;
        }
        for (size_t j = 0; j < keep.size(); ++j) {
            const std::string& r = *keep[j];
            int flag = (j == 0) ? 0 : 256;      // first survivor primary, rest secondary
            out << sam_field(r, 0) << "\t" << flag << "\t" << sam_field(r, 2)
                << "\t" << sam_field(r, 3) << "\t" << sam_field(r, 4) << "\t"
                << sam_field(r, 5) << "\t*\t0\t0\t" << sam_field(r, 9) << "\t*\n";
        }
    };

    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '@') { out << line << "\n"; continue; }
        size_t tab = line.find('\t');
        std::string q = (tab == std::string::npos) ? line : line.substr(0, tab);
        if (q != cur_q) { flush(); group.clear(); cur_q = q; }
        group.push_back(line);
    }
    flush();
    in.close(); out.close();
    std::rename(tmp.c_str(), path.c_str());
}

static void cmd_align(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: pseudoalign align <index> <reads.fastq> [-l mean_frag] [-s sd_frag] [--pseudobam <out.sam>] [--min-count <c>] [--sam-raw]\n";
        std::exit(1);
    }

    std::string index_path = argv[2];
    std::string fastq      = argv[3];
    double mean_fl = 200.0, sd_fl = 30.0;
    double min_count = 0.0;   // post-EM prune threshold (0 = disabled)
    bool sam_raw = false;     // --sam-raw: emit full EC instead of EM-aware SAM
    std::string sam_path;

    for (int i = 4; i < argc; ++i) {
        std::string flag = argv[i];
        if ((flag == "-l" || flag == "-s") && i + 1 < argc) {
            if (flag == "-l") mean_fl = std::stod(argv[++i]);
            else              sd_fl   = std::stod(argv[++i]);
        } else if (flag == "--pseudobam" && i + 1 < argc) {
            sam_path = argv[++i];
        } else if (flag == "--min-count" && i + 1 < argc) {
            min_count = std::stod(argv[++i]);
        } else if (flag == "--sam-raw") {
            sam_raw = true;
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

    // Optional post-EM prune: zero transcripts whose estimated count is below
    // --min-count and renormalise TPM over the survivors. Off by default
    // (min_count == 0). NOTE: kallisto itself keeps fractional counts, so this
    // does NOT make abundances closer to kallisto -- it cuts low-count
    // over-detections at the cost of some true low-expressers (lower Pearson r).
    // Provided for downstream analyses that want a cleaner transcript list.
    if (min_count > 0.0) {
        uint64_t pruned = 0;
        double tpm_kept = 0.0;
        for (size_t t = 0; t < em.counts.size(); ++t) {
            if (em.counts[t] > 0.0 && em.counts[t] < min_count) {
                em.counts[t] = 0.0;
                em.tpm[t]    = 0.0;
                ++pruned;
            } else {
                tpm_kept += em.tpm[t];
            }
        }
        if (tpm_kept > 0.0)
            for (double& x : em.tpm) x = x / tpm_kept * 1e6;  // renormalise to sum 1e6
        std::cout << "Post-EM prune    : min-count=" << min_count
                  << "  -> pruned " << pruned << " transcripts\n\n";
    }

    std::string tsv_path = "results.tsv";
    std::ofstream tsv(tsv_path);
    tsv << "transcript\test_counts\ttpm\n";
    for (size_t t = 0; t < g.transcript_names.size(); t++) {
        tsv << g.transcript_names[t] << "\t"
            << std::fixed << std::setprecision(1) << em.counts[t] << "\t"
            << std::setprecision(0) << em.tpm[t] << "\n";
    }
    std::cout << "TSV written to: " << tsv_path << "\n";

    // EM-aware pseudobam: restrict each read's SAM alignments to expressed
    // transcripts (estimated count >= active threshold), matching kallisto.
    // --sam-raw disables this and keeps the full equivalence class.
    if (!sam_path.empty() && !sam_raw) {
        double active_min = (min_count > 0.0) ? min_count : 1.0;
        std::unordered_set<std::string> active;
        for (size_t t = 0; t < g.transcript_names.size(); ++t)
            if (em.counts[t] >= active_min) active.insert(g.transcript_names[t]);
        filter_sam_emaware(sam_path, active);
        std::cout << "EM-aware pseudobam: kept " << active.size()
                  << " expressed transcripts (count >= " << active_min << ")\n";
    }

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
