#include "fasta.hpp"
#include <fstream>
#include <stdexcept>

void parse_fasta(const std::string& path,
                 std::function<void(const std::string&, const std::string&)> callback)
{
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open FASTA file: " + path);

    std::string line, id, seq;

    auto flush = [&]() {
        if (!id.empty() && !seq.empty())
            callback(id, seq);
    };

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush();
            id.clear(); seq.clear();
            size_t end = line.find_first_of(" \t", 1);
            id = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
        } else {
            seq += line;
        }
    }
    flush();
}
