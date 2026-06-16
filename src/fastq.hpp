#pragma once
#include <fstream>
#include <functional>
#include <string>
#include <stdexcept>

void parse_fastq(const std::string& path,
                 std::function<void(const std::string& name, const std::string& seq)> callback)
{
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open FASTQ file: " + path);

    std::string line;
    std::string name;
    int lnum = 0;
    while (std::getline(f, line)) {
        ++lnum;
        int pos = lnum % 4;
        if (pos == 1) {
            // strip '@' and take only up to first whitespace
            size_t end = line.find_first_of(" \t", 1);
            name = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
        } else if (pos == 2) {
            callback(name, line);
        }
    }
}
