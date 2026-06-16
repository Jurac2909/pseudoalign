#pragma once
#include <string>
#include <functional>

void parse_fasta(const std::string& path,
                 std::function<void(const std::string& id,
                                    const std::string& seq)> callback);
