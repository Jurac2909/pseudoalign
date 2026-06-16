#!/bin/bash
set -e

clang++ -std=c++17 -O2 -Isrc \
    src/fasta.cpp src/dbg.cpp src/pseudoalign.cpp \
    tests/test_pseudoalign.cpp \
    -o test_pseudoalign

./test_pseudoalign "$@"
