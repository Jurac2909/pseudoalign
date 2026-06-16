#!/bin/bash
set -e

clang++ -std=c++17 -O2 -Isrc \
    src/fasta.cpp src/dbg.cpp src/pseudoalign.cpp src/em.cpp \
    main.cpp \
    -o pseudoalign

echo "Built: ./pseudoalign"
echo "Usage:"
echo "  ./pseudoalign build <transcripts.fa> <index> [k]"
echo "  ./pseudoalign align  <index> <reads.fastq>"
