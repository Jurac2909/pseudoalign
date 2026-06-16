#pragma once
#include <cstdint>
#include <string>

using Kmer = uint64_t;

inline uint8_t base_to_bits(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 255;
    }
}

inline char bits_to_base(uint8_t b) {
    switch (b & 3) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        default: return 'T';
    }
}

inline uint8_t complement_bits(uint8_t b) {
    switch (b) {
        case 0: return 3;
        case 1: return 2;
        case 2: return 1;
        default: return 0;
    }
}

inline bool encode_kmer(const char* s, int k, Kmer& out) {
    out = 0;
    for (int i = 0; i < k; ++i) {
        uint8_t b = base_to_bits(s[i]);
        if (b == 255) return false;
        out = (out << 2) | b;
    }
    return true;
}

inline std::string decode_kmer(Kmer km, int k) {
    std::string s(k, 'A');
    for (int i = k - 1; i >= 0; --i) {
        s[i] = bits_to_base(km & 3);
        km >>= 2;
    }
    return s;
}

inline Kmer reverse_complement(Kmer km, int k) {
    Kmer rc = 0;
    for (int i = 0; i < k; ++i) {
        rc = (rc << 2) | complement_bits(km & 3);
        km >>= 2;
    }
    return rc;
}

inline Kmer canonical(Kmer km, int k) {
    Kmer rc = reverse_complement(km, k);
    return km < rc ? km : rc;
}

inline bool roll_kmer(Kmer prev, char next_base, int k, Kmer& out) {
    uint8_t b = base_to_bits(next_base);
    if (b == 255) return false;
    Kmer mask = (k < 32) ? ((Kmer(1) << (2 * k)) - 1) : ~Kmer(0);
    out = ((prev << 2) | b) & mask;
    return true;
}
