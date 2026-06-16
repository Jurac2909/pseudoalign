#pragma once
#include "kmer.hpp"
#include <vector>
#include <cstdint>

struct KmerEntry {
    uint32_t ec_id;
    uint32_t contig_id;
    uint32_t pos;
    uint32_t contig_len;
};

// Open-addressing hash map with linear probing.
// Uses UINT64_MAX as the empty-slot sentinel (canonical k-mers for k<=31
// fit in 62 bits, so this value is never a valid key).
struct KmerMap {
    static constexpr Kmer EMPTY = ~Kmer(0);

    struct Slot {
        Kmer key = EMPTY;
        KmerEntry value = {};
    };

    std::vector<Slot> table;
    size_t count = 0;

    KmerMap() = default;

    explicit KmerMap(size_t expected) {
        size_t cap = 1;
        while (cap < expected * 10 / 7) cap <<= 1;
        table.assign(cap, Slot{});
    }

    size_t size() const { return count; }

    KmerEntry* find(Kmer key) {
        if (table.empty()) return nullptr;
        size_t mask = table.size() - 1;
        size_t idx = mix(key) & mask;
        while (table[idx].key != EMPTY) {
            if (table[idx].key == key) return &table[idx].value;
            idx = (idx + 1) & mask;
        }
        return nullptr;
    }

    const KmerEntry* find(Kmer key) const {
        if (table.empty()) return nullptr;
        size_t mask = table.size() - 1;
        size_t idx = mix(key) & mask;
        while (table[idx].key != EMPTY) {
            if (table[idx].key == key) return &table[idx].value;
            idx = (idx + 1) & mask;
        }
        return nullptr;
    }

    void insert(Kmer key, const KmerEntry& val) {
        size_t mask = table.size() - 1;
        size_t idx = mix(key) & mask;
        while (table[idx].key != EMPTY && table[idx].key != key)
            idx = (idx + 1) & mask;
        if (table[idx].key == EMPTY) ++count;
        table[idx] = {key, val};
    }

    Slot* begin() { return table.data(); }
    Slot* end()   { return table.data() + table.size(); }
    const Slot* begin() const { return table.data(); }
    const Slot* end()   const { return table.data() + table.size(); }

private:
    static size_t mix(Kmer k) {
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccdULL;
        k ^= k >> 33;
        return static_cast<size_t>(k);
    }
};

// Lighter map used only during index construction: Kmer -> uint32_t (ec_id).
// 12 bytes per slot vs 24 for KmerMap.
struct KmerMap32 {
    static constexpr Kmer EMPTY = ~Kmer(0);

    struct Slot {
        Kmer     key   = EMPTY;
        uint32_t value = 0;
    };

    std::vector<Slot> table;
    size_t count = 0;

    KmerMap32() = default;

    explicit KmerMap32(size_t expected) {
        size_t cap = 1;
        while (cap < expected * 10 / 7) cap <<= 1;
        table.assign(cap, Slot{});
    }

    uint32_t* find(Kmer key) {
        if (table.empty()) return nullptr;
        size_t mask = table.size() - 1;
        size_t idx  = mix(key) & mask;
        while (table[idx].key != EMPTY) {
            if (table[idx].key == key) return &table[idx].value;
            idx = (idx + 1) & mask;
        }
        return nullptr;
    }

    // insert may grow the table; any pointer returned by find() before
    // this call must NOT be used afterward.
    void insert(Kmer key, uint32_t val) {
        if ((count + 1) * 10 > table.size() * 7) grow();
        size_t mask = table.size() - 1;
        size_t idx  = mix(key) & mask;
        while (table[idx].key != EMPTY && table[idx].key != key)
            idx = (idx + 1) & mask;
        if (table[idx].key == EMPTY) ++count;
        table[idx] = {key, val};
    }

    void grow() {
        std::vector<Slot> bigger(table.size() * 2, Slot{});
        size_t mask = bigger.size() - 1;
        for (const auto& s : table) {
            if (s.key == EMPTY) continue;
            size_t idx = mix(s.key) & mask;
            while (bigger[idx].key != EMPTY) idx = (idx + 1) & mask;
            bigger[idx] = s;
        }
        table = std::move(bigger);
    }

private:
    static size_t mix(Kmer k) {
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccdULL;
        k ^= k >> 33;
        return static_cast<size_t>(k);
    }
};
