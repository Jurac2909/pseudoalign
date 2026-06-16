# pseudokall

A de Bruijn graph pseudoalignment tool built in C++17, implementing the core ideas from kallisto.

## Build

```bash
clang++ -std=c++17 -O2 -Isrc src/fasta.cpp src/dbg.cpp tools/dump_dot.cpp -o dump_dot
clang++ -std=c++17 -O2 -Isrc src/fasta.cpp src/dbg.cpp tests/test_dbg.cpp -o test_dbg
```

## Usage

```bash
./dump_dot transcripts.fa <k> [out.dot]
dot -Tpng out.dot -o graph.png && open graph.png
```

## Source files

### `src/kmer.hpp`

All k-mer operations. A k-mer is stored as a `uint64_t` with each base packed into 2 bits (A=00, C=01, G=10, T=11), supporting k up to 31.

| Function | Description |
|----------|-------------|
| `encode_kmer(s, k, out)` | Encodes a string of length k into a 64-bit integer. Returns false if any non-ACGT character is found. |
| `decode_kmer(km, k)` | Converts a packed k-mer back to a string. |
| `reverse_complement(km, k)` | Returns the reverse complement of a k-mer. Each base is complemented by flipping both bits (A↔T, C↔G), then the order of bases is reversed. |
| `canonical(km, k)` | Returns the lexicographically smaller of a k-mer and its reverse complement. Used so that a sequence and its reverse complement are treated as the same node in the graph. |
| `roll_kmer(prev, next_base, k, out)` | Slides the k-mer window one base to the right: drops the leftmost base and appends the new base on the right. O(1) per step. |

---

### `src/fasta.hpp` / `src/fasta.cpp`

| Function | Description |
|----------|-------------|
| `parse_fasta(path, callback)` | Streams a FASTA file and calls `callback(id, sequence)` for each record. Handles multi-line sequences. The sequence ID is everything between `>` and the first whitespace on the header line. |

---

### `src/dbg.hpp` / `src/dbg.cpp`

The de Bruijn graph. Built in two passes over the reference transcriptome.

#### Data structures

**`Unitig`** — a group of k-mers that all appear in exactly the same set of transcripts. Fields:
- `id` — index into `DBG::unitigs`
- `tids` — sorted list of transcript IDs that contain this k-mer group (the equivalence class)

**`DBG`** — the graph itself. Fields:
- `k` — k-mer length
- `unitigs` — list of all unitigs / equivalence classes
- `kmer_to_unitig` — hash map from canonical k-mer to unitig ID
- `transcript_names` — transcript IDs in the order they were read from the FASTA
- `transcript_lengths` — lengths of each transcript

#### Methods

| Method | Description |
|--------|-------------|
| `build(fasta_path, k)` | Builds the graph from a reference FASTA file. **Pass 1:** slides a window of width k across every transcript and records, for each canonical k-mer, the set of transcripts that contain it. **Pass 2:** groups k-mers that share the same transcript set into a single unitig, assigns each a unique ID, and populates `kmer_to_unitig`. |
| `print_stats()` | Prints the number of transcripts, distinct k-mers, unitigs, and a histogram of how many unitigs are shared by exactly N transcripts. |
| `print_kmers()` | Prints every equivalence class with its transcript set and all k-mer sequences that belong to it, sorted alphabetically. |
| `to_dot(fasta_path, dot_path)` | Writes a Graphviz DOT file. Nodes are unitigs labelled with their transcript set and k-mer count. Edges are computed by re-walking the FASTA and recording which unitig transitions occur in each transcript. Node colour: blue = unique to one transcript, green = shared by all, orange = shared by some. |

#### How `build` works

**Pass 1** produces a temporary map `kmer_to_tids`:
```
canonical k-mer  →  sorted list of transcript IDs that contain it
```
A k-mer appearing multiple times in the same transcript is recorded only once per transcript.

**Pass 2** interns transcript sets into unitigs using `TidSetHash` (FNV-1a over the raw bytes of the transcript ID list). Every k-mer with the same transcript set is assigned the same unitig ID. The temporary map is then discarded.

---

### `tools/dump_dot.cpp`

Command-line tool that builds the graph from a FASTA file and writes a DOT file.

```
Usage: dump_dot <transcripts.fa> <k> [out.dot]
```

Calls `build`, `print_stats`, `print_kmers`, and `to_dot` in sequence.

---

### `tests/test_dbg.cpp`

Unit tests for all three modules. Accepts a path to a FASTA file as its first argument (defaults to `tests/data/isoform.fa`).

| Test | What it checks |
|------|---------------|
| `test_kmer_encode_decode` | encode → decode round-trip; non-ACGT input returns false |
| `test_kmer_reverse_complement` | known RC pairs; RC(RC(x)) == x |
| `test_kmer_canonical` | canonical(x) == canonical(RC(x)); canonical ≤ its RC |
| `test_kmer_rolling` | every step of `roll_kmer` matches a fresh `encode_kmer` |
| `test_fasta_basic` | correct IDs and non-empty sequences |
| `test_fasta_multiline` | multi-line sequences are concatenated; description is stripped |
| `test_dbg_transcript_count` | correct number of transcripts and lengths |
| `test_dbg_unitig_invariants` | every unitig has a non-empty sorted tids list; all indices valid |
| `test_dbg_ec_count` | all expected equivalence classes are present |
| `test_dbg_shared_exon_a` | k-mers from the shared 5' exon map to {TX1, TX2} |
| `test_dbg_unique_exon_b1` | k-mers from TX1's unique exon map to {TX1} only |
| `test_dbg_unique_exon_b2` | k-mers from TX2's unique exon map to {TX2} only |
| `test_dbg_shared_exon_c` | k-mers from the shared 3' exon map to {TX1, TX2, TX3} |
| `test_dbg_rc_same_unitig` | a k-mer and its RC resolve to the same unitig |
| `test_dbg_short_transcript` | a transcript shorter than k contributes no k-mers |
| `test_dbg_repeated_kmer` | a k-mer repeated within one transcript is counted once |
| `test_dbg_absent_kmer` | a k-mer not in the reference returns empty |

---

### `tests/data/isoform.fa`

Synthetic test transcriptome with an isoform structure (k=7):

```
TX1 = [EXON_A 15bp] + [EXON_B1 12bp] + [EXON_C 16bp]
TX2 = [EXON_A 15bp] + [EXON_B2 12bp] + [EXON_C 16bp]
TX3 = [EXON_C 16bp]
TX4 = ACGT  (shorter than k, contributes no k-mers)
```

Expected graph structure:

```
{TX1,TX2} → {TX1} → {TX1,TX2,TX3}
          → {TX2} ↗
```
