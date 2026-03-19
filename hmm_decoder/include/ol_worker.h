#ifndef OL_WORKER_H
#define OL_WORKER_H

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include "data_structures.h"

// Forward declaration
class OLWorker;

// Lookup table for base encoding
class OLDatabase {
public:
    OLDatabase(const std::string& prefix_starts_path,
               const std::string& suffixes_path,
               const std::string& counts_path,
               int k = 21,
               bool canonical = true,
               uint32_t linear_cutoff = 32);
    
    // Query a single k-mer (42-bit encoded)
    int query_kmer(uint64_t& kmer42) const;
    
    // Get k-mer size
    int get_k() const { return k_; }
    
    // Allow OLWorker to access private members
    friend class OLWorker;

private:
    int k_;
    bool canonical_;
    uint32_t linear_cutoff_;
    
    // Database arrays (loaded once per worker)
    std::vector<uint32_t> starts_;    // 2^26 + 1 entries
    std::vector<uint16_t> suffixes_;  // N entries (sorted within prefix)
    std::vector<uint16_t> counts_;    // N entries
    
    // Lookup tables
    uint8_t LUT_[256];
    uint16_t RC16_[1u << 16];  // reverse-complement for 8 bases
    uint16_t RC10_[1u << 10];  // reverse-complement for 5 bases
    
    // Initialize lookup tables
    void init_lut();
    void init_rc_tables();
    
    // Fast reverse complement for 42-bit k-mer (k=21)
    uint64_t revcomp_42_fast(uint64_t x42) const;
    
    // Adaptive search (linear for small buckets, binary for large)
    int find_suffix_adaptive(const uint16_t* suf, uint32_t n, uint16_t key) const;
    
    // Binary search helper
    const uint16_t* lower_bound_u16(const uint16_t* a, const uint16_t* b, uint16_t x) const;
};

class OLWorker {
public:
    OLWorker(const std::string& prefix_starts_path,
             const std::string& suffixes_path,
             const std::string& counts_path,
             int k = 21,
             bool canonical = true);

    // Process a batch of sequences and return k-mer counts
    std::vector<ReadKmerData> process_batch(const std::vector<SequenceRecord>& batch);

private:
    std::unique_ptr<OLDatabase> db_;
    int k_;
    bool canonical_;
    
    // Helper function to query a single sequence
    ReadKmerData process_sequence(const SequenceRecord& seq);
    
    // Encode a k-mer string to 42-bit integer
    bool encode_kmer(const std::string& seq, size_t pos, uint64_t& kmer42) const;
};

#endif // OL_WORKER_H