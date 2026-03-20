#include "ol_worker.h"
#include <algorithm>
#include <iostream>


OLDatabase::OLDatabase(const std::string& prefix_starts_path,
                       const std::string& suffixes_path,
                       const std::string& counts_path,
                       int k,
                       bool canonical,
                       uint32_t linear_cutoff)
    : k_(k), canonical_(canonical), linear_cutoff_(linear_cutoff) {
    
    if (k_ != 21) {
        throw std::runtime_error("OLDatabase currently only supports k=21");
    }
    
    init_lut();
    if (canonical_) {
        init_rc_tables();
    }
    
    // Read prefix starts
    std::ifstream in_starts(prefix_starts_path, std::ios::binary);
    if (!in_starts) {
        throw std::runtime_error("Cannot open prefix_starts file: " + prefix_starts_path);
    }
    in_starts.seekg(0, std::ios::end);
    size_t starts_size = in_starts.tellg() / sizeof(uint32_t);
    constexpr uint32_t P = (1u << 26);
    if (starts_size != P + 1) {
        throw std::runtime_error("prefix_starts.bin wrong length (expected 2^26+1 uint32)");
    }
    starts_.resize(starts_size);
    in_starts.seekg(0, std::ios::beg);
    in_starts.read(reinterpret_cast<char*>(starts_.data()), starts_size * sizeof(uint32_t));
    in_starts.close();
    
    // Read suffixes
    std::ifstream in_suf(suffixes_path, std::ios::binary);
    if (!in_suf) {
        throw std::runtime_error("Cannot open suffixes file: " + suffixes_path);
    }
    in_suf.seekg(0, std::ios::end);
    size_t suf_size = in_suf.tellg() / sizeof(uint16_t);
    suffixes_.resize(suf_size);
    in_suf.seekg(0, std::ios::beg);
    in_suf.read(reinterpret_cast<char*>(suffixes_.data()), suf_size * sizeof(uint16_t));
    in_suf.close();
    
    // Read counts
    std::ifstream in_cnt(counts_path, std::ios::binary);
    if (!in_cnt) {
        throw std::runtime_error("Cannot open counts file: " + counts_path);
    }
    in_cnt.seekg(0, std::ios::end);
    size_t cnt_size = in_cnt.tellg() / sizeof(uint16_t);
    counts_.resize(cnt_size);
    in_cnt.seekg(0, std::ios::beg);
    in_cnt.read(reinterpret_cast<char*>(counts_.data()), cnt_size * sizeof(uint16_t));
    in_cnt.close();
    
    if (suffixes_.size() != counts_.size()) {
        throw std::runtime_error("Suffixes and counts size mismatch");
    }
}

void OLDatabase::init_lut() {
    for (int i = 0; i < 256; ++i) LUT_[i] = 255;
    LUT_[(unsigned)'A'] = LUT_[(unsigned)'a'] = 0;
    LUT_[(unsigned)'C'] = LUT_[(unsigned)'c'] = 1;
    LUT_[(unsigned)'G'] = LUT_[(unsigned)'g'] = 2;
    LUT_[(unsigned)'T'] = LUT_[(unsigned)'t'] = 3;
}

void OLDatabase::init_rc_tables() {
    // Initialize RC16 for 8 bases (16 bits)
    for (uint32_t x = 0; x < (1u << 16); ++x) {
        uint32_t y = 0;
        for (int i = 0; i < 8; ++i) {
            uint32_t b = (x >> (2 * i)) & 3u;
            b ^= 3u;           // complement
            y = (y << 2) | b;  // reverse
        }
        RC16_[x] = (uint16_t)y;
    }
    
    // Initialize RC10 for 5 bases (10 bits)
    for (uint32_t x = 0; x < (1u << 10); ++x) {
        uint32_t y = 0;
        for (int i = 0; i < 5; ++i) {
            uint32_t b = (x >> (2 * i)) & 3u;
            b ^= 3u;
            y = (y << 2) | b;
        }
        RC10_[x] = (uint16_t)y;
    }
}

uint64_t OLDatabase::revcomp_42_fast(uint64_t x42) const {
    // k=21 => 42 bits split: low16 (8 bases), mid16 (8 bases), high10 (5 bases)
    uint16_t lo = (uint16_t)(x42 & 0xFFFFULL);
    uint16_t mi = (uint16_t)((x42 >> 16) & 0xFFFFULL);
    uint16_t hi = (uint16_t)((x42 >> 32) & 0x03FFULL);
    
    uint64_t y = 0;
    y |= (uint64_t)RC16_[lo] << 26;  // 2*(21-8)=26
    y |= (uint64_t)RC16_[mi] << 10;  // 2*(21-16)=10
    y |= (uint64_t)RC10_[hi];        // bottom 10 bits
    return y;
}

const uint16_t* OLDatabase::lower_bound_u16(const uint16_t* a, const uint16_t* b, uint16_t x) const {
    const uint16_t* first = a;
    size_t n = (size_t)(b - a);
    while (n) {
        size_t half = n >> 1;
        const uint16_t* mid = first + half;
        if (*mid < x) {
            first = mid + 1;
            n -= half + 1;
        } else {
            n = half;
        }
    }
    return first;
}

int OLDatabase::find_suffix_adaptive(const uint16_t* suf, uint32_t n, uint16_t key) const {
    if (n == 0) return -1;
    
    if (n <= linear_cutoff_) {
        // Linear search for small buckets
        for (uint32_t i = 0; i < n; ++i) {
            if (suf[i] == key) return (int)i;
        }
        return -1;
    }
    
    // Binary search for large buckets
    const uint16_t* it = lower_bound_u16(suf, suf + n, key);
    return (it != suf + n && *it == key) ? (int)(it - suf) : -1;
}

int OLDatabase::query_kmer(uint64_t& kmer42) const {
    if (canonical_) {
        uint64_t rc = revcomp_42_fast(kmer42);
        if (rc < kmer42) kmer42 = rc;
    }
    // Extract prefix (top 26 bits) and suffix (bottom 16 bits)
    const uint32_t prefix = (uint32_t)((kmer42 >> 16) & ((1u << 26) - 1u));
    const uint16_t sufkey = (uint16_t)(kmer42 & 0xFFFFu);
    
    // Get bucket boundaries
    const uint32_t Lidx = starts_[prefix];
    const uint32_t Ridx = starts_[prefix + 1];
    
    if (Lidx >= Ridx || Ridx > (uint32_t)suffixes_.size()) {
        return 0;  // Empty bucket or invalid
    }
    
    const uint32_t len = Ridx - Lidx;
    if (len == 0) return 0;
    
    // Search for suffix in bucket
    int off = find_suffix_adaptive(suffixes_.data() + Lidx, len, sufkey);
    if (off >= 0) {
        return (int)counts_[(size_t)Lidx + (size_t)off];
    }
    
    return 0;  // Not found
}


OLWorker::OLWorker(const std::string& prefix_starts_path,
                   const std::string& suffixes_path,
                   const std::string& counts_path,
                   int k,
                   bool canonical)
    : k_(k), canonical_(canonical) {
    
    db_ = std::make_unique<OLDatabase>(prefix_starts_path, suffixes_path, 
                                       counts_path, k, canonical);
}

bool OLWorker::encode_kmer(const std::string& seq, size_t pos, uint64_t& kmer42) const {
    if (pos + k_ > seq.length()) return false;
    
    kmer42 = 0;
    for (int i = 0; i < k_; ++i) {
        uint8_t b = db_->LUT_[(unsigned char)seq[pos + i]];
        if (b == 255) return false;  // Invalid base
        kmer42 = (kmer42 << 2) | (uint64_t)b;
    }
    return true;
}

ReadKmerData OLWorker::process_sequence(const SequenceRecord& seq) {
    ReadKmerData data;
    data.read_name = seq.name;
    
    const std::string& sequence = seq.sequence;
    
    if (sequence.length() < (size_t)k_) {
        return data;  // Sequence too short
    }
    
    // Rolling k-mer encoding
    constexpr uint64_t MASK42 = (1ULL << 42) - 1ULL;
    uint64_t rolling = 0;
    uint64_t query_kmer = 0;
    int valid_len = 0;
    
    for (size_t i = 0; i < sequence.length(); ++i) {
        uint8_t b = db_->LUT_[(unsigned char)sequence[i]];
        
        if (b == 255) {
            // Invalid base - reset
            rolling = 0;
            valid_len = 0;
            continue;
        }
        
        rolling = ((rolling << 2) | (uint64_t)b) & MASK42;
        if (valid_len < k_) {
            ++valid_len;
            if (valid_len < k_) continue;
        }

        query_kmer = rolling;
        
        // Query the database
        uint64_t count = db_->query_kmer(query_kmer);
        
        // Store k-mer and count
        data.kmer_positions.push_back(query_kmer);
        data.counts.push_back(count);
    }
    
    return data;
}

std::vector<ReadKmerData> OLWorker::process_batch(const std::vector<SequenceRecord>& batch) {
    std::vector<ReadKmerData> results;
    results.reserve(batch.size());
    
    // Process each sequence in the batch
    // Can be parallelized with OpenMP if desired
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < batch.size(); ++i) {
        ReadKmerData data = process_sequence(batch[i]);
        
        #pragma omp critical
        {
            results.push_back(std::move(data));
        }
    }
    
    return results;
}
