#ifndef HMM_PIPELINE_H
#define HMM_PIPELINE_H

#include <string>
#include <vector>
#include <unordered_map>
#include <thread>
#include <chrono>
#include <algorithm>
#include "data_structures.h"
#include "thread_safe_queue.h"
#include "pipeline_workers.h"

class HMMPipeline {
public:
    HMMPipeline(const HMMParams& params,
                const std::string& fasta_file,
                const std::string& prefix_starts_path,
                const std::string& suffixes_path,
                const std::string& counts_path,
                const std::string& output_filename,
                int total_threads = 32,
                int batch_size = 5000);
    
    void run();

private:
    HMMParams params_;
    std::string fasta_file_;
    std::string prefix_starts_path_;
    std::string suffixes_path_;
    std::string counts_path_;
    std::string output_filename_;
    int num_ol_workers_;
    int num_hmm_workers_;
    int batch_size_;
};

#endif