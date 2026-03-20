#ifndef PIPELINE_WORKERS_H
#define PIPELINE_WORKERS_H

#include <string>
#include <vector>
#include <memory>
#include <atomic>
#include <iostream>
#include "thread_safe_queue.h"
#include "data_structures.h"
#include "fasta_reader.h"
#include "ol_worker.h"
#include "hmm_functions.h"
#include "output_writer.h"

class FastaReaderWorker {
public:
    FastaReaderWorker(const std::string& fasta_file, 
                     ThreadSafeQueue<std::vector<SequenceRecord>>& output_queue,
                     int batch_size = 5000);
    void run();
    size_t get_total_reads() const { return total_reads_; }

private:
    std::string fasta_file_;
    ThreadSafeQueue<std::vector<SequenceRecord>>& output_queue_;
    int batch_size_;
    std::atomic<size_t> total_reads_{0};
};

class OLWorkerThread {
public:
    OLWorkerThread(ThreadSafeQueue<std::vector<SequenceRecord>>& input_queue,
                   ThreadSafeQueue<std::shared_ptr<ReadKmerData>>& output_queue,
                   const std::string& prefix_starts_path,
                   const std::string& suffixes_path,
                   const std::string& counts_path,
                   int worker_id);
    void run();

private:
    ThreadSafeQueue<std::vector<SequenceRecord>>& input_queue_;
    ThreadSafeQueue<std::shared_ptr<ReadKmerData>>& output_queue_;
    OLWorker worker_;
    int worker_id_;
};

class HMMDecoderWorker {
public:
    HMMDecoderWorker(ThreadSafeQueue<std::shared_ptr<ReadKmerData>>& input_queue,
                    ThreadSafeQueue<std::pair<std::string, std::vector<HetLocation>>>& output_queue,
                    const HMMParams& params,
                    int worker_id);
    void run();

private:
    ThreadSafeQueue<std::shared_ptr<ReadKmerData>>& input_queue_;
    ThreadSafeQueue<std::pair<std::string, std::vector<HetLocation>>>& output_queue_;
    const HMMParams& params_;
    int worker_id_;
};

class OutputFileWriter {
public:
    OutputFileWriter(ThreadSafeQueue<std::pair<std::string, std::vector<HetLocation>>>& input_queue,
                    const std::string& filename);
    
    void run();

private:
    ThreadSafeQueue<std::pair<std::string, std::vector<HetLocation>>>& input_queue_;
    const std::string& filename_;
};

#endif