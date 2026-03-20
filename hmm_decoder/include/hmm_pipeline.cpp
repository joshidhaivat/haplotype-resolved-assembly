#include "hmm_pipeline.h"
#include <iostream>

HMMPipeline::HMMPipeline(const HMMParams& params,
            const std::string& fasta_file,
            const std::string& prefix_starts_path,
            const std::string& suffixes_path,
            const std::string& counts_path,
            const std::string& output_filename,
            int total_threads,
            int batch_size)
    : params_(params), fasta_file_(fasta_file), 
      prefix_starts_path_(prefix_starts_path),
      suffixes_path_(suffixes_path),
      counts_path_(counts_path),
      output_filename_(output_filename),
      batch_size_(batch_size) {
    
    // Allocate threads
    num_ol_workers_ = std::max(1, std::min(1, total_threads / 10));
    num_hmm_workers_ = std::max(1, total_threads - num_ol_workers_ - 2);
    
    std::cout << "Pipeline configuration:" << std::endl;
    std::cout << "  Total threads: " << total_threads << std::endl;
    std::cout << "  OL workers: " << num_ol_workers_ << std::endl;
    std::cout << "  HMM workers: " << num_hmm_workers_ << std::endl;
    std::cout << "  Batch size: " << batch_size_ << std::endl;
}

void HMMPipeline::run() {
    std::cout << "\nStarting pipeline..." << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    size_t queue_capacity = std::max(3, num_ol_workers_) + 2;  // e.g., 8 workers + 2 = 10 batches max
    
    ThreadSafeQueue<std::vector<SequenceRecord>> fasta_to_ol(queue_capacity);
    ThreadSafeQueue<std::shared_ptr<ReadKmerData>> ol_to_hmm(2 * queue_capacity * batch_size_);
    ThreadSafeQueue<std::pair<std::string, std::vector<HetLocation>>> hmm_to_output(0);
    
    std::cout << "  Queue capacity: " << queue_capacity << " batches (backpressure enabled)" << std::endl;
    
    // Start FASTA reader thread
    std::cout << "\nStage 1: Reading FASTA file..." << std::endl;
    FastaReaderWorker fasta_reader(fasta_file_, fasta_to_ol, batch_size_);
    std::thread fasta_thread([&fasta_reader]() { fasta_reader.run(); });
    
    // Start OL worker threads
    std::cout << "\nStage 2: Starting " << num_ol_workers_ << " OL workers..." << std::endl;
    std::vector<std::thread> ol_threads;
    for (int i = 0; i < num_ol_workers_; ++i) {
        ol_threads.emplace_back([this, &fasta_to_ol, &ol_to_hmm, i]() {
            OLWorkerThread worker(fasta_to_ol, ol_to_hmm, 
                                 prefix_starts_path_, suffixes_path_, counts_path_, i);
            worker.run();
        });
    }
    
    // Start HMM decoder threads
    std::cout << "\nStage 3: Starting " << num_hmm_workers_ << " HMM workers..." << std::endl;
    std::vector<std::thread> hmm_threads;
    for (int i = 0; i < num_hmm_workers_; ++i) {
        hmm_threads.emplace_back([this, &ol_to_hmm, &hmm_to_output, i]() {
            HMMDecoderWorker worker(ol_to_hmm, hmm_to_output, params_, i);
            worker.run();
        });
    }

    // Start Output writer threads
    std::cout << "\nStage 4: Starting 1 output writer worker..." << std::endl;
    OutputFileWriter file_writer(hmm_to_output, output_filename_);
    std::thread output_thread([&file_writer]() { file_writer.run(); });

    // Monitor queue sizes periodically
    std::atomic<bool> monitoring{false};
    std::thread monitor_thread([&]() {
        while (monitoring) {
            std::this_thread::sleep_for(std::chrono::seconds(5));
            if (!monitoring) break;
            
            size_t q1 = fasta_to_ol.size();
            size_t q2 = ol_to_hmm.size();
            size_t q3 = hmm_to_output.size();
            
            if (q1 > 0 || q2 > 0 || q3 > 0) {
                std::cout << "  Queue sizes: FASTA→OL=" << q1 
                         << ", OL→HMM=" << q2 
                         << ", HMM→OUT=" << q3 << std::endl;
            }
        }
    });
    
    // Wait for FASTA reader to finish
    fasta_thread.join();
    std::cout << "\nFASTA reading complete" << std::endl;
    
    // Wait for OL workers to finish, then signal HMM workers
    for (auto& t : ol_threads) {
        t.join();
    }
    ol_to_hmm.finish();
    std::cout << "\nOL processing complete" << std::endl;
    
    // Wait for HMM workers to finish
    for (auto& t : hmm_threads) {
        t.join();
    }
    hmm_to_output.finish();
    std::cout << "\nHMM decoding complete" << std::endl;

    /*
    // Collect results
    std::cout << "\nCollecting results..." << std::endl;
    std::unordered_map<std::string, std::vector<HetLocation>> results;
    std::pair<std::string, std::vector<HetLocation>> result;
    
    while (hmm_to_output.pop(result)) {
        results[result.first] = result.second;
    }
    */

    // Wait for Output writer to finish
    output_thread.join();
    std::cout << "\nOutput writing complete" << std::endl;

    // Stop monitoring
    monitoring = false;
    monitor_thread.join();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\n✓ Pipeline complete!" << std::endl;
    std::cout << "  Total time: " << duration.count() << " seconds" << std::endl;
}