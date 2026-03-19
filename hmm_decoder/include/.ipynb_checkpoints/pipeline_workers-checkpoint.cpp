#include "pipeline_workers.h"
//#include "debug_logger.h"
#include <chrono>

// FastaReaderWorker implementation - STREAMING VERSION
FastaReaderWorker::FastaReaderWorker(const std::string& fasta_file, 
                 ThreadSafeQueue<std::vector<SequenceRecord>>& output_queue,
                 int batch_size)
    : fasta_file_(fasta_file), output_queue_(output_queue), batch_size_(batch_size) {}

void FastaReaderWorker::run() {
    try {
        std::cout << "Starting streaming FASTA reader..." << std::endl;
        
        // Use streaming reader with callback
        FastaReader::stream_batches(fasta_file_, batch_size_, 
            [this](std::vector<SequenceRecord>&& batch) {
                total_reads_ += batch.size();
                
                // Report progress every 10 batches
                static std::atomic<int> batch_count{0};
                int count = ++batch_count;
                if (count % 10 == 0) {
                    std::cout << "  Read " << total_reads_ << " sequences..." << std::endl;
                }
                
                // Push batch to queue
                output_queue_.push(std::move(batch));
            });
        
        std::cout << "Finished reading " << total_reads_ << " sequences from FASTA file" << std::endl;
        output_queue_.finish();
        
    } catch (const std::exception& e) {
        std::cerr << "FastaReader error: " << e.what() << std::endl;
        output_queue_.finish();
    }
}

// OLWorkerThread implementation
OLWorkerThread::OLWorkerThread(
    ThreadSafeQueue<std::vector<SequenceRecord>>& input_queue,
    ThreadSafeQueue<std::shared_ptr<ReadKmerData>>& output_queue,
    const std::string& prefix_starts_path,
    const std::string& suffixes_path,
    const std::string& counts_path,
    int worker_id)
    : input_queue_(input_queue), output_queue_(output_queue),
      worker_(prefix_starts_path, suffixes_path, counts_path), worker_id_(worker_id) {
    // Database is loaded once in the constructor of OLWorker
    std::cout << "OL worker " << worker_id_ << " initialized and database loaded" << std::endl;
}

void OLWorkerThread::run() {
    try {
        std::vector<SequenceRecord> batch;
        int batches_processed = 0;

        auto start_time = std::chrono::high_resolution_clock::now();
        uint64_t total_kmers = 0;
        // uint64_t local_kmers = 0;
        
        while (input_queue_.pop(batch)) {
            
            // auto local_start_time = std::chrono::high_resolution_clock::now();
            // local_kmers = 0;
            
            auto kmer_data = worker_.process_batch(batch);
            for (auto& data : kmer_data) {
                
                // local_kmers += data.counts.size();
                total_kmers += data.counts.size();

                /*
                LOG_STREAM("OL counts for read: " << data.read_name);
                //for (size_t ii = 0; ii < data.counts.size(); ++ii) {
                    //LOG_STREAM(ii << "\t" << data.kmer_positions[ii] << "\t" << data.counts[ii]);
                }
                */
                
                output_queue_.push(std::make_shared<ReadKmerData>(std::move(data)));
            }
            batches_processed++;

            /*
            //Benchmark kmers/sec
            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed_sec = std::chrono::duration<double>(current_time - local_start_time).count();
            
            long kmers_per_sec = static_cast<long>(local_kmers / elapsed_sec);
            
            std::cout << "  [Worker " << worker_id_ << "] "
                     << kmers_per_sec << " kmers/sec" << std::endl;
            */
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_elapsed = std::chrono::duration<double>(end_time - start_time).count();
        long final_kmers_per_sec = static_cast<long>(total_kmers / total_elapsed);
        
        std::cout << "[Worker " << worker_id_ << " DONE] "
                 << final_kmers_per_sec << " kmers/sec average" << std::endl;
                 
    } catch (const std::exception& e) {
        std::cerr << "OLWorker " << worker_id_ << " error: " << e.what() << std::endl;
    }
}

// HMMDecoderWorker implementation
HMMDecoderWorker::HMMDecoderWorker(
    ThreadSafeQueue<std::shared_ptr<ReadKmerData>>& input_queue,
    ThreadSafeQueue<std::pair<std::string, std::vector<HetLocation>>>& output_queue,
    const HMMParams& params,
    int worker_id)
    : input_queue_(input_queue), output_queue_(output_queue), 
      params_(params), worker_id_(worker_id) {}

void HMMDecoderWorker::run() {
    try {
        std::shared_ptr<ReadKmerData> kmer_data_ptr;
        int reads_processed = 0;
        
        while (input_queue_.pop(kmer_data_ptr)) {
            auto& kmer_data = *kmer_data_ptr;
            try {
                // Run HMM decoding
                std::vector<int> states;
                std::vector<double> probs;
                Matrix gamma;
                std::vector<int> repeat_mask;
                
                log_viterbi(params_, kmer_data.counts, states, probs, 
                           gamma, repeat_mask);
                
                kmer_data.states = states;
                kmer_data.probs = probs;
                kmer_data.repeat_mask = repeat_mask;
                
                // Extract heterozygous locations
                auto het_locs = extract_het_locs(states);
                
                // Compute window probabilities
                for (auto& loc : het_locs) {
                    std::vector<int> indices;
                    for (int idx = loc.kmer_start; idx <= loc.kmer_end && idx < static_cast<int>(states.size()); ++idx) {
                        indices.push_back(idx);
                    }

                    for (int idx : loc.het_kmer_indices) indices.push_back(idx);
                    
                    if (!indices.empty()) {
                        loc.window_prob = compute_window_prob(gamma, indices);
                        
                        // Check if in repeat region
                        int repeat_count = 0;
                        for (int idx = loc.kmer_start; idx <= loc.kmer_end && idx < static_cast<int>(repeat_mask.size()); ++idx) {
                            repeat_count += repeat_mask[idx];
                        }
                        loc.in_repeat = (repeat_count > 0.5 * (loc.kmer_end - loc.kmer_start + 1)) ? 1 : 0;
                        
                        // Get probabilities for het k-mers
                        for (int idx : loc.het_kmer_indices) {
                            loc.het_kmers.push_back(kmer_data.kmer_positions[idx]);
                            if (idx < static_cast<int>(probs.size())) {
                                loc.het_probs.push_back(std::exp(probs[idx]));
                            }
                        }
                    }
                }
                
                // Output result
                output_queue_.push({std::move(kmer_data.read_name), std::move(het_locs)});
                reads_processed++;
                
            } catch (const std::exception& e) {
                std::cerr << "Error processing read " << kmer_data.read_name 
                         << ": " << e.what() << std::endl;
            }

            /*
            if (reads_processed % 5000 == 0) {
                std::cout << "  HMM worker " << worker_id_ 
                         << " processed " << reads_processed << " reads" << std::endl;
            }
            */
        }
        
        std::cout << "HMM worker " << worker_id_ << " finished ("
                 << reads_processed << " reads)" << std::endl;
                 
    } catch (const std::exception& e) {
        std::cerr << "\nHMMDecoder " << worker_id_ << " error: " << e.what() << std::endl;
    }
}

// Output File Writer implementation
OutputFileWriter::OutputFileWriter(
    ThreadSafeQueue<std::pair<std::string, std::vector<HetLocation>>>& input_queue,
    const std::string& filename)
    : input_queue_(input_queue), filename_(filename){}

void OutputFileWriter::run() {
    // Open file once at the start
    std::ofstream output_file(filename_, std::ios::out);
    if (!output_file.is_open()) {
        throw std::runtime_error("Cannot open output file: " + filename_);
    }
    
    std::pair<std::string, std::vector<HetLocation>> result;
    while (input_queue_.pop(result)) {
        // Write immediately - no batching needed
        write_results_to_text(result, output_file);
    }
    
    output_file.close();
}