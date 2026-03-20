#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <string>
#include <vector>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <thread>
#include <atomic>
#include "data_structures.h"

class FastaReader {
public:
    static void stream_batches_parallel(const std::string& filename, 
                                       int batch_size,
                                       std::function<void(std::vector<SequenceRecord>&&)> callback) {
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open FASTA file: " + filename);
        }

        // Use larger internal buffer for faster I/O
        constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
        std::vector<char> buffer(BUFFER_SIZE);
        file.rdbuf()->pubsetbuf(buffer.data(), BUFFER_SIZE);

        std::string line;
        std::vector<SequenceRecord> batch;
        batch.reserve(batch_size);
        
        SequenceRecord current;
        bool in_sequence = false;
        size_t total_reads = 0;
        int batches_sent = 0;

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            if (line[0] == '>' || line[0] == '@') {
                
                if (!current.sequence.empty()) {
                    batch.push_back(std::move(current));
                    total_reads++;
                    
                    // Send batch when full
                    if (batch.size() >= (size_t)batch_size) {
                        callback(std::move(batch));
                        batches_sent++;
                        batch.clear();
                        batch.reserve(batch_size);
                    }
                }
                
                current.name = line.substr(1);
                current.sequence.clear();
                current.sequence.reserve(20000);
                in_sequence = true;
            } else if (in_sequence) {
                current.sequence += line;
                in_sequence = false;
            }
        }

        if (!current.sequence.empty()) {
            batch.push_back(std::move(current));
            total_reads++;
        }
        
        if (!batch.empty()) {
            callback(std::move(batch));
            batches_sent++;
        }
    }

    static void stream_batches(const std::string& filename, 
                               int batch_size,
                               std::function<void(std::vector<SequenceRecord>&&)> callback) {
        stream_batches_parallel(filename, batch_size, callback);
    }

    static std::vector<SequenceRecord> read_all(const std::string& filename) {
        std::vector<SequenceRecord> records;
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open FASTA file: " + filename);
        }

        std::string line;
        SequenceRecord current;
        bool in_sequence = false;

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            if (line[0] == '>' || line[0] == '@') {
                if (!current.sequence.empty()) {
                    records.push_back(current);
                }
                
                current.name = line.substr(1);
                current.sequence.clear();
                in_sequence = true;
            } else if (in_sequence) {
                current.sequence += line;
                in_sequence = false;
            }
        }

        if (!current.sequence.empty()) {
            records.push_back(current);
        }

        return records;
    }
};

#endif