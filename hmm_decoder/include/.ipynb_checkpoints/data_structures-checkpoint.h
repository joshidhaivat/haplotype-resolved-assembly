#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include "matrix.h"

struct SequenceRecord {
    std::string name;
    std::string sequence;
};

struct ReadKmerData {
    std::string read_name;
    std::vector<unsigned long int> kmer_positions;
    std::vector<unsigned long int> counts;
    std::vector<int> states;
    std::vector<double> probs;
    std::vector<int> repeat_mask;
};

struct HetLocation {
    int start = -1;
    int end = -1;
    std::vector<unsigned long int> het_kmers;
    std::vector<int> het_kmer_indices;
    std::vector<double> het_probs;
    int total_kmers = -1;
    double window_prob = 0.0;
    int in_repeat = 0;
    int kmer_start = -1;
    int kmer_end = -1;
};

struct HMMParams {
    std::vector<double> initial;
    Matrix transition;
    std::vector<double> emission_mean;
    std::vector<double> emission_var;
    size_t num_states;

    HMMParams() : num_states(0) {}

    bool load_from_text(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open parameter file: " << filename << std::endl;
            return false;
        }

        std::string line;
        
        if (!std::getline(file, line)) return false;
        num_states = std::stoull(line);
        
        initial.resize(num_states);
        for (size_t i = 0; i < num_states; ++i) {
            if (!std::getline(file, line)) return false;
            initial[i] = std::stod(line);
        }
        
        transition.resize(num_states, num_states);
        for (size_t i = 0; i < num_states; ++i) {
            for (size_t j = 0; j < num_states; ++j) {
                if (!std::getline(file, line)) return false;
                transition(i, j) = std::stod(line);
            }
        }
        
        emission_mean.resize(num_states);
        for (size_t i = 0; i < num_states; ++i) {
            if (!std::getline(file, line)) return false;
            emission_mean[i] = std::stod(line);
        }
        
        emission_var.resize(num_states);
        for (size_t i = 0; i < num_states; ++i) {
            if (!std::getline(file, line)) return false;
            emission_var[i] = std::stod(line);
        }
        
        return true;
    }

    void print() const {
        std::cout << "HMM Parameters:\n";
        std::cout << "  Number of states: " << num_states << "\n";
        std::cout << "  Initial probs (exp): [";
        for (size_t i = 0; i < num_states; ++i) {
            std::cout << std::exp(initial[i]);
            if (i < num_states - 1) std::cout << ", ";
        }
        std::cout << "]\n";
        std::cout << "  Emission means: [";
        for (size_t i = 0; i < num_states; ++i) {
            std::cout << emission_mean[i];
            if (i < num_states - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
};

#endif // DATA_STRUCTURES_H
