#include "output_writer.h"

void write_results_to_text(
    const std::unordered_map<std::string, std::vector<HetLocation>>& results,
    const std::string& filename) {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open output file: " + filename);
    }
    
    for (const auto& [read_name, het_locs] : results) {
        file << read_name;
        
        std::ostringstream locs_ss, kmers_ss, probs_ss, total_ss, window_ss, repeat_ss;
        
        for (size_t i = 0; i < het_locs.size(); ++i) {
            if (i > 0) {
                locs_ss << ";";
                kmers_ss << ";";
                probs_ss << ";";
                total_ss << ";";
                window_ss << ";";
                repeat_ss << ";";
            }
            
            locs_ss << het_locs[i].start << "," << het_locs[i].end;
            
            for (size_t j = 0; j < het_locs[i].het_kmers.size(); ++j) {
                if (j > 0) kmers_ss << ",";
                kmers_ss << het_locs[i].het_kmers[j];
            }
            
            for (size_t j = 0; j < het_locs[i].het_probs.size(); ++j) {
                if (j > 0) probs_ss << ",";
                probs_ss << std::fixed << std::setprecision(6) << het_locs[i].het_probs[j];
            }
            
            total_ss << het_locs[i].total_kmers;
            window_ss << std::fixed << std::setprecision(6) << het_locs[i].window_prob;
            repeat_ss << het_locs[i].in_repeat;
        }
        
        file << "\t" << locs_ss.str()
             << "\t" << kmers_ss.str()
             << "\t" << probs_ss.str()
             << "\t" << total_ss.str()
             << "\t" << window_ss.str()
             << "\t" << repeat_ss.str()
             << "\n";
    }
}
