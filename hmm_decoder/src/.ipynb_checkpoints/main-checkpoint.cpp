#include <iostream>
#include <string>
#include <thread>
#include <stdexcept>
#include <cstdlib>
#include <sys/stat.h>
#include "data_structures.h"
#include "hmm_pipeline.h"
//#include "debug_logger.h"

const std::string JELLYFISH = "/home1/dhaivatj/softwares/jellyfish/jellyfish-2.3.0/bin/jellyfish";

// Check if a file exists
bool file_exists(const std::string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

// Build OL database from Jellyfish database
bool build_ol_database(const std::string& jellyfish_db,
                       const std::string& jftool_path,
                       const std::string& prefix_starts_path,
                       const std::string& suffixes_path,
                       const std::string& counts_path) {
    
    std::cout << "\n=== Building OL Database ===" << std::endl;
    std::cout << "  Input: " << jellyfish_db << std::endl;
    std::cout << "  Output: " << prefix_starts_path << ", " << suffixes_path << ", " << counts_path << std::endl;
    
    // Check if Jellyfish database exists
    if (!file_exists(jellyfish_db)) {
        std::cerr << "Error: Jellyfish database not found: " << jellyfish_db << std::endl;
        return false;
    }
    
    // Check if JFToOl tool exists
    if (!file_exists(jftool_path)) {
        std::cerr << "Error: JFToOl converter not found: " << jftool_path << std::endl;
        std::cerr << "Please compile it with:" << std::endl;
        std::cerr << "  g++ -O3 -march=native -std=c++17 -DNDEBUG JFToOl.cpp -o JFToOl" << std::endl;
        return false;
    }
    
    // Build the command: jellyfish dump -c <db> | ./JFToOl <outputs>
    std::string command = JELLYFISH + " dump -c " + jellyfish_db + 
                         " | " + jftool_path + 
                         " " + prefix_starts_path +
                         " " + suffixes_path +
                         " " + counts_path;
    
    std::cout << "  Running: " << command << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    int result = system(command.c_str());
    auto end = std::chrono::high_resolution_clock::now();
    
    if (result != 0) {
        std::cerr << "Error: Failed to build OL database (exit code: " << result << ")" << std::endl;
        return false;
    }
    
    // Verify output files were created
    if (!file_exists(prefix_starts_path) || !file_exists(suffixes_path) || !file_exists(counts_path)) {
        std::cerr << "Error: OL database files were not created" << std::endl;
        return false;
    }
    
    double elapsed = std::chrono::duration<double>(end - start).count();
    std::cout << "  OL database built successfully in " << elapsed << " seconds" << std::endl;
    
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] 
                 << " <params_text> <fasta_file> <jellyfish_db.jf>"
                 << " [output_file] [threads] [batch_size] [debug_logger_file] [--jftool=./JFToOl] [--rebuild]" 
                 << std::endl;
        std::cerr << "\nArguments:" << std::endl;
        std::cerr << "  params_text       : HMM parameters file" << std::endl;
        std::cerr << "  fasta_file        : Input FASTA/FASTQ file with reads" << std::endl;
        std::cerr << "  jellyfish_db.jf   : Jellyfish k-mer database" << std::endl;
        std::cerr << "  output_file       : Output results file (default: hmm_results.txt)" << std::endl;
        std::cerr << "  threads           : Number of threads (default: auto-detect)" << std::endl;
        std::cerr << "  batch_size        : Batch size for processing (default: 5000)" << std::endl;
        std::cerr << "\nOptions:" << std::endl;
        std::cerr << "  --jftool=PATH     : Path to JFToOl converter (default: ./JFToOl)" << std::endl;
        std::cerr << "  --rebuild         : Force rebuild of OL database even if it exists" << std::endl;
        std::cerr << "\nExample:" << std::endl;
        std::cerr << "  " << argv[0] << " params.txt reads.fasta counts.jf results.txt 32 5000" << std::endl;
        std::cerr << "\nNote: The program will automatically:" << std::endl;
        std::cerr << "  1. Check if OL database exists (counts.jf.prefix, counts.jf.suffix, counts.jf.counts)" << std::endl;
        std::cerr << "  2. Build it if missing using: jellyfish dump -c | JFToOl" << std::endl;
        std::cerr << "  3. Run the HMM pipeline" << std::endl;
        return 1;
    }

    std::string params_file = argv[1];
    std::string fasta_file = argv[2];
    std::string jellyfish_db = argv[3];
    std::string output_file = (argc > 4 && argv[4][0] != '-') ? argv[4] : "hmm_results.txt";
    
    // Parse remaining arguments
    int threads = std::thread::hardware_concurrency();
    int batch_size = 1000;
    std::string debug_file = "";
    std::string jftool_path = "./bin/JFToOl";
    bool force_rebuild = false;
    bool debugger = false;
    
    int next_arg = 5;
    if (argc > 4 && argv[4][0] != '-') next_arg = 5;
    else next_arg = 4;
    
    for (int i = next_arg; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--jftool=", 0) == 0) {
            jftool_path = arg.substr(9);
        } else if (arg == "--rebuild") {
            force_rebuild = true;
        } else if (arg[0] != '-') {
            // Positional arguments
            if (i == next_arg && next_arg == 4) {
                output_file = arg;
            } else if ((i == next_arg && next_arg == 5) || (i == next_arg + 1 && next_arg == 4)) {
                threads = std::stoi(arg);
            } else if ((i == next_arg + 1 && next_arg == 5) || (i == next_arg + 2 && next_arg == 4)) {
                batch_size = std::stoi(arg);
            } else if ((i == next_arg + 2 && next_arg == 5) || (i == next_arg + 3 && next_arg == 4)) {
                debug_file = arg;
            }
        }
    }

    /*
    if (debug_file != "") debugger = true;

    if (debugger) {
        LOG_INIT(debug_file);
        LOG_STREAM("Program started with " << argc << " arguments");
        LOG_VAR(params_file);
        LOG_VAR(fasta_file);
        LOG_VAR(jellyfish_db);
        LOG_VAR(output_file);
        LOG_VAR(threads);
        LOG_VAR(batch_size);
        LOG_VAR(debug_file);
        LOG_VAR(debugger);
    }
    */

    // Derive OL database paths from Jellyfish database path
    std::string prefix_starts_path = jellyfish_db + ".prefix";
    std::string suffixes_path = jellyfish_db + ".suffix";
    std::string counts_path = jellyfish_db + ".counts";

    try {
        // Check if OL database exists
        bool ol_exists = file_exists(prefix_starts_path) && 
                        file_exists(suffixes_path) && 
                        file_exists(counts_path);
        
        if (!ol_exists || force_rebuild) {
            if (force_rebuild && ol_exists) {
                std::cout << "Forcing rebuild of OL database..." << std::endl;
            } else {
                std::cout << "OL database not found. Building from Jellyfish database..." << std::endl;
            }
            
            if (!build_ol_database(jellyfish_db, jftool_path, 
                                  prefix_starts_path, suffixes_path, counts_path)) {
                return 1;
            }
        } else {
            std::cout << "Using existing OL database:" << std::endl;
            std::cout << "  " << prefix_starts_path << std::endl;
            std::cout << "  " << suffixes_path << std::endl;
            std::cout << "  " << counts_path << std::endl;
        }
        
        // Load HMM parameters
        std::cout << "\n=== Loading HMM Parameters ===" << std::endl;
        std::cout << "  File: " << params_file << std::endl;
        HMMParams params;
        if (!params.load_from_text(params_file)) {
            throw std::runtime_error("Failed to load HMM parameters");
        }
        params.print();

        // Create and run pipeline
        std::cout << "\n=== Running HMM Pipeline ===" << std::endl;
        std::cout << "  FASTA: " << fasta_file << std::endl;
        std::cout << "  Threads: " << threads << std::endl;
        std::cout << "  Batch size: " << batch_size << std::endl;
        
        HMMPipeline pipeline(params, fasta_file, 
                           prefix_starts_path, suffixes_path, counts_path, output_file,
                           threads, batch_size);
        //LOG_DEBUG("Running pipeline...");
        pipeline.run();
        //LOG_DEBUG("Pipeline completed successfully");

        std::cout << "\n=== Results written to file ===" << std::endl;
        std::cout << "  Output: " << output_file << std::endl;
        
        std::cout << "\n✓ All done!" << std::endl;

    } catch (const std::exception& e) {
        //LOG_STREAM("EXCEPTION: " << e.what());
        std::cerr << "\n❌ Error: " << e.what() << std::endl;
        //LOG_SHUTDOWN();
        return 1;
    }

    //LOG_SHUTDOWN();
    return 0;
}