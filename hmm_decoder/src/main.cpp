#include <iostream>
#include <string>
#include <thread>
#include <stdexcept>
#include <cstdlib>
#include <sys/stat.h>
#include "data_structures.h"
#include "hmm_pipeline.h"
//#include "debug_logger.h"

// Check if a file exists
bool file_exists(const std::string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

// Build OL database from Jellyfish database
bool build_ol_database(const std::string& jellyfish_db,
                       const std::string& jftool_path,
                       const std::string& jellyfish_path,
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
    std::string command = jellyfish_path + " dump -c " + jellyfish_db + 
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

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " -p <params> -r <reads> -d <jellyfish_db.jf> [options]\n"
              << "\nRequired arguments:\n"
              << "  -p, --params <file>       HMM parameters file\n"
              << "  -r, --reads <file>        Input FASTA/FASTQ file with reads\n"
              << "  -d, --database <file>     Jellyfish k-mer database (.jf)\n"
              << "\nOptional arguments:\n"
              << "  -o, --out <file>          Output results file (default: hmm_results.txt)\n"
              << "  -t, --threads <int>       Number of threads (default: auto-detect)\n"
              << "  -b, --batch <int>         Batch size for processing (default: 1000)\n"
              << "  -j, --jellyfish <path>    Path to jellyfish binary (default: jellyfish)\n"
              << "  -J, --jftool <path>       Path to JFToOl converter (default: ./bin/JFToOl)\n"
              << "  --rebuild                 Force rebuild of OL database even if it exists\n"
              << "  -h, --help                Show this help message\n"
              << "\nExamples:\n"
              << "  " << prog << " -p params.txt -r reads.fasta -d counts.jf\n"
              << "  " << prog << " -p params.txt -r reads.fasta -d counts.jf -o results.txt -t 32 -b 5000\n"
              << "  " << prog << " -p params.txt -r reads.fasta -d counts.jf -j /opt/jellyfish/bin/jellyfish --rebuild\n"
              << "\nNote: The program will automatically:\n"
              << "  1. Check if OL database exists (<db>.prefix, <db>.suffix, <db>.counts)\n"
              << "  2. Build it if missing using: jellyfish dump -c | JFToOl\n"
              << "  3. Run the HMM pipeline\n";
}

int main(int argc, char* argv[]) {
    // Defaults
    std::string params_file;
    std::string fasta_file;
    std::string jellyfish_db;
    std::string output_file = "hmm_results.txt";
    int threads = std::thread::hardware_concurrency();
    int batch_size = 1000;
    std::string jellyfish_path = "jellyfish";
    std::string jftool_path = "./bin/JFToOl";
    bool force_rebuild = false;

    // Show help if no arguments
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    // Parse arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // Helper: consume the next token as the value for a flag
        auto next_val = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << flag << " requires an argument\n";
                std::exit(1);
            }
            return argv[++i];
        };

        if (arg == "-p" || arg == "--params") {
            params_file = next_val(arg);
        } else if (arg == "-r" || arg == "--reads") {
            fasta_file = next_val(arg);
        } else if (arg == "-d" || arg == "--database") {
            jellyfish_db = next_val(arg);
        } else if (arg == "-o" || arg == "--out") {
            output_file = next_val(arg);
        } else if (arg == "-t" || arg == "--threads") {
            threads = std::stoi(next_val(arg));
        } else if (arg == "-b" || arg == "--batch") {
            batch_size = std::stoi(next_val(arg));
        } else if (arg == "-j" || arg == "--jellyfish") {
            jellyfish_path = next_val(arg);
        } else if (arg == "-J" || arg == "--jftool") {
            jftool_path = next_val(arg);
        } else if (arg == "--rebuild") {
            force_rebuild = true;
        } else if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Error: Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    // Validate required arguments
    if (params_file.empty() || fasta_file.empty() || jellyfish_db.empty()) {
        std::cerr << "Error: --params, --reads, and --database are required.\n\n";
        print_usage(argv[0]);
        return 1;
    }

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
            
            if (!build_ol_database(jellyfish_db, jftool_path, jellyfish_path,
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