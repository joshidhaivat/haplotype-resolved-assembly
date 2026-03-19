hmm_decoder/
├── include/
│   ├── kmer_utils.hpp           # K-mer encoding/decoding utilities
│   ├── thread_safe_queue.hpp    # Thread-safe queue template
│   ├── matrix.hpp                # Matrix operations for HMM
│   ├── hmm_math.hpp              # HMM mathematical functions
│   ├── hmm_params.hpp            # HMM parameters structure
│   ├── jellyfish_database.hpp   # K-mer database (sorted array)
│   ├── fasta_reader.hpp          # FASTA file parsing
│   ├── hmm_decoder.hpp           # HMM Viterbi decoding
│   ├── pipeline.hpp              # Pipeline orchestration
│   └── output_writer.hpp         # Result output formatting
│
├── src/
│   └── main.cpp                  # Main entry point
│
├── CMakeLists.txt                # CMake build configuration
├── Makefile                      # Simple makefile (alternative)
├── README.md                     # Documentation
└── examples/
    ├── convert_params.py         # Python parameter converter
    └── run_example.sh            # Example usage script
