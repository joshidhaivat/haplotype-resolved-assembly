# hmm_decoder

A multi-threaded C++ pipeline for detecting heterozygous genomic sites from nanopore sequencing reads using Hidden Markov Model (HMM) decoding on k-mer count profiles.

## Overview

`hmm_decoder` identifies heterozygous regions within long reads (e.g., Oxford Nanopore) by analyzing per-position k-mer count profiles against a whole-genome k-mer database. The core idea is that k-mers spanning heterozygous sites will appear at roughly half the homozygous coverage depth, and an HMM can segment each read into homozygous vs. heterozygous states based on these count fluctuations.

### Algorithm

The pipeline operates in four streaming stages connected by thread-safe queues:

**Stage 1 — FASTA Reading.** Reads are streamed from a FASTA/FASTQ file in configurable batches (default 1000).

**Stage 2 — K-mer Count Lookup (OL Workers).** Each read's sequence is scanned with a sliding window of k=21. Every 21-mer is encoded as a 42-bit integer, and looked up in a memory-mapped Ordered List (OL) database. The OL database splits each 21-mer into a 13-bit prefix (2^26 entries after 2-bit encoding) and an 8-base suffix (16-bit), enabling O(1) prefix lookup followed by a short binary or linear search within the suffix bucket. This produces a vector of k-mer counts along the read.

**Stage 3 — HMM Decoding (HMM Workers).** The count vector is decoded by a 5-state Gaussian HMM operating in log-space:

| State | Interpretation | Typical emission mean |
|-------|---------------|----------------------|
| 0 | Error | ~0 |
| 1 | Heterozygous | ~½ × homozygous mean |
| 2 | Homozygous | sample-specific coverage |
| 3 | Repeat (moderate) | > homozygous mean |
| 4 | Repeat (high) | >> homozygous mean |

Emission probabilities use a discretized Gaussian model (CDF-difference for integer counts). Transition probabilities are modulated by the magnitude of count-gradient between adjacent positions: large jumps widen the off-diagonal transition mass, while flat regions lock the identity matrix — implemented via a cached `compute_logAt_fast`.

Decoding proceeds in **two Viterbi passes**:

1. **Pass 1** produces an initial state path used to estimate local coverage via a sliding window (w=1000). Positions dominated by states 3/4 are flagged as repeat regions.
2. **Pass 2** re-runs Viterbi with locally adjusted emission means (the homozygous mean is shifted to the windowed average of state-1/2 counts), improving accuracy in regions with coverage fluctuation.

After the second Viterbi pass, forward–backward inference computes posterior state probabilities (gamma matrix). Heterozygous segments are extracted by a sliding-window scan (k=21, threshold fraction f=0.5) that merges overlapping windows and reports per-segment confidence.

**Stage 4 — Output Writing.** Results are written per-read as a tab-separated line (see [Output Format](#output-format)).

### Pipeline Architecture

```
┌──────────┐    ┌───────────┐    ┌──────────────┐    ┌──────────────┐
│  FASTA   │───>│ OL Worker │───>│  HMM Decoder │───>│ Output Writer│
│  Reader  │    │  (1+ thr) │    │  (N threads) │    │  (1 thread)  │
│ (1 thr)  │    │           │    │              │    │              │
└──────────┘    └───────────┘    └──────────────┘    └──────────────┘
           Queue            Queue               Queue
        (bounded)        (bounded)           (unbounded)
```

Thread allocation is automatic: given T total threads, 1 is assigned to FASTA reading, 1 to OL lookup, 1 to output writing, and the remainder (T-3) to HMM decoding, which is the compute bottleneck. Bounded queues provide backpressure to prevent memory blow-up.

## Dependencies

- **C++17** compiler (GCC >= 7 or Clang >= 5)
- **OpenMP** (optional, for parallelized emission probability computation)
- **pthreads**
- **zlib** (`-lz`)
- **Jellyfish 2.x** — used only at the OL database build step (`jellyfish dump -c`). Not needed at runtime if the OL database files already exist.

## Building

```bash
git clone <repo-url>
cd hmm_decoder
make
```

This produces two binaries in `bin/`:

- `hmm_decoder` — the main pipeline executable
- `JFToOl` — a converter that transforms Jellyfish's text dump into the OL database format

To clean and rebuild:

```bash
make rebuild
```

## Usage

```
hmm_decoder -p <params> -r <reads> -d <jellyfish_db.jf> [options]
```

### Required Arguments

| Short | Long | Description |
|-------|------|-------------|
| `-p` | `--params <file>` | HMM parameters file (see [Parameter File Format](#parameter-file-format)) |
| `-r` | `--reads <file>` | Input FASTA or FASTQ file containing long reads |
| `-d` | `--database <file>` | Jellyfish k-mer count database (`.jf` file) |

### Optional Arguments

| Short | Long | Description | Default |
|-------|------|-------------|---------|
| `-o` | `--out <file>` | Output results file | `hmm_results.txt` |
| `-t` | `--threads <int>` | Number of threads | all available cores |
| `-b` | `--batch <int>` | Number of reads per processing batch | `1000` |
| `-j` | `--jellyfish <path>` | Path to the `jellyfish` binary | `jellyfish` (via `$PATH`) |
| `-J` | `--jftool <path>` | Path to `JFToOl` converter binary | `./bin/JFToOl` |
| | `--rebuild` | Force rebuild of the OL database even if files already exist | off |
| `-h` | `--help` | Show help message and exit | |

### Examples

Minimal run (auto-detect threads, default output file):

```bash
./bin/hmm_decoder -p params.txt -r reads.fasta -d counts.jf
```

Full options with custom jellyfish path:

```bash
./bin/hmm_decoder -p params.txt -r reads.fasta -d counts.jf \
    -o results.txt -t 32 -b 5000 \
    -j /opt/jellyfish-2.3.0/bin/jellyfish \
    -J ./bin/JFToOl
```

Force OL database rebuild:

```bash
./bin/hmm_decoder -p params.txt -r reads.fasta -d counts.jf -o results.txt -t 32 --rebuild
```

### OL Database Auto-Build

On first run, `hmm_decoder` checks whether the OL database files exist alongside the Jellyfish database:

```
<jellyfish_db>.prefix   — prefix start offsets (uint32, 2^26 + 1 entries)
<jellyfish_db>.suffix   — suffix array (uint16)
<jellyfish_db>.counts   — count array (uint16, capped at 65535)
```

If any are missing, the program automatically builds them by running:

```bash
jellyfish dump -c <db.jf> | JFToOl <db.jf>.prefix <db.jf>.suffix <db.jf>.counts
```

On subsequent runs, the existing OL database is reused. Use `--rebuild` to force regeneration.

## Input Formats

### Parameter File Format

A plain-text file with **one value per line**. The layout is:

```
<num_states>
<initial_prob_state_0>         (log-space)
<initial_prob_state_1>
...
<initial_prob_state_N-1>
<transition(0,0)>              (log-space, row-major order)
<transition(0,1)>
...
<transition(N-1,N-1)>
<emission_mean_state_0>        (natural scale)
<emission_mean_state_1>
...
<emission_mean_state_N-1>
<emission_var_state_0>         (natural scale)
<emission_var_state_1>
...
<emission_var_state_N-1>
```

For the standard 5-state model, this file has 1 + 5 + 25 + 5 + 5 = **41 lines**. All initial and transition probability values must be in **log-space** (i.e., `ln(p)`). Emission means and variances are in natural (linear) scale.

**Example** (`params.txt` for a 5-state model at ~30x coverage):

```
5
-11.5129
-2.3026
-0.1054
-4.6052
-6.9078
-0.1054
-2.3026
-4.6052
-6.9078
-9.2103
-2.3026
-0.1054
-2.3026
-6.9078
-9.2103
-4.6052
-2.3026
-0.1054
-2.3026
-6.9078
-6.9078
-6.9078
-2.3026
-0.1054
-2.3026
-9.2103
-9.2103
-6.9078
-2.3026
-0.1054
0.5
15.0
30.0
60.0
120.0
1.0
25.0
50.0
200.0
800.0
```

### FASTA/FASTQ Input

```
>read_001
ACGTACGTACGTACGTACGT...
>read_002
TGCATGCATGCATGCATGCA...
```

## Output Format

The output is a **tab-separated** file with one row per read. Each row contains 7 fields:

| Column | Field | Description |
|--------|-------|-------------|
| 1 | `read_name` | Read identifier (from FASTA header, without `>`) |
| 2 | `locations` | Heterozygous segment coordinates as `start,end` pairs; semicolon-separated across segments |
| 3 | `het_kmers` | 42-bit encoded k-mer positions within each segment (comma-separated within a segment, semicolon-separated across segments) |
| 4 | `het_probs` | Posterior probability (exp of log-gamma) for each heterozygous k-mer, same order as column 3 |
| 5 | `total_kmers` | Number of k-mers in each segment's detection window; semicolon-separated |
| 6 | `window_prob` | Sigmoid-aggregated window-level heterozygosity confidence (0-1); semicolon-separated |
| 7 | `in_repeat` | Repeat flag per segment: `1` if >50% of window positions were in a repeat region, `0` otherwise |

**Example** (single read with 4 detected heterozygous segments):

```
read_001	108,149;247,291;2479,2519;5427,5476	171628...,2241...;...	0.9826,...;0.1466,...	21;21;21;21	0.7459;0.6528;0.7412;0.5830	1;1;1;1
```

Coordinates in column 2 are **0-indexed base-pair positions** along the read. The `end` coordinate includes the k-1 extension to cover the full span of the last k-mer in the window.

## Project Structure

```
hmm_decoder/
├── src/
│   └── main.cpp                 # Entry point, CLI parsing, OL database build
├── include/
│   ├── data_structures.h        # Core types: SequenceRecord, ReadKmerData, HetLocation, HMMParams
│   ├── matrix.h                 # Row-major dense matrix class
│   ├── math_utils.h             # fast_erf, logsumexp, gradient_coeff, quantize_delta
│   ├── fasta_reader.h           # Streaming FASTA/FASTQ reader with batch callbacks
│   ├── thread_safe_queue.h      # Bounded thread-safe queue with backpressure
│   ├── ol_worker.h/.cpp         # OL database loader and per-read k-mer query engine
│   ├── hmm_functions.h/.cpp     # HMM core: emission, forward-backward, gamma, two-pass Viterbi, het extraction
│   ├── pipeline_workers.h/.cpp  # Stage worker classes (FASTA reader, OL, HMM decoder, output writer)
│   ├── hmm_pipeline.h/.cpp      # Pipeline orchestration and thread management
│   ├── output_writer.h/.cpp     # Tab-separated result formatter
│   └── JFToOl.cpp               # Standalone Jellyfish-to-OL database converter (compiled as separate binary)
├── bin/                         # Compiled binaries (created by make)
├── obj/                         # Object files (created by make)
├── test/                        # Test data
│   ├── test_read.fasta          # Example FASTA input
│   └── out.txt                  # Example output
├── Makefile
└── README.md
```

## License

TBD

## Citation

If you use this software, please cite:

> (Preprint forthcoming)
