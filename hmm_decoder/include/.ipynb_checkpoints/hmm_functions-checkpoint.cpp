#include "hmm_functions.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <cstring>
//#include "debug_logger.h"

// ============================================================================
// OPTIMIZATION 1: Inline logsumexp for fixed M=5
// ============================================================================
// Before: Generic vector-based logsumexp with iterator overhead
// After: Unrolled loop for exactly 5 states (no iterator, no branching)

inline double logsumexp_fixed5(const double* vals) {
    // Find max (unrolled)
    double max_val = vals[0];
    max_val = (vals[1] > max_val) ? vals[1] : max_val;
    max_val = (vals[2] > max_val) ? vals[2] : max_val;
    max_val = (vals[3] > max_val) ? vals[3] : max_val;
    max_val = (vals[4] > max_val) ? vals[4] : max_val;
    
    if (std::isinf(max_val) && max_val < 0) return max_val;
    
    // Sum (unrolled, no loop)
    double sum = std::exp(vals[0] - max_val)
               + std::exp(vals[1] - max_val)
               + std::exp(vals[2] - max_val)
               + std::exp(vals[3] - max_val)
               + std::exp(vals[4] - max_val);
    
    return max_val + std::log(sum);
}

// ============================================================================
// OPTIMIZATION 2: Pre-compute reciprocals (multiply faster than divide)
// ============================================================================

void multivariate_gaussian(const std::vector<unsigned long int>& x,
                          const std::vector<double>& mu,
                          const std::vector<double>& var,
                          Matrix& result,
                          bool integer_mode,
                          double threshold) {
    
    size_t n = x.size();
    size_t m = mu.size();
    result.resize(n, m);

    if (integer_mode) {
        // Pre-calculate constants (moved outside loops)
        double inv_sqrt_2var[5];
        for (size_t j = 0; j < m; ++j) {
            inv_sqrt_2var[j] = 1.0 / std::sqrt(2.0 * var[j]);
        }
        
        #pragma omp parallel for collapse(2) schedule(static)
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                double x_val = std::min(static_cast<double>(x[i]), threshold);
                // Multiply instead of divide (faster)
                double lower = 0.5 * (1.0 + fast_erf((x_val - 0.5 - mu[j]) * inv_sqrt_2var[j]));
                double upper = 0.5 * (1.0 + fast_erf((x_val + 0.5 - mu[j]) * inv_sqrt_2var[j]));
                double p = std::max(upper - lower, 1e-300);
                result(i, j) = std::log(p);
            }
        }
    }
}

void compute_logAt_fast(const Matrix& logA, double delta_t, Matrix& out) {
    size_t m = logA.rows;
    out.resize(m, m);

    double log_delta = std::log(delta_t);
    double log1m_delta = std::log1p(-delta_t);
    double logm = std::log(static_cast<double>(m - 1));

    std::vector<double> diag = logA.get_diagonal();

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < m; ++i) {
        double add_term = diag[i] + log_delta - logm;
        for (size_t j = 0; j < m; ++j) {
            if (i == j) {
                out(i, j) = diag[i] + log1m_delta;
            } else {
                out(i, j) = logsumexp2(logA(i, j), add_term);
            }
        }
    }
}

// ============================================================================
// OPTIMIZATION 3: Stack arrays instead of heap vectors for small fixed sizes
// ============================================================================

void log_forward(const Matrix& A, const Matrix& B, 
                const std::vector<double>& initial, Matrix& alpha) {
    
    size_t T = B.rows;
    size_t M = A.rows;
    alpha.resize(T, M);

    for (size_t j = 0; j < M; ++j) {
        alpha(0, j) = initial[j] + B(0, j);
    }

    // Use stack array instead of vector (faster allocation)
    double temp_vals[5];  // M=5 is constant
    
    for (size_t t = 1; t < T; ++t) {
        for (size_t j = 0; j < M; ++j) {
            for (size_t i = 0; i < M; ++i) {
                temp_vals[i] = A(i, j) + alpha(t-1, i);
            }
            alpha(t, j) = B(t, j) + logsumexp_fixed5(temp_vals);
        }
    }
}

void log_backward(const Matrix& A, const Matrix& B, Matrix& beta) {
    size_t T = B.rows;
    size_t M = A.rows;
    beta.resize(T, M);

    for (size_t j = 0; j < M; ++j) {
        beta(T-1, j) = 0.0;
    }

    double temp_vals[5];  // Stack array
    
    for (int t = T - 2; t >= 0; --t) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < M; ++j) {
                temp_vals[j] = A(i, j) + B(t+1, j) + beta(t+1, j);
            }
            beta(t, i) = logsumexp_fixed5(temp_vals);
        }
    }
}

void compute_gamma(const Matrix& A, const Matrix& B,
                  const std::vector<double>& initial, Matrix& gamma) {
    
    Matrix alpha, beta;
    log_forward(A, B, initial, alpha);
    log_backward(A, B, beta);

    size_t T = alpha.rows;
    size_t M = alpha.cols;
    gamma.resize(T, M);

    double log_probs[5];  // Stack array
    
    for (size_t t = 0; t < T; ++t) {
        for (size_t j = 0; j < M; ++j) {
            log_probs[j] = alpha(t, j) + beta(t, j);
        }
        double log_sum = logsumexp_fixed5(log_probs);
        
        for (size_t j = 0; j < M; ++j) {
            gamma(t, j) = log_probs[j] - log_sum;
        }
    }
}

std::vector<HetLocation> extract_het_locs(const std::vector<int>& states, 
                                         double f, int k, 
                                         int extrabuffer) {

    //LOG_FUNC_ENTER();
    
    std::vector<HetLocation> result;
    std::deque<std::tuple<int, int, int>> temp_out;
    int count[5];
    
    for (int i = 0; i < 5; ++i) count[i] = 0;
    
    for (int i = 0; i < k - 1 && i < static_cast<int>(states.size()); ++i) {
        count[states[i]]++;
    }
    
    for (size_t i = k - 1; i < states.size(); ++i) {
        if (i > static_cast<size_t>(k - 1)) {
            count[states[i - k]]--;
        }
        count[states[i]]++;

        //LOG_STREAM("i = " << i << "/" << states.size() << " count = " << count[0] << ", " << count[1] << ", " 
        //            << count[2] << ", " << count[3] << ", " << count[4]);
        
        if (count[1] >= f * k) {
            int start = i - k + 1;
            int end = i;
            int c = count[1];
            
            if (!temp_out.empty() && start <= std::get<1>(temp_out.back())) {
                auto t = temp_out.back();
                temp_out.pop_back();
                if (states[std::get<0>(t)] == 1) {
                    start = std::get<0>(t);
                    if (states[i] != 1) {
                        end = std::get<1>(t);
                        c = std::get<2>(t);
                    } else {
                        c = k;
                    }
                }
            }
            temp_out.push_back(std::make_tuple(start, end, c));
            //LOG_STREAM("last added to temp_out : (" << start << ", " << end << ", " << c << ")");
        }
    }
    
    if (temp_out.empty()) return result;
    
    auto get_het_indices = [&states](int start, int end) {
        std::vector<int> indices;
        indices.reserve(end - start + 1);  // Pre-allocate
        for (int i = start; i <= end; ++i) {
            if (states[i] == 1) indices.push_back(i);
        }
        return indices;
    };
    
    auto entry = temp_out.front();
    temp_out.pop_front();
    int buffer = std::max(0, k - std::get<2>(entry));

    //LOG_STREAM("Entry: (" << std::get<0>(entry) << ", " << std::get<1>(entry) << ", " << std::get<2>(entry) << ")");
    
    HetLocation loc;
    loc.start = std::max(0, std::get<0>(entry) - std::max(extrabuffer, buffer));
    loc.end = std::min(std::get<1>(entry) + std::max(extrabuffer, buffer) - buffer, 
                      static_cast<int>(states.size())) + k - 1;
    loc.het_kmer_indices = get_het_indices(std::get<0>(entry), std::get<1>(entry));
    loc.kmer_start = std::get<0>(entry) - buffer;
    loc.kmer_end = std::get<1>(entry);
    loc.total_kmers = std::get<1>(entry) - std::get<0>(entry) + 1;
    result.push_back(loc);

    //LOG_STREAM("last element of result : (" << result.back().start << ", " << result.back().end << ", "
    //    << result.back().kmer_start << ", " << result.back().kmer_end << ", "
    //    << result.back().total_kmers << ")");
    
    while (!temp_out.empty()) {
        entry = temp_out.front();
        temp_out.pop_front();
        buffer = std::max(0, k - std::get<2>(entry));

        //LOG_STREAM("Entry: (" << std::get<0>(entry) << ", " << std::get<1>(entry) << ", " << std::get<2>(entry) << ")");
        
        int start = std::max(0, std::get<0>(entry) - std::max(extrabuffer, buffer));
        int end = std::min(std::get<1>(entry) + std::max(extrabuffer, buffer) - buffer, 
                          static_cast<int>(states.size())) + k - 1;
        
        auto new_indices = get_het_indices(std::get<0>(entry), std::get<1>(entry));
        
        if (start <= result.back().end) {
            result.back().end = end;
            result.back().het_kmer_indices.insert(
                result.back().het_kmer_indices.end(),
                new_indices.begin(), new_indices.end()
            );
            result.back().kmer_start = std::min(std::get<0>(entry) - buffer, result.back().kmer_start);
            result.back().kmer_end = std::max(std::get<1>(entry), result.back().kmer_end);
            result.back().total_kmers += std::get<1>(entry) - std::get<0>(entry) + 1;
        } else {
            HetLocation new_loc;
            new_loc.start = start;
            new_loc.end = end;
            new_loc.het_kmer_indices = new_indices;
            new_loc.kmer_start = std::get<0>(entry) - buffer;
            new_loc.kmer_end = std::get<1>(entry);
            new_loc.total_kmers = std::get<1>(entry) - std::get<0>(entry) + 1;
            result.push_back(new_loc);
        }
        //LOG_STREAM("last element of result : (" << result.back().start << ", " << result.back().end << ", "
        //    << result.back().kmer_start << ", " << result.back().kmer_end << ", "
        //    << result.back().total_kmers << ")");
    }

    //LOG_FUNC_EXIT();
    
    return result;
}

// ============================================================================
// OPTIMIZATION 4: Direct array cache instead of unordered_map
// OPTIMIZATION 5: Reuse matrices between passes
// ============================================================================

void log_viterbi(const HMMParams& params,
                const std::vector<unsigned long int>& V,
                std::vector<int>& states_out,
                std::vector<double>& probs_out,
                Matrix& gamma_out,
                std::vector<int>& repeat_mask_out) {
    
    size_t T = V.size();
    size_t M = params.num_states;
    
    if (T == 0) throw std::runtime_error("Empty observation sequence");
    
    // ========================================================================
    // Pre-compute delta values (shared by both passes)
    // ========================================================================
    std::vector<int> dV(T - 1);
    std::vector<double> delta_t(T - 1);
    for (size_t t = 0; t < T - 1; ++t) {
        dV[t] = std::max(V[t+1], V[t]) - std::min(V[t+1], V[t]);
        delta_t[t] = gradient_coeff(dV[t]);
    }
    
    // ========================================================================
    // Calculate initial emission probabilities (will be modified in pass 2)
    // ========================================================================
    Matrix B;
    multivariate_gaussian(V, params.emission_mean, params.emission_var, B, true);
    
    // ========================================================================
    // Allocate matrices once (reused in both passes)
    // ========================================================================
    Matrix states(T, M);
    Matrix trace(T, M);
    
    // ========================================================================
    // Static transition matrix for small deltas (shared)
    // ========================================================================
    Matrix Astat(M, M, -std::numeric_limits<double>::infinity());
    for (size_t i = 0; i < M; ++i) {
        Astat(i, i) = 0.0;
    }
    
    // ========================================================================
    // OPTIMIZATION: Direct-indexed cache instead of unordered_map
    // ========================================================================
    constexpr int CACHE_SIZE = 1024;
    std::vector<Matrix> At_cache(CACHE_SIZE);
    std::vector<bool> cache_valid(CACHE_SIZE, false);
    
    double mean_threshold = 0.15 * params.emission_mean[2];
    
    // ========================================================================
    // Helper lambda for Viterbi pass (reused twice)
    // ========================================================================
    auto viterbi_pass = [&](bool is_second_pass) {
        // Initialize
        for (size_t j = 0; j < M; ++j) {
            states(0, j) = params.initial[j] + B(0, j);
            trace(0, j) = 0;
        }
        
        // double temp_vals[5];  // Stack array for M=5
        
        for (size_t t = 1; t < T; ++t) {
            // Get transition matrix
            const Matrix* A_current;
            if (dV[t-1] > mean_threshold) {
                int qd = quantize_delta(delta_t[t-1]);
                // Clamp to valid cache range
                qd = std::min(std::max(qd, 0), CACHE_SIZE - 1);
                
                if (!cache_valid[qd]) {
                    compute_logAt_fast(params.transition, delta_t[t-1], At_cache[qd]);
                    cache_valid[qd] = true;
                }
                A_current = &At_cache[qd];
            } else {
                A_current = &Astat;
            }
            
            // Viterbi step
            for (size_t j = 0; j < M; ++j) {
                double max_val = -std::numeric_limits<double>::infinity();
                int max_idx = 0;
                
                for (size_t i = 0; i < M; ++i) {
                    double val = (*A_current)(i, j) + states(t-1, i);
                    if (val > max_val) {
                        max_val = val;
                        max_idx = i;
                    }
                }
                
                states(t, j) = max_val + B(t, j);
                trace(t, j) = max_idx;
            }
        }
    };
    
    // ========================================================================
    // FIRST PASS: Initial Viterbi
    // ========================================================================
    viterbi_pass(false);
    
    // Backtrack
    std::vector<int> out(T);
    double max_val = states(T-1, 0);
    int max_idx = 0;
    for (size_t j = 1; j < M; ++j) {
        if (states(T-1, j) > max_val) {
            max_val = states(T-1, j);
            max_idx = j;
        }
    }
    out[T-1] = max_idx;
    
    for (int t = T - 2; t >= 0; --t) {
        out[t] = static_cast<int>(trace(t+1, out[t+1]));
    }
    
    // ========================================================================
    // SECOND PASS: Adjust mu based on first pass results
    // ========================================================================
    
    // Sliding window setup
    int w = 1000;
    int count = 0;
    int total = 0;
    int count_repeat = 0;
    std::vector<int> repeat_mask(T, 0);
    bool update_B = false;
    
    for (int i = 0; i < std::min(w, static_cast<int>(T)); ++i) {
        if (out[i] == 2 || out[i] == 1) {
            count++;
            total += V[i];
        }
        if (out[i] == 3 || out[i] == 4) {
            count_repeat++;
        }
    }
    
    double temp_mu = params.emission_mean[2];
    double factor[4];
    for (double& x : factor){
        x = 1.0;
    }
    
    // Update B matrix based on sliding window analysis
    for (size_t t = 1; t < T; ++t) {
        
        temp_mu = params.emission_mean[2]; // Initialize each time to global mean
        
        if (t > static_cast<size_t>(w/2) && T - t > static_cast<size_t>(w/2)) {
            int remove_idx = t - w/2 - 1;
            int add_idx = t + w/2 - 1;
            
            if (out[remove_idx] == 2 || out[remove_idx] == 1) {
                count--;
                total -= V[remove_idx];
            } else if (out[remove_idx] == 3 || out[remove_idx] == 4) {
                count_repeat--;
            }
            
            if (out[add_idx] == 2 || out[add_idx] == 1) {
                count++;
                total += V[add_idx];
            } else if (out[add_idx] == 3 || out[add_idx] == 4) {
                count_repeat++;
            }
        }

        update_B = false;
        
        if (count >= 0.5 * std::min(w, static_cast<int>(T))) {
            temp_mu = std::ceil(static_cast<double>(total) / count);
            factor[0] = 0.5;
            // factor[1] = 1.0;
            factor[2] = 1.0;
            factor[3] = 1.0;
            update_B = true;
        } else if (count_repeat >= 0.5 * std::min(w, static_cast<int>(T))) {
            repeat_mask[t] = 1;
            factor[0] = 0.4;
            // factor[1] = 1.0;
            factor[2] = 0.8;
            factor[3] = 0.8;
            update_B = true;
        }

        if (update_B) {
            double x_val = std::min(static_cast<double>(V[t]), 1000.0);
            std::vector<double> local_mu = {std::floor(factor[0]*temp_mu), std::floor(factor[1]*temp_mu)};
            std::vector<double> local_var = {
                factor[2]*params.emission_var[1], 
                factor[3]*params.emission_var[2]
            };
            
            // Pre-compute reciprocal
            double inv_sqrt_2var[2];
            for (int j = 0; j < 2; ++j) {
                inv_sqrt_2var[j] = 1.0 / std::sqrt(2.0 * local_var[j]);
            }
            
            for (int j = 0; j < 2; ++j) {
                double lower = 0.5 * (1.0 + fast_erf((x_val - 0.5 - local_mu[j]) * inv_sqrt_2var[j]));
                double upper = 0.5 * (1.0 + fast_erf((x_val + 0.5 - local_mu[j]) * inv_sqrt_2var[j]));
                double p = std::max(upper - lower, 1e-300);
                B(t, j + 1) = std::log(p);
            }
        }
    }
    
    // ========================================================================
    // Run Viterbi again with adjusted B matrix
    // ========================================================================
    viterbi_pass(true);
    
    // ========================================================================
    // Compute gamma (using optimized version)
    // ========================================================================
    compute_gamma(params.transition, B, params.initial, gamma_out);
    
    // ========================================================================
    // Final backtrack
    // ========================================================================
    max_val = states(T-1, 0);
    max_idx = 0;
    for (size_t j = 1; j < M; ++j) {
        if (states(T-1, j) > max_val) {
            max_val = states(T-1, j);
            max_idx = j;
        }
    }
    out[T-1] = max_idx;
    
    std::vector<double> probs(T);
    probs[T-1] = gamma_out(T-1, max_idx);
    
    for (int t = T - 2; t >= 0; --t) {
        out[t] = static_cast<int>(trace(t+1, out[t+1]));
        probs[t] = gamma_out(t, out[t]);
    }
    
    states_out = out;
    probs_out = probs;
    repeat_mask_out = repeat_mask;
}
