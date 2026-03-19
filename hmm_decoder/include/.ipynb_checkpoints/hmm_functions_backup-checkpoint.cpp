#include "hmm_functions.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

void multivariate_gaussian(const std::vector<int>& x,
                          const std::vector<double>& mu,
                          const std::vector<double>& var,
                          Matrix& result,
                          bool integer_mode,
                          double threshold) {
    
    size_t n = x.size();
    size_t m = mu.size();
    result.resize(n, m);

    if (integer_mode) {
        #pragma omp parallel for collapse(2) schedule(static)
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                double x_val = std::min(static_cast<double>(x[i]), threshold);
                double sigma = std::sqrt(2.0 * var[j]);
                double lower = 0.5 * (1.0 + fast_erf((x_val - 0.5 - mu[j]) / sigma));
                double upper = 0.5 * (1.0 + fast_erf((x_val + 0.5 - mu[j]) / sigma));
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

void log_forward(const Matrix& A, const Matrix& B, 
                const std::vector<double>& initial, Matrix& alpha) {
    
    size_t T = B.rows;
    size_t M = A.rows;
    alpha.resize(T, M);

    for (size_t j = 0; j < M; ++j) {
        alpha(0, j) = initial[j] + B(0, j);
    }

    std::vector<double> temp_vals(M);
    for (size_t t = 1; t < T; ++t) {
        for (size_t j = 0; j < M; ++j) {
            for (size_t i = 0; i < M; ++i) {
                temp_vals[i] = A(i, j) + alpha(t-1, i);
            }
            alpha(t, j) = B(t, j) + logsumexp(temp_vals);
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

    std::vector<double> temp_vals(M);
    for (int t = T - 2; t >= 0; --t) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < M; ++j) {
                temp_vals[j] = A(i, j) + B(t+1, j) + beta(t+1, j);
            }
            beta(t, i) = logsumexp(temp_vals);
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

    std::vector<double> log_probs(M);
    for (size_t t = 0; t < T; ++t) {
        for (size_t j = 0; j < M; ++j) {
            log_probs[j] = alpha(t, j) + beta(t, j);
        }
        double log_sum = logsumexp(log_probs);
        
        for (size_t j = 0; j < M; ++j) {
            gamma(t, j) = log_probs[j] - log_sum;
        }
    }
}

std::vector<HetLocation> extract_het_locs(const std::vector<int>& states, 
                                         double f, int k, 
                                         int extrabuffer) {
    
    std::vector<HetLocation> result;
    std::deque<std::tuple<int, int, int>> temp_out;
    
    std::map<int, int> count;
    for (int i = 0; i < 5; ++i) count[i] = 0;
    
    for (int i = 0; i < k - 1 && i < static_cast<int>(states.size()); ++i) {
        count[states[i]]++;
    }
    
    for (size_t i = k - 1; i < states.size(); ++i) {
        if (i > static_cast<size_t>(k - 1)) {
            count[states[i - k]]--;
        }
        count[states[i]]++;
        
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
        }
    }
    
    if (temp_out.empty()) return result;
    
    auto get_het_indices = [&states](int start, int end) {
        std::vector<int> indices;
        for (int i = start; i <= end; ++i) {
            if (states[i] == 1) indices.push_back(i);
        }
        return indices;
    };
    
    auto entry = temp_out.front();
    temp_out.pop_front();
    int buffer = k - std::get<2>(entry);
    
    HetLocation loc;
    loc.start = std::max(0, std::get<0>(entry) - std::max(extrabuffer, buffer));
    loc.end = std::min(std::get<1>(entry) + std::max(extrabuffer, buffer) - buffer, 
                      static_cast<int>(states.size())) + k - 1;
    loc.het_kmer_indices = get_het_indices(std::get<0>(entry), std::get<1>(entry));
    loc.total_kmers = std::get<1>(entry) - std::get<0>(entry) + 1;
    result.push_back(loc);
    
    while (!temp_out.empty()) {
        entry = temp_out.front();
        temp_out.pop_front();
        buffer = k - std::get<2>(entry);
        
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
            result.back().total_kmers += std::get<1>(entry) - std::get<0>(entry) + 1;
        } else {
            HetLocation new_loc;
            new_loc.start = start;
            new_loc.end = end;
            new_loc.het_kmer_indices = new_indices;
            new_loc.total_kmers = std::get<1>(entry) - std::get<0>(entry) + 1;
            result.push_back(new_loc);
        }
    }
    
    return result;
}

void log_viterbi(const HMMParams& params,
                const std::vector<int>& V,
                std::vector<int>& states_out,
                std::vector<double>& probs_out,
                Matrix& gamma_out,
                std::vector<int>& repeat_mask_out) {
    
    size_t T = V.size();
    size_t M = params.num_states;
    
    if (T == 0) throw std::runtime_error("Empty observation sequence");
    
    std::vector<int> dV(T - 1);
    std::vector<double> delta_t(T - 1);
    for (size_t t = 0; t < T - 1; ++t) {
        dV[t] = std::abs(V[t+1] - V[t]);
        delta_t[t] = gradient_coeff(dV[t]);
    }
    
    Matrix B;
    multivariate_gaussian(V, params.emission_mean, params.emission_var, B, true);
    
    Matrix states(T, M);
    Matrix trace(T, M);
    
    for (size_t j = 0; j < M; ++j) {
        states(0, j) = params.initial[j] + B(0, j);
        trace(0, j) = 0;
    }
    
    Matrix Astat(M, M, -std::numeric_limits<double>::infinity());
    for (size_t i = 0; i < M; ++i) {
        Astat(i, i) = 0.0;
    }
    
    std::unordered_map<int, Matrix> At_cache;
    constexpr int MAX_CACHE_SIZE = 1024;
    
    double mean_threshold = 0.15 * params.emission_mean[2];
    
    // First pass
    for (size_t t = 1; t < T; ++t) {
        Matrix A_current;
        
        if (dV[t-1] > mean_threshold) {
            int qd = quantize_delta(delta_t[t-1]);
            auto it = At_cache.find(qd);
            if (it != At_cache.end()) {
                A_current = it->second;
            } else {
                compute_logAt_fast(params.transition, delta_t[t-1], A_current);
                if (At_cache.size() < MAX_CACHE_SIZE) {
                    At_cache[qd] = A_current;
                }
            }
        } else {
            A_current = Astat;
        }
        
        for (size_t j = 0; j < M; ++j) {
            double max_val = -std::numeric_limits<double>::infinity();
            int max_idx = 0;
            
            for (size_t i = 0; i < M; ++i) {
                double val = A_current(i, j) + states(t-1, i);
                if (val > max_val) {
                    max_val = val;
                    max_idx = i;
                }
            }
            
            states(t, j) = max_val + B(t, j);
            trace(t, j) = max_idx;
        }
    }
    
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
    
    // Second pass
    states.resize(T, M, 0.0);
    trace.resize(T, M, 0.0);
    
    std::vector<int> repeat_mask(T, 0);
    
    int w = 1000;
    int count = 0;
    int total = 0;
    int count_repeat = 0;
    
    for (int i = 0; i < std::min(w, static_cast<int>(T)); ++i) {
        if (out[i] == 2 || out[i] == 1) {
            count++;
            total += V[i];
        }
        if (out[i] == 3 || out[i] == 4) {
            count_repeat++;
        }
    }
    
    for (size_t j = 0; j < M; ++j) {
        states(0, j) = params.initial[j] + B(0, j);
        trace(0, j) = 0;
    }
    
    double temp_mu = params.emission_mean[2];
    
    for (size_t t = 1; t < T; ++t) {
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
        
        if (count >= 0.5 * std::min(w, static_cast<int>(T))) {
            temp_mu = std::ceil(static_cast<double>(total) / count);
            
            double x_val = std::min(static_cast<double>(V[t]), 1000.0);
            std::vector<double> local_mu = {temp_mu / 2.0, temp_mu};
            std::vector<double> local_var = {
                params.emission_var[1], 
                params.emission_var[2]
            };
            
            for (int j = 0; j < 2; ++j) {
                double sigma = std::sqrt(2.0 * local_var[j]);
                double lower = 0.5 * (1.0 + fast_erf((x_val - 0.5 - local_mu[j]) / sigma));
                double upper = 0.5 * (1.0 + fast_erf((x_val + 0.5 - local_mu[j]) / sigma));
                double p = std::max(upper - lower, 1e-300);
                B(t, j + 1) = std::log(p);
            }
        } else if (count_repeat >= 0.5 * std::min(w, static_cast<int>(T))) {
            repeat_mask[t] = 1;
        }
        
        Matrix A_current;
        if (dV[t-1] > 0.15 * temp_mu) {
            int qd = quantize_delta(delta_t[t-1]);
            auto it = At_cache.find(qd);
            if (it != At_cache.end()) {
                A_current = it->second;
            } else {
                compute_logAt_fast(params.transition, delta_t[t-1], A_current);
                if (At_cache.size() < MAX_CACHE_SIZE) {
                    At_cache[qd] = A_current;
                }
            }
        } else {
            A_current = Astat;
        }
        
        for (size_t j = 0; j < M; ++j) {
            double max_val = -std::numeric_limits<double>::infinity();
            int max_idx = 0;
            
            for (size_t i = 0; i < M; ++i) {
                double val = A_current(i, j) + states(t-1, i);
                if (val > max_val) {
                    max_val = val;
                    max_idx = i;
                }
            }
            
            states(t, j) = max_val + B(t, j);
            trace(t, j) = max_idx;
        }
    }
    
    compute_gamma(params.transition, B, params.initial, gamma_out);
    
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
