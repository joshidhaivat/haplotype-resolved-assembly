#ifndef HMM_FUNCTIONS_H
#define HMM_FUNCTIONS_H

#include <vector>
#include <deque>
#include <map>
#include <tuple>
#include <unordered_map>
#include <limits>
#include "matrix.h"
#include "data_structures.h"
#include "math_utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

inline double compute_window_prob(const Matrix& gamma, const std::vector<int>& indices,
                                  const double sigmoid_factor = 0.1, const double threshold = 0.95) {
    /*  =============  OLD METHOD  ====================
        std::vector<double> state1_probs;
        for (int idx : indices) {
            state1_probs.push_back(gamma(idx, 1));
        }
        double log_sum_state1 = logsumexp(state1_probs);
        
        std::vector<double> all_probs;
        for (int idx : indices) {
            for (size_t s = 0; s < params_.num_states; ++s) {
                all_probs.push_back(gamma(idx, s));
            }
        }
        double log_sum_all = logsumexp(all_probs);
        
        window_prob = std::exp(log_sum_state1 - log_sum_all);
    */
    double sum = 0.0;
    double logp = 0.0;
    double upper_threshold = std::log(threshold);
    double lower_threshold = std::log(1 - threshold);
    
    for (int idx : indices) {
        logp = gamma(idx, 1);

        // Inline logit_from_logp logic
        if (logp <= lower_threshold) {
            sum += (lower_threshold - upper_threshold);
        }
        else if (logp >= upper_threshold) {
            sum += (upper_threshold - lower_threshold);
        }
        else {
            // Stable log(p / (1 - p))
            sum += logp - std::log(-std::expm1(logp));
        }
    }

    // Stable sigmoid
    sum = sigmoid_factor * sum;
    if (sum >= 0.0) {
        return 1.0 / (1.0 + std::exp(-sum));
    } else {
        double e = std::exp(sum);
        return e / (1.0 + e);
    }
}

inline double logsumexp_fixed5(const double* vals);

void multivariate_gaussian(const std::vector<unsigned long int>& x,
                          const std::vector<double>& mu,
                          const std::vector<double>& var,
                          Matrix& result,
                          bool integer_mode = true,
                          double threshold = 1000.0);

void compute_logAt_fast(const Matrix& logA, double delta_t, Matrix& out);

void log_forward(const Matrix& A, const Matrix& B, 
                const std::vector<double>& initial, Matrix& alpha);

void log_backward(const Matrix& A, const Matrix& B, Matrix& beta);

void compute_gamma(const Matrix& A, const Matrix& B,
                  const std::vector<double>& initial, Matrix& gamma);

std::vector<HetLocation> extract_het_locs(const std::vector<int>& states, 
                                         double f = 0.5, int k = 21, 
                                         int extrabuffer = 0);

void log_viterbi(const HMMParams& params,
                const std::vector<unsigned long int>& V,
                std::vector<int>& states_out,
                std::vector<double>& probs_out,
                Matrix& gamma_out,
                std::vector<int>& repeat_mask_out);

#endif // HMM_FUNCTIONS_H