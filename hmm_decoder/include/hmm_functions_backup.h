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

void multivariate_gaussian(const std::vector<int>& x,
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
                const std::vector<int>& V,
                std::vector<int>& states_out,
                std::vector<double>& probs_out,
                Matrix& gamma_out,
                std::vector<int>& repeat_mask_out);

#endif // HMM_FUNCTIONS_H
