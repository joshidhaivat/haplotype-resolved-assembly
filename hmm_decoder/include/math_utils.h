#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

// Fast erf approximation (Abramowitz and Stegun)
inline double fast_erf(double x) {
    constexpr double a1 = 0.254829592;
    constexpr double a2 = -0.284496736;
    constexpr double a3 = 1.421413741;
    constexpr double a4 = -1.453152027;
    constexpr double a5 = 1.061405429;
    constexpr double p = 0.3275911;

    int sign = (x < 0) ? -1 : 1;
    x = std::abs(x);
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);
    return sign * y;
}

// Log-sum-exp for numerical stability
inline double logsumexp(const std::vector<double>& values) {
    if (values.empty()) return -std::numeric_limits<double>::infinity();
    
    double max_val = *std::max_element(values.begin(), values.end());
    if (std::isinf(max_val) && max_val < 0) return max_val;
    
    double sum = 0.0;
    for (double v : values) {
        sum += std::exp(v - max_val);
    }
    return max_val + std::log(sum);
}

inline double logsumexp2(double a, double b) {
    double max_val = std::max(a, b);
    if (std::isinf(max_val)) return max_val;
    return max_val + std::log(std::exp(a - max_val) + std::exp(b - max_val));
}

// Gradient coefficient
inline double gradient_coeff(double x, double alpha = 0.5) {
    return 0.9 * ((2.0 / (1.0 + std::exp(-alpha * std::abs(x)))) - 1.0);
}

// Quantize delta for caching
inline int quantize_delta(double delta, double q = 1e-4) {
    return static_cast<int>(std::round(delta / q));
}

#endif
