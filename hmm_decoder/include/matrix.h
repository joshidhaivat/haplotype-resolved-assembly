#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <algorithm>

class Matrix {
public:
    std::vector<double> data;
    size_t rows, cols;

    Matrix(size_t r = 0, size_t c = 0, double init_val = 0.0) 
        : data(r * c, init_val), rows(r), cols(c) {}

    inline double& operator()(size_t i, size_t j) {
        return data[i * cols + j];
    }

    inline const double& operator()(size_t i, size_t j) const {
        return data[i * cols + j];
    }

    void resize(size_t r, size_t c, double val = 0.0) {
        rows = r;
        cols = c;
        data.assign(r * c, val);
    }

    std::vector<double> get_diagonal() const {
        size_t n = std::min(rows, cols);
        std::vector<double> diag(n);
        for (size_t i = 0; i < n; ++i) {
            diag[i] = data[i * cols + i];
        }
        return diag;
    }
};

#endif // MATRIX_H
