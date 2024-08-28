#ifndef SPARSE_MATRIX_UTIL_H
#define SPARSE_MATRIX_UTIL_H

#include "sparse_matrix_base.hpp"

#include <iostream>
#include <sstream>
#include <vector>

namespace ldpc {
namespace sparse_matrix_util {

template <class SPARSE_MATRIX_CLASS>
std::stringstream print_sparse_matrix(SPARSE_MATRIX_CLASS& matrix, bool SILENT = false) {
    std::stringstream ss;
    int               m = matrix.m;
    int               n = matrix.n;
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            //     // std::cout<<j<<" "<<i<<std::endl;
            auto e = matrix.get_entry(j, i);

            if (e.at_end()) {
                ss << unsigned(0);
            } else {
                ss << e.str();
                // std::cout<<e->row_index<<" "<<e->col_index<<" "<<unsigned(e->value)<<std::endl;
            }

            if (i != (n - 1)) {
                ss << " ";
            }
        }
        if (j != (m - 1)) {
            ss << "\n";
        }
    }
    if (!SILENT) {
        std::cout << ss.str() << std::endl;
    }
    return ss;
}

template <class SPARSE_MATRIX_CLASS>
SPARSE_MATRIX_CLASS copy_cols(SPARSE_MATRIX_CLASS& mat, const std::vector<int>& cols) {
    int m              = 0;
    int n              = 0;
    int i              = 0;
    int j              = 0;
    m                  = mat.m;
    n                  = cols.size();
    auto copy_mat      = SPARSE_MATRIX_CLASS(m, n);
    int  new_col_index = -1;
    for (auto col_index : cols) {
        new_col_index += 1;
        for (auto& e : mat.iterate_column(col_index)) {
            copy_mat.insert_entry(e.row_index, new_col_index, e.value);
        }
    }
    return copy_mat;
}

template <class T>
void print_vector(const T& input) {
    int length = input.size();
    std::cout << "[";
    for (int i = 0; i < length; i++) {
        if (std::is_same<T, std::vector<uint8_t>>::value) {
            std::cout << unsigned(input[i]);
        } else {
            std::cout << input[i];
        }
        if (i != (length - 1)) {
            std::cout << " ";
        }
    }
    std::cout << "]" << std::endl;
}

template <class T>
void print_array(T array, int length) {
    for (int i = 0; i < length; i++) {
        if (std::is_same<T, uint8_t*>::value) {
            std::cout << unsigned(array[i]) << " ";
        } else {
            std::cout << array[i] << " ";
        }
    }
    std::cout << std::endl;
}

} // namespace sparse_matrix_util
} // namespace ldpc

#endif
