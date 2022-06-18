//
// Created by luca on 09/06/22.
//
#ifndef QUNIONFIND_UTILS_HPP
#define QUNIONFIND_UTILS_HPP
#include "TreeNode.hpp"

#include <NTL/GF2.h>
#include <NTL/mat_GF2.h>
#include <cassert>
#include <iostream>
#include <ostream>
#include <random>
#include <set>
#include <vector>

class Utils {
public:
    /**
     * Computes and returns the matrix obtained by appending the column vector to the input matrix, result = (matrix|vector)
     * @param matrix
     * @param vector
     * @return
     */
    static std::vector<std::vector<bool>> getAugmentedMatrix(const std::vector<std::vector<bool>>& matrix, const std::vector<bool>& vector) {
        std::vector<std::vector<bool>> result(matrix.size());

        for (size_t i = 0; i < matrix.size(); i++) {
            result.at(i) = std::vector<bool>(matrix.at(i).size() + 1);
            for (std::size_t j = 0; j < matrix.at(0).size(); j++) {
                result[i][j] = matrix[i][j];
            }
            result[i][matrix.at(0).size()] = vector.at(i);
        }
        return result;
    }

    /**
     * Checks if the give matrix (assumed to be square) is invertible by a determinant computation
     * @param matrix
     * @return
     */
    static bool isInvertible(const std::vector<std::vector<bool>>& matrix) {
        if (!matrix.empty()) {
            if (matrix.size() != matrix.at(0).size()) {
                return false;
            }
        }
        return computeDeterminant(matrix) == 0;
    }

    /**
     * Standard Determinant computation
     * @param matrix
     * @return
     */
    static std::size_t computeDeterminant(const std::vector<std::vector<bool>>& matrix) {
        std::size_t result = 0U;
        std::size_t n      = matrix.size();
        if (n == 1) {
            return matrix[0][0];
        }
        if (n == 2) {
            return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
        }
        std::vector<std::vector<bool>> tmp;
        bool                           sign = 1;
        for (std::size_t i = 0; i < n; i++) {
            tmp = getSubMatrix(matrix, 0, i);
            result += sign && matrix[0][i] && computeDeterminant(tmp);
            sign = !sign;
        }
        return result;
    }

    /**
     * Returns the submatrix where the row and column given by the index parameters are erased
     * @param matrix
     * @param rowIdx
     * @param colIdx
     * @return
     */
    static std::vector<std::vector<bool>> getSubMatrix(const std::vector<std::vector<bool>>& matrix, const std::size_t rowIdx, const std::size_t colIdx) {
        std::size_t                    rCnt = 0, cCnt = 0;
        std::size_t                    n = matrix.size();
        std::vector<std::vector<bool>> result(n);

        for (std::size_t row = 0; row < n; row++) {
            std::vector<bool> tmp(n);
            for (std::size_t col = 0; col < n; col++) {
                if (row != rowIdx && col != colIdx) {
                    tmp[cCnt++] = matrix[row][col];
                }
            }
            result.at(rCnt++) = tmp;
        }
        return result;
    }

    /**
     * Solves linear equation mod 2 given by Mx = vec
     * @param M
     * @param vec
     * @return solution x to Mx=vec (mod 2)
     */
    static std::vector<bool> solveSystem(const std::vector<std::vector<bool>>& M, const std::vector<bool>& vec) {
        std::vector<bool>              result(M.at(0).size());
        auto                           n = vec.size();
        std::vector<std::vector<bool>> matrix(n);
        if (M.size() < M.at(0).size()) { // if system is overdetermined we only consider the first nxn block n=dim(vec) of the matrix
            for (size_t i = 0; i < n; i++) {
                std::vector<bool> tmp(n);
                for (size_t j = 0; j < n; j++) {
                    tmp[j] = M[i][j];
                }
                matrix.at(i) = tmp;
            }
        } else {
            matrix = M;
        }
        auto         matr = getNtlMatrix(matrix);
        NTL::vec_GF2 x;
        x.SetLength(matrix.size());
        NTL::vec_GF2 b;
        b.SetLength(vec.size());
        for (size_t i = 0; i < vec.size(); i++) {
            b[i] = vec[i];
        }
        NTL::GF2 d;
        NTL::solve(d, matr, x, b);
        for (size_t i = 0; i < x.length(); i++) {
            if (x[i] == 1) {
                result[i] = 1;
            } else {
                result[i] = 0;
            }
        }
        return result;
    }

    /**
     * Checks if the given vector is in the rowspace of matrix M
     * @param M
     * @param vec
     * @return
     */
    static bool isVectorInRowspace(const std::vector<std::vector<bool>>& M, const std::vector<bool>& vec) {
        auto matrix = getTranspose(M);
        auto augm   = getAugmentedMatrix(matrix, vec);
        matrix      = gauss(augm);
        printGF2matrix(matrix);
        std::vector<bool> vector(vec.size());

        for (size_t i = 0; i < matrix.size(); i++) {
            vector.at(i) = matrix[i][matrix.at(i).size() - 1];
        }
        printGF2vector(vector);
        // check consistency
        for (size_t i = 0; i < vector.size(); i++) {
            if (vector[i]) {
                for (size_t j = 0; j < matrix.at(i).size(); j++) {
                    if (std::none_of(matrix.at(i).begin(), matrix.at(i).end() - 1, [](const bool val) { return val; })) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    static NTL::Mat<NTL::GF2> getNtlMatrix(const std::vector<std::vector<bool>>& matrix) {
        NTL::Mat<NTL::GF2> result;
        result.SetDims(matrix.size(), matrix.at(0).size());
        for (long i = 0; i < matrix.size(); i++) {
            for (long j = 0; j < matrix.at(0).size(); j++) {
                if (matrix[i][j]) {
                    result[i][j] = 1;
                } else {
                    result[i][j] = 0;
                }
            }
        }
        return result;
    }

    static std::vector<std::vector<bool>> getMatrixFromNtl(const NTL::Mat<NTL::GF2>& matrix) {
        std::vector<std::vector<bool>> result(matrix.NumRows());

        for (long i = 0; i < matrix.NumRows(); i++) {
            result.at(i) = std::vector<bool>(matrix.NumCols());
            for (long j = 0; j < matrix.NumCols(); j++) {
                if (matrix[i][j] == 1) {
                    result[i][j] = true;
                } else {
                    result[i][j] = false;
                }
            }
        }
        return result;
    }

    // produces row echelon form
    static std::vector<std::vector<bool>> gauss(const std::vector<std::vector<bool>>& matrix) {
        auto mat = getNtlMatrix(matrix);
        NTL::gauss(mat);
        return getMatrixFromNtl(mat);
    }

    /**
     * Computes the transpose of the given matrix
     * @param matrix
     * @return
     */
    static std::vector<std::vector<bool>> getTranspose(const std::vector<std::vector<bool>>& matrix) {
        std::vector<std::vector<bool>> transp(matrix.at(0).size());
        for (auto& i: transp) {
            i = std::vector<bool>(matrix.size());
        }
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix.at(i).size(); j++) {
                transp[j][i] = matrix[i][j];
            }
        }
        return transp;
    }

    /**
     * Standard matrix multiplication
     * @param m1
     * @param m2
     * @return
     */
    static std::vector<std::vector<bool>> rectMatrixMultiply(const std::vector<std::vector<bool>>& m1, const std::vector<std::vector<bool>>& m2) {
        std::vector<std::vector<bool>> result(m1.size());
        for (std::size_t i = 0; i < m1.size(); i++) {
            result.at(i) = std::vector<bool>(m2.at(0).size());
            for (std::size_t j = 0; j < m2.at(0).size(); j++) {
                result[i][j] = false;

                for (std::size_t k = 0; k < m2.size(); k++) {
                    result[i][j] = result[i][j] ^ (m1[i][k] && m2[k][j]);
                }
            }
        }
        return result;
    }

    static void swapRows(std::vector<std::vector<bool>>& matrix, const std::size_t row1, const std::size_t row2) {
        for (std::size_t col = 0; col < matrix.at(0).size(); col++) {
            std::swap(matrix.at(row1).at(col), matrix.at(row2).at(col));
        }
    }

    static void printGF2matrix(const std::vector<std::vector<bool>>& matrix) {
        for (const auto& i: matrix) {
            for (bool j: i) {
                std::cout << j << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    static void printGF2vector(const std::vector<bool>& vector) {
        if (vector.empty()) {
            std::cout << "[]";
            return;
        }
        for (const auto& i: vector) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
    /**
     * Checks if vec is in the rowspace of matrix M by gaussian elimination (computing reduced echelon form)
     * @param matrix matrix, contains reduced echelon form at end of the function
     * @param vec column vector, not altered by function
     * @return solution to the system of equations, or an empty vector if there is no unique solution

    static std::vector<bool> rowEchelonReduce(std::vector<std::vector<bool>>& matrix, const std::vector<bool>& vect) { //https://stackoverflow.com/questions/11483925/how-to-implementing-gaussian-elimination-for-binary-equations
        auto vec = vect;
        assert(matrix.size() == vec.size());
        std::size_t m = matrix.size();
        std::size_t n = matrix.at(0).size();

        for (std::size_t k = 0, h = 0; h < m && k < n;) {
            long pIdx = -1;
            for (size_t i = h; i < m; i++) {
                if (matrix[i][k]) {
                    pIdx = i;
                    break;
                }
            }
            if (pIdx == -1) {
                k++;
            } else {
                swapRows(matrix, h, pIdx);
                std::swap(vec.at(pIdx), vec.at(h));

                for (size_t i = h + 1; i < m; i++) {
                    matrix[i][k] = false;
                    for (size_t j = k + 1; j < n; j++) {
                        matrix[i][j] = matrix[i][j] - (matrix[h][j] * (matrix[i][k] / matrix[h][k]));
                    }
                }
                h++, k++;
            }
        }
                //alternatively
                for (size_t i = 0; i < m; i++) {
            std::size_t maxi = 0;
            for (size_t k = i; k < n; k++) {
                if (matrix[k][i]) {
                    maxi = k;
                    break;
                }
            }
            if (matrix[maxi][i]) {
                swapRows(matrix, i, maxi);
                for (size_t u = i + 1; u < n; u++) {
                    for (size_t j = 0; j < m; j++) {
                        matrix[u][j] = matrix[u][j] ^ (matrix[i][j] * matrix[u][i]);
                    }
                }
            } else {
                return std::vector<bool>{};
            }
        }

        return vec;
    }
*/
    /**
     * Returns a bitstring representing am n-qubit Pauli error (all Z or all X)
     * The qubits have iid error probabilities given by the parameter
     * @param n
     * @param physicalErrRate
     * @return
     */
    static std::vector<bool>
    sampleErrorIidPauliNoise(const std::size_t n, const double physicalErrRate) {
        std::random_device rd;
        std::mt19937       gen(rd());
        std::vector<bool>  result;

        // Setup the weights, iid noise for each bit
        std::discrete_distribution<> d({1 - physicalErrRate, physicalErrRate});
        for (std::size_t i = 0; i < n; i++) {
            result.emplace_back(d(gen));
        }
        return result;
    }

    /**
         *
         * @param error bool vector representing error
         * @param residual estimate vector that contains residual error at end of function
    */
    static void computeResidualErr(const std::vector<bool>& error, std::vector<bool>& residual) {
        for (std::size_t j = 0; j < residual.size(); j++) {
            residual.at(j) = residual.at(j) ^ error.at(j);
        }
    }
};
#endif //QUNIONFIND_UTILS_HPP
