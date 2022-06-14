//
// Created by luca on 09/06/22.
//
#ifndef QUNIONFIND_UTILS_HPP
#define QUNIONFIND_UTILS_HPP
#include "TreeNode.hpp"

#include <cassert>
#include <iostream>
#include <ostream>
#include <random>
#include <set>
#include <vector>

class Utils {
public:
    static std::vector<bool> solveSystem(const std::vector<std::vector<bool>>& M, const std::vector<bool>& vec) { //TODO fix
        auto              matrix = M;
        std::vector<bool> v      = {1, 1, 1, 1, 1, 1, 1};
        auto              sol    = Utils::rowEchelonReduce(matrix, vec);
        std::vector<bool> result(matrix.size());

        // back substitution
        for (std::size_t i = result.size() - 1; i > (-1); i--) {
            //result[i] = matrix[i][i] && sol[i];
            for (std::size_t j = (i - 1); j > (-1); j--) {
                matrix[i][i] = matrix[i][i] ^ (matrix[j][i] && result[i]);
            }
        }
        return result;
    }

    static bool isVectorInRowspace(const std::vector<std::vector<bool>>& M, const std::vector<bool>& sol) {
        auto matrix = M;
        auto vector = Utils::rowEchelonReduce(matrix, sol);
        // check consistency
        if (vector.empty()) {
            return false;
        }
        for (size_t i = 0; i < vector.size(); i++) {
            if (vector[i]) {
                for (size_t j = 0; j < matrix.at(i).size(); j++) {
                    if (std::none_of(matrix.at(i).begin(), matrix.at(i).end(), [](const bool val) { return val; })) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

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
     */
    static std::vector<bool> rowEchelonReduce(std::vector<std::vector<bool>>& matrix, const std::vector<bool>& vect) { //https://stackoverflow.com/questions/11483925/how-to-implementing-gaussian-elimination-for-binary-equations
        auto vec = vect;
        assert(matrix.size() == vec.size());
        std::size_t nrRows = matrix.size();
        std::size_t nrCols = matrix.at(0).size();
        std::size_t row    = 0;

        for (std::size_t col = 0; col < nrCols && row < nrRows; col++, row++) {
            if (matrix[col][col] == 0) {
                for (std::size_t i = 0; i < nrRows; i++) {
                    if (matrix[i][col] != 0) {
                        std::cout << "swapping " << i << "and" << row << std::endl;
                        Utils::swapRows(matrix, row, i);
                        std::swap(vec.at(i), vec.at(row));
                    }
                }
            }
            if (matrix[col][col] == 0) {
                return std::vector<bool>{};
            }
            for (std::size_t j = 0; j < nrRows; j++) {
                if (j != col) {
                    if (matrix[j][col]) {
                        std::size_t k;
                        for (k = col; k < nrCols; k++) {
                            matrix[j][k] = matrix[j][k] ^ matrix[row][k];
                        }
                        vec[j] = vec[j] ^ vec[row];
                    }
                }
            }
        }
        return vec;
    }

    static std::vector<bool> sampleErrorIidPauliNoise(const std::size_t n, const double physicalErrRate) {
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
