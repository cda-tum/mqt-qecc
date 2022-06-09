//
// Created by luca on 09/06/22.
//
#include "TreeNode.hpp"

#include <iostream>
#include <ostream>
#include <set>
#include <vector>

#ifndef QUNIONFIND_UTILS_HPP
#define QUNIONFIND_UTILS_HPP
class Utils {
public:
    static std::vector<std::vector<bool>> rectMatrixMultiply(const std::vector<std::vector<bool>>& m1, const std::vector<std::vector<bool>>& m2) {
        std::vector<std::vector<bool>> result(m1.size());
        for (std::size_t i = 0; i < m1.size(); i++) {
            result.at(i) = std::vector<bool>(m2.at(0).size());
            for (std::size_t j = 0; j < m2.at(0).size(); j++) {
                result[i][j] = false;

                for (std::size_t k = 0; k < m2.size(); k++) {
                    result[i][j] = result[i][j] + (m1[i][k] * m2[k][j]);
                }
            }
        }
        return result;
    }

    static void swapRows(std::vector<std::vector<bool>>& matrix, const std::size_t row1, const std::size_t row2) {
        for (std::size_t col = 0; col <= matrix.at(0).size(); col++) {
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

    /*
     * Checks if vec is in the rowspace of matrix M by gaussian elimination (computing reduced echelon form)
     */
    static bool checkVectorInRowspace(std::vector<std::vector<bool>> M, std::vector<bool> vec) { //https://stackoverflow.com/questions/11483925/how-to-implementing-gaussian-elimination-for-binary-equations
        std::size_t                    nrCols = M.size();
        std::size_t                    nrRows = M.at(0).size();
        std::size_t                    row    = 0;
        std::vector<std::vector<bool>> transp(nrRows);
        for (auto& i: transp) {
            i = std::vector<bool>(nrCols);
        }
        for (size_t i = 0; i < M.size(); i++) {
            for (size_t j = 0; j < M.at(i).size(); j++) {
                transp[j][i] = M[i][j];
            }
        }
        M = transp;
        for (std::size_t col = 0; col < nrCols && row < nrRows; col++, row++) {
            if (M[row][col] == 0) {
                for (std::size_t i = 0; i < nrRows; i++) {
                    if (M[i][col] != 0) {
                        Utils::swapRows(M, i, row);
                        std::swap(vec.at(i), vec.at(row));
                    }
                }
            }
            if (M[row][col] == 0) {
                return false;
            }
            for (std::size_t j = 0; j < nrRows; ++j) {
                if (j != col) {
                    if (M[j][col]) {
                        std::size_t k;
                        for (k = col; k < nrCols; ++k) {
                            M[j][k] = M[j][k] ^ M[col][k];
                        }
                        vec[k] = vec[k] ^ vec[col];
                    }
                }
            }
        }
        for (size_t i = 0; i < vec.size(); i++) {
            if (vec[i]) {
                for (size_t j = 0; j < M.at(i).size(); j++) {
                    if (M[i][j]) {
                        return false;
                    }
                }
            }
        }
        return true;
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
};
#endif //QUNIONFIND_UTILS_HPP
