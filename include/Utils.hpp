//
// Created by lucas on 09/06/22.
//
#ifndef QUNIONFIND_UTILS_HPP
#define QUNIONFIND_UTILS_HPP

#include "TreeNode.hpp"

#include <cassert>
#include <flint/nmod_matxx.h>
#include <fstream>
#include <iostream>
#include <ostream>
#include <random>
#include <set>
#include <vector>
extern "C" {
#include <flint/nmod_mat.h>
};

typedef std::vector<std::vector<bool>> gf2Mat;
typedef std::vector<bool>              gf2Vec;

class Utils {
public:
    /**
     * Uses flint's integers mod n matrix package nnmod_mat to solve the system given by Mx=b
     * Returns x if there is a solution, or an empty vector if there is no solution
     * By the behaviour of flint's solve function, if there are multiple valid solutions one is returned
     * @param M
     * @param vec
     * @return
     */
    static gf2Vec solveSystem(const gf2Mat& M, const gf2Vec& vec) {
        gf2Vec     result{};
        long       rows = M.size();
        long       cols = M.at(0).size();
        nmod_mat_t mat;
        nmod_mat_t x;
        nmod_mat_t b;
        mp_limb_t  mod = 2U;
        nmod_mat_init(mat, rows, cols, mod);
        nmod_mat_init(x, cols, 1, mod);
        nmod_mat_init(b, rows, 1, mod);

        for (long i = 0; i < nmod_mat_nrows(mat); i++) {
            for (long j = 0; j < nmod_mat_ncols(mat); j++) {
                nmod_mat_set_entry(mat, i, j, M[i][j]);
            }
        }
        auto bColIdx = nmod_mat_ncols(b) - 1;
        for (long i = 0; i < nmod_mat_nrows(b); i++) {
            mp_limb_t tmp = vec[i];
            nmod_mat_set_entry(b, i, bColIdx, tmp);
        }
        int sol = nmod_mat_can_solve(x, mat, b);
        std::cout << "mat: " << std::endl;
        nmod_mat_print_pretty(mat);
        std::cout << "b: " << std::endl;
        nmod_mat_print_pretty(b);

        if (sol == 1) {
            std::cout << "solution exists:" << std::endl;
            nmod_mat_print_pretty(x);
            result       = gf2Vec(nmod_mat_nrows(x));
            auto xColIdx = nmod_mat_ncols(x) - 1;
            for (long i = 0; i < nmod_mat_nrows(x); i++) {
                result.at(i) = nmod_mat_get_entry(x, i, xColIdx);
            }
        } else {
            std::cout << "no sol" << std::endl;
        }

        return result;
    }

    static gf2Mat gauss(const gf2Mat& matrix) {
        gf2Mat result(matrix.at(0).size());
        auto   mat = getFlintMatrix(matrix);
        mat.set_rref(); // reduced row echelon form
        return getMatrixFromFlint(mat);
    }

    static flint::nmod_matxx getFlintMatrix(const gf2Mat& matrix) {
        auto ctxx   = flint::nmodxx_ctx(2);
        auto result = flint::nmod_matxx(matrix.size(), matrix.at(0).size(), 2);
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix.at(0).size(); j++) {
                if (matrix[i][j]) {
                    result.at(i, j) = flint::nmodxx::red(1, ctxx);
                } else {
                    result.at(i, j) = flint::nmodxx::red(0, ctxx);
                }
            }
        }
        return result;
    }

    static gf2Mat getMatrixFromFlint(const flint::nmod_matxx& matrix) {
        auto   ctxx = flint::nmodxx_ctx(2);
        gf2Mat result(matrix.rows());
        auto   a = flint::nmodxx::red(1, ctxx);

        for (long i = 0; i < matrix.rows(); i++) {
            result.at(i) = gf2Vec(matrix.cols());
            for (long j = 0; j < matrix.cols(); j++) {
                if (matrix.at(i, j) == a) {
                    result[i][j] = true;
                } else {
                    result[i][j] = false;
                }
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
    static bool isVectorInRowspace(const gf2Mat& M, const gf2Vec& vec) {
        if (std::none_of(vec.begin(), vec.end(), [](const bool val) { return val; })) { // all zeros vector trivial
            return true;
        }
        gf2Mat matrix;
        if (vec.size() == M.at(0).size()) {
            matrix = getTranspose(M); // v is in rowspace of M <=> v is in col space of M^T
        } else {
            throw std::errc::invalid_argument;
        }
        auto augm = getAugmentedMatrix(matrix, vec);
        matrix    = gauss(augm);
        printGF2matrix(matrix);
        gf2Vec vector(vec.size());

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

    /**
     * Computes and returns the matrix obtained by appending the column vector to the input matrix, result = (matrix|vector)
     * @param matrix
     * @param vector
     * @return
     */
    static gf2Mat getAugmentedMatrix(const gf2Mat& matrix, const gf2Vec& vector) {
        gf2Mat result(matrix.size());

        for (size_t i = 0; i < matrix.size(); i++) {
            result.at(i) = gf2Vec(matrix.at(i).size() + 1);
            for (std::size_t j = 0; j < matrix.at(0).size(); j++) {
                result[i][j] = matrix[i][j];
            }
            result[i][matrix.at(0).size()] = vector.at(i);
        }
        return result;
    }

    /**
     * Computes the transpose of the given matrix
     * @param matrix
     * @return
     */
    static gf2Mat getTranspose(const gf2Mat& matrix) {
        gf2Mat transp(matrix.at(0).size());
        for (auto& i: transp) {
            i = gf2Vec(matrix.size());
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
    static gf2Mat rectMatrixMultiply(const gf2Mat& m1, const gf2Mat& m2) {
        auto mat1   = getFlintMatrix(m1);
        auto mat2   = getFlintMatrix(m2);
        auto result = flint::nmod_matxx(mat1.rows(), mat2.cols(), 2);
        result      = mat1.mul_classical(mat2);
        return getMatrixFromFlint(result);
    }

    static void swapRows(gf2Mat& matrix, const std::size_t row1, const std::size_t row2) {
        for (std::size_t col = 0; col < matrix.at(0).size(); col++) {
            std::swap(matrix.at(row1).at(col), matrix.at(row2).at(col));
        }
    }

    static void printGF2matrix(const gf2Mat& matrix) {
        std::cout << getStringFrom(matrix);
    }

    static void printGF2vector(const gf2Vec& vector) {
        std::cout << getStringFrom(vector);
    }

    static std::string getStringFrom(const gf2Mat& matrix) {
        if (matrix.empty()) {
            return "[]";
        }
        std::stringstream s;
        s << matrix.size() << "x" << matrix.at(0).size() << "matrix [" << std::endl;
        for (const auto& i: matrix) {
            s << "[ ";
            for (bool j: i) {
                s << j << " ";
            }
            s << "]";
            s << std::endl;
        }
        s << "]";
        return s.str();
    }

    static std::string getStringFrom(const gf2Vec& vector) {
        if (vector.empty()) {
            return "[]";
        }
        std::stringstream s;
        s << "[ ";
        for (const auto& i: vector) {
            s << i << " ";
        }
        s << "]";
        return s.str();
    }

    /**
     * Returns a bitstring representing am n-qubit Pauli error (all Z or all X)
     * The qubits have iid error probabilities given by the parameter
     * @param n
     * @param physicalErrRate
     * @return
     */
    static gf2Vec sampleErrorIidPauliNoise(const std::size_t n, const double physicalErrRate) {
        std::random_device rd;
        std::mt19937       gen(rd());
        gf2Vec             result;

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
    static void computeResidualErr(const gf2Vec& error, gf2Vec& residual) {
        for (std::size_t j = 0; j < residual.size(); j++) {
            residual.at(j) = residual.at(j) ^ error.at(j);
        }
    }

    static gf2Mat importGf2MatrixFromFile(const std::string& filepath) {
        std::string   line;
        int           word;
        std::ifstream inFile(filepath);
        gf2Mat        result;

        if (inFile) {
            while (getline(inFile, line, '\n')) {
                gf2Vec             tempVec;
                std::istringstream instream(line);
                while (instream >> word) {
                    tempVec.push_back(word);
                }
                result.emplace_back(tempVec);
            }
        } else {
            std::cerr << "File " << filepath << " cannot be opened." << std::endl;
        }

        inFile.close();
        std::cout << "imported: " << std::endl;
        printGF2matrix(result);
        std::cout << std::endl;
        return result;
    }
};
#endif //QUNIONFIND_UTILS_HPP
