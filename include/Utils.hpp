//
// Created by lucas on 09/06/22.
//
#ifndef QUNIONFIND_UTILS_HPP
#define QUNIONFIND_UTILS_HPP

#include "TreeNode.hpp"
#include "nlohmann/json.hpp"

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
}

using gf2Mat = std::vector<std::vector<bool>>;
using gf2Vec = std::vector<bool>;

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
        assertMatrixPresent(M);
        assertVectorPresent(vec);
        if (M.size() > std::numeric_limits<long int>::max() || M.front().size() > std::numeric_limits<long int>::max()) {
            throw QeccException("size of matrix too large for flint");
        }

        gf2Vec     result{};
        slong      rows = M.size();
        slong      cols = M.front().size();
        nmod_mat_t mat;
        nmod_mat_t x;
        nmod_mat_t b;
        mp_limb_t  mod = 2;
        // initializes mat to rows x cols matrix with coefficients mod 2
        nmod_mat_init(mat, rows, cols, mod);
        nmod_mat_init(x, cols, 1, mod);
        nmod_mat_init(b, rows, 1, mod);

        for (slong i = 0; i < nmod_mat_nrows(mat); i++) {
            for (slong j = 0; j < nmod_mat_ncols(mat); j++) {
                nmod_mat_set_entry(mat, i, j, M.at(i).at(j) ? 1U : 0U);
            }
        }
        slong bColIdx = nmod_mat_ncols(b) - 1;
        for (slong i = 0; i < nmod_mat_nrows(b); i++) {
            nmod_mat_set_entry(b, i, bColIdx, vec.at(i) ? 1U : 0U);
        }
        int sol = nmod_mat_can_solve(x, mat, b);
        //std::cout << "mat: " << std::endl;
        nmod_mat_print_pretty(mat);
        //std::cout << "b: " << std::endl;
        nmod_mat_print_pretty(b);

        if (sol == 1) {
            //std::cout << "solution exists:" << std::endl;
            nmod_mat_print_pretty(x);
            result       = gf2Vec(nmod_mat_nrows(x));
            auto xColIdx = nmod_mat_ncols(x) - 1;
            for (long i = 0; i < nmod_mat_nrows(x); i++) {
                result.at(i) = nmod_mat_get_entry(x, i, xColIdx);
            }
        } else {
            //std::cout << "no sol" << std::endl;
        }
        nmod_mat_clear(mat);
        nmod_mat_clear(x);
        nmod_mat_clear(b);
        return result;
    }

    static gf2Mat gauss(const gf2Mat& matrix) {
        assertMatrixPresent(matrix);
        gf2Mat result(matrix.at(0).size());
        auto   mat = getFlintMatrix(matrix);
        mat.set_rref(); // reduced row echelon form
        return getMatrixFromFlint(mat);
    }

    static flint::nmod_matxx getFlintMatrix(const gf2Mat& matrix) {
        assertMatrixPresent(matrix);
        const slong      rows = static_cast<slong>(matrix.size());
        const slong      cols = static_cast<slong>(matrix.front().size());
        const mp_limb_t modul = 2;
        const auto& ctxx   = flint::nmodxx_ctx(modul);
        auto result = flint::nmod_matxx(matrix.size(), matrix.front().size(), modul);
        for (slong i = 0; i < rows; i++) {
            for (slong j = 0; j < cols; j++) {
                if (matrix.at(i).at(j)) {
                    const mp_limb_t one = 1U;
                    result.at(i, j) = flint::nmodxx::red(one, ctxx);
                } else {
                    const mp_limb_t zero = 0;
                    result.at(i, j) = flint::nmodxx::red(zero, ctxx);
                }
            }
        }
        return result;
    }

    static gf2Mat getMatrixFromFlint(const flint::nmod_matxx& matrix) {
        const auto&   ctxx = flint::nmodxx_ctx(2);
        gf2Mat result(matrix.rows());
        const auto&   a = flint::nmodxx::red(1, ctxx);

        for (slong i = 0; i < matrix.rows(); i++) {
            result.at(i) = gf2Vec(matrix.cols());
            for (slong j = 0; j < matrix.cols(); j++) {
                if (matrix.at(i, j) == a) {
                    result.at(i).at(j) = true;
                } else {
                    result.at(i).at(j) = false;
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
        assertMatrixPresent(M);
        assertVectorPresent(vec);
        if (std::none_of(vec.begin(), vec.end(), [](const bool val) { return val; })) { // all zeros vector trivial
            return true;
        }
        gf2Mat matrix;
        if (vec.size() == M.at(0).size()) {
            matrix = getTranspose(M); // v is in rowspace of M <=> v is in col space of M^T
        } else {
            throw QeccException("Cannot check if in rowspace, dimensions of matrix and vector do not match");
        }
        const auto& augm = getAugmentedMatrix(matrix, vec);
        matrix    = gauss(augm);
        gf2Vec vector(vec.size());

        for (std::size_t i = 0; i < matrix.size(); i++) {
            vector.at(i) = matrix[i][matrix.at(i).size() - 1];
        }
        // check consistency
        for (std::size_t i = 0; i < vector.size(); i++) {
            if (vector[i]) {
                for (std::size_t j = 0; j < matrix.at(i).size(); j++) {
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
        assertMatrixPresent(matrix);
        assertVectorPresent(vector);
        gf2Mat result(matrix.size());

        for (std::size_t i = 0; i < matrix.size(); i++) {
            result.at(i) = gf2Vec(matrix.at(i).size() + 1);
            for (std::size_t j = 0; j < matrix.at(0).size(); j++) {
                result.at(i).at(j) = matrix.at(i).at(j);
            }
            result.at(i).at(matrix.at(0).size()) = vector.at(i);
        }
        return result;
    }

    /**
     * Computes the transpose of the given matrix
     * @param matrix
     * @return
     */
    static gf2Mat getTranspose(const gf2Mat& matrix) {
        assertMatrixPresent(matrix);
        gf2Mat transp(matrix.at(0).size());
        for (auto& i: transp) {
            i = gf2Vec(matrix.size());
        }
        for (std::size_t i = 0; i < matrix.size(); i++) {
            for (std::size_t j = 0; j < matrix.at(i).size(); j++) {
                transp.at(j).at(i) = matrix.at(i).at(j);
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
        assertMatrixPresent(m1);
        assertMatrixPresent(m2);
        const auto& mat1   = getFlintMatrix(m1);
        const auto& mat2   = getFlintMatrix(m2);
        auto result = flint::nmod_matxx(mat1.rows(), mat2.cols(), 2);
        result      = mat1.mul_classical(mat2);
        return getMatrixFromFlint(result);
    }

    static void assertMatrixPresent(const gf2Mat& matrix) {
        if (matrix.empty() || matrix.at(0).empty()) {
            throw QeccException("Matrix is empty");
        }
    }

    static void assertVectorPresent(const gf2Vec& vector) {
        if (vector.empty()) {
            throw QeccException("Vector is empty");
        }
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
        const auto&              nrows = matrix.size();
        const auto&              ncols = matrix.at(0).size();
        std::stringstream s;
        s << nrows << "x" << ncols << "matrix [" << std::endl;
        for (std::size_t i = 0; i < nrows; i++) {
            s << "[";
            for (std::size_t j = 0; j < ncols; j++) {
                s << matrix.at(i).at(j);
                if (j != ncols - 1) {
                    s << ",";
                }
            }
            s << "]";
            if (i != nrows - 1) {
                s << ",";
            }
            s << std::endl;
        }
        s << "]";
        return s.str();
    }

    static std::string getStringFrom(const gf2Vec& vector) {
        if (vector.empty()) {
            return "[]";
        }
        const auto&              nelems = vector.size();
        std::stringstream s;
        s << "[";
        for (std::size_t j = 0; j < nelems; j++) {
            s << vector.at(j);
            if (j != nelems - 1) {
                s << ",";
            }
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

        // Set up the weights, iid noise for each bit
        std::discrete_distribution d({1 - physicalErrRate, physicalErrRate});
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
            residual.at(j) = (residual.at(j) != error.at(j));
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
        return result;
    }
    static void printTimePerSampleRun(const std::map<std::string, std::size_t, std::less<>>& avgSampleRuns){
        nlohmann::json avgData = avgSampleRuns;
        std::cout << "trial:timesum = " << avgData.dump(2U) << std::endl;
    }

    static void readInFilePathsFromDirectory(const std::string& inPath, std::vector<std::string>& codePaths){
        for (const auto& file: std::filesystem::directory_iterator(inPath)) {
            codePaths.emplace_back(file.path());
        }
    }
};
#endif //QUNIONFIND_UTILS_HPP
