//
// Created by lucas on 09/06/22.
//
#ifndef QUNIONFIND_UTILS_HPP
#define QUNIONFIND_UTILS_HPP

#include "QeccException.hpp"
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

using gf2Mat = std::vector<std::vector<bool>>;
using gf2Vec = std::vector<bool>;

class Utils {
public:
    // NOLINTBEGIN(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    /**
     * Uses flint's integers mod n matrix package nnmod_mat to solve the system given by Mx=b
     * Returns x if there is a solution, or an empty vector if there is no solution
     * By the behaviour of flint's solve function, if there are multiple valid solutions one is returned
     * @param inmat
     * @param vec
     * @return
     */
    static gf2Vec solveSystem(const gf2Mat& inmat, const gf2Vec& vec) {
        assertMatrixPresent(inmat);
        assertVectorPresent(vec);
        if (inmat.size() > std::numeric_limits<long int>::max() || inmat.front().size() > std::numeric_limits<long int>::max()) { // NOLINT(google-runtime-int)
            throw QeccException("size of matrix too large for flint");
        }
        if (inmat.size() != vec.size()) { // NOLINT(readability-else-after-return)
            std::cerr << "Cannot solve system, dimensions do not match" << std::endl;
            throw QeccException("Cannot solve system, dimensions do not match");
        }

        gf2Vec          result{};
        slong           rows = static_cast<std::int64_t>(inmat.size());
        slong           cols = static_cast<std::int64_t>(inmat.front().size());
        nmod_mat_t      mat;
        nmod_mat_t      x;
        nmod_mat_t      b;
        mp_limb_t const mod = 2;
        // initializes mat to rows x cols matrix with coefficients mod 2
        nmod_mat_init(mat, rows, cols, mod); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
        nmod_mat_init(x, cols, 1, mod);
        nmod_mat_init(b, rows, 1, mod);

        for (slong i = 0; i < nmod_mat_nrows(mat); i++) {
            for (slong j = 0; j < nmod_mat_ncols(mat); j++) {
                nmod_mat_set_entry(mat, i, j, inmat.at(static_cast<std::uint64_t>(i)).at(static_cast<std::uint64_t>(j)) ? 1U : 0U);
            }
        }
        slong bColIdx = nmod_mat_ncols(b) - 1;
        for (slong i = 0; i < nmod_mat_nrows(b); i++) {
            nmod_mat_set_entry(b, i, bColIdx, vec.at(static_cast<std::uint64_t>(i)) ? 1U : 0U);
        }
        int const sol = nmod_mat_can_solve(x, mat, b);

        if (sol == 1) {
            result       = gf2Vec(static_cast<std::uint64_t>(nmod_mat_nrows(x)));
            auto xColIdx = nmod_mat_ncols(x) - 1;
            for (auto i = 0; i < nmod_mat_nrows(x); i++) {
                result.at(static_cast<std::size_t>(i)) = nmod_mat_get_entry(x, i, xColIdx); // NOLINT(readability-implicit-bool-conversion)
            }
        } else {
            // no solution
        }
        nmod_mat_clear(mat);
        nmod_mat_clear(x);
        nmod_mat_clear(b);
        return result;
    }
    // NOLINTEND(cppcoreguidelines-pro-bounds-array-to-pointer-decay)

    static flint::nmod_matxx gauss(const gf2Mat& matrix) {
        assertMatrixPresent(matrix);
        auto res = getFlintMatrix(matrix);
        res.set_rref(); // reduced row echelon form
        return res;
    }

    static flint::nmod_matxx getFlintMatrix(const gf2Mat& matrix) {
        assertMatrixPresent(matrix);
        const slong     rows    = static_cast<slong>(matrix.size());
        const slong     cols    = static_cast<slong>(matrix.front().size());
        const mp_limb_t modulus = 2;
        const auto&     ctxx    = flint::nmodxx_ctx(modulus);
        auto            result  = flint::nmod_matxx(matrix.size(), matrix.front().size(), modulus);
        for (slong i = 0; i < rows; i++) {
            for (slong j = 0; j < cols; j++) {
                if (matrix.at(static_cast<std::uint64_t>(i)).at(static_cast<std::uint64_t>(j))) {
                    const mp_limb_t one = 1U;
                    result.at(i, j)     = flint::nmodxx::red(one, ctxx);
                } else {
                    const mp_limb_t zero = 0;
                    result.at(i, j)      = flint::nmodxx::red(zero, ctxx);
                }
            }
        }
        return result;
    }

    static gf2Mat getMatrixFromFlint(const flint::nmod_matxx& matrix) {
        const auto& ctxx = flint::nmodxx_ctx(2);
        gf2Mat      result(static_cast<std::uint64_t>(matrix.rows()));
        const auto& a = flint::nmodxx::red(1, ctxx);

        for (slong i = 0; i < matrix.rows(); i++) {
            result.at(static_cast<std::uint64_t>(i)) = gf2Vec(static_cast<std::uint64_t>(matrix.cols()));
            for (slong j = 0; j < matrix.cols(); j++) {
                if (matrix.at(i, j) == a) {
                    result.at(static_cast<std::uint64_t>(i)).at(static_cast<std::uint64_t>(j)) = true;
                } else {
                    result.at(static_cast<std::uint64_t>(i)).at(static_cast<std::uint64_t>(j)) = false;
                }
            }
        }
        return result;
    }

    /**
     * Checks if the given vector is in the rowspace of matrix M
     * @param inmat
     * @param vec
     * @return
     */
    static bool isVectorInRowspace(const gf2Mat& inmat, const gf2Vec& vec) {
        assertMatrixPresent(inmat);
        assertVectorPresent(vec);
        if (std::none_of(vec.begin(), vec.end(), [](const bool val) { return val; })) { // all zeros vector trivial
            return true;
        }
        gf2Mat matrix;
        if (vec.size() == inmat.at(0).size()) {
            matrix = getTranspose(inmat); // v is in rowspace of M <=> v is in col space of M^T
        } else {
            throw QeccException("Cannot check if in rowspace, dimensions of matrix and vector do not match");
        }

        for (std::size_t i = 0; i < matrix.size(); i++) {
            matrix.at(i).emplace_back(vec.at(i));
        }
        auto reduced = gauss(matrix);
        // flint::print_pretty(reduced);
        //  check consistency, inconsistent <=> vec not in rowspace
        for (slong i = 0; i < reduced.rows(); i++) {
            if (reduced.at(i, reduced.cols() - 1)._limb() == 1) {
                bool inconsistent = true;
                for (slong k = 0; k < (reduced.cols() - 1); k++) {
                    if (reduced.at(i, k)._limb() == 1) {
                        inconsistent = false;
                    }
                }
                if (inconsistent) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Computes the transpose of the given matrix
     * @param matrix
     * @return
     */
    static gf2Mat
    getTranspose(const gf2Mat& matrix) {
        assertMatrixPresent(matrix);
        gf2Mat transp(matrix.at(0).size());
        for (auto& i : transp) {
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
     * Computes matrix vector product and stores it in result vector
     * @param m1
     * @param vec
     * @param result
     */
    static void rectMatrixMultiply(const gf2Mat& m1, const gf2Vec& vec, gf2Vec& result) {
        assertMatrixPresent(m1);
        assertVectorPresent(vec);
        if (m1.front().size() != vec.size() || m1.size() > result.capacity()) {
            throw QeccException("Cannot multiply, dimensions wrong");
        }
        for (std::size_t i = 0; i < m1.size(); i++) {
            const auto& row = m1.at(i);
            for (std::size_t k = 0; k < vec.size(); k++) {
                result.at(i) = result.at(i) ^ (row.at(k) && vec.at(k));
            }
        }
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
            std::vector<bool>::swap(matrix.at(row1).at(col), matrix.at(row2).at(col));
        }
    }

    [[maybe_unused]] static void printGF2matrix(const gf2Mat& matrix) {
        std::cout << getStringFrom(matrix);
    }

    static void printGF2vector(const gf2Vec& vector) {
        std::cout << getStringFrom(vector);
    }

    static std::string getStringFrom(const gf2Mat& matrix) {
        if (matrix.empty()) {
            return "[]";
        }
        const auto&       nrows = matrix.size();
        const auto&       ncols = matrix.at(0).size();
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
        const auto&       nelems = vector.size();
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
        std::mt19937_64    gen(rd());
        gf2Vec             result;
        result.reserve(n);

        // Set up the weights, iid noise for each bit
        std::bernoulli_distribution d(physicalErrRate);
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
        int           word; // NOLINT(cppcoreguidelines-init-variables)
        std::ifstream inFile(filepath);
        gf2Mat        result;

        if (inFile) {
            while (getline(inFile, line, '\n')) {
                gf2Vec             tempVec;
                std::istringstream instream(line);
                while (instream >> word) {
                    tempVec.push_back(static_cast<bool>(word));
                }
                result.emplace_back(tempVec);
            }
        } else {
            std::cerr << "File " << filepath << " cannot be opened." << std::endl;
        }

        inFile.close();
        return result;
    }
    [[maybe_unused]] static void printTimePerSampleRun(const std::map<std::string, std::size_t, std::less<>>& avgSampleRuns) {
        const nlohmann::json avgData = avgSampleRuns;
        std::cout << "trial:timesum = " << avgData.dump(2U) << std::endl;
    }

    [[maybe_unused]] static void readInFilePathsFromDirectory(const std::string& inPath, std::vector<std::string>& codePaths) {
        for (const auto& file : std::filesystem::directory_iterator(inPath)) {
            codePaths.emplace_back(file.path());
        }
    }
};
#endif // QUNIONFIND_UTILS_HPP
