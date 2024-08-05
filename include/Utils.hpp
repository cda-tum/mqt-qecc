#pragma once

#include "QeccException.hpp"
#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <gf2dense.hpp>
#include <iostream>
#include <ostream>
#include <random>
#include <vector>

using gf2Mat = std::vector<std::vector<bool>>;
using gf2Vec = std::vector<bool>;

class Utils {
public:
    /**
     * Checks if the given vector is in the rowspace of matrix M
     * @param inmat
     * @param vec
     * @return
     */
    static bool isVectorInRowspace(const gf2Mat &inmat, const gf2Vec &vec) {
        assertMatrixPresent(inmat);
        assertVectorPresent(vec);
        if (std::none_of(vec.begin(), vec.end(), [](const bool val) { return val; })) { // all zeros vector trivial
            return true;
        }
        gf2Mat matrix = {};
        if (vec.size() == inmat.at(0).size()) {
            matrix = getTranspose(inmat); // v is in rowspace of M <=> v is in col space of M^T
        } else {
            throw QeccException("Cannot check if in rowspace, dimensions of matrix and vector do not match");
        }
        std::vector<std::vector<int>> matrixCsc = {};
        matrixCsc = Utils::toCsc(matrix);
        std::vector<int> idxs = {};
        for (auto i = 0; i < vec.size(); i++) {
            if (vec.at(i)) {
                idxs.push_back(i);
            }
        }
        auto pluDecomp = ldpc::gf2dense::PluDecomposition(matrix.size(), matrix.at(0).size(), matrixCsc);
        pluDecomp.rref();

        matrixCsc.push_back(idxs);
        auto pluExt = ldpc::gf2dense::PluDecomposition(matrix.size(), matrix.at(0).size() + 1, matrixCsc);
        pluExt.rref();

        return pluExt.matrix_rank == pluDecomp.matrix_rank;
    }

    /**
     * Computes the transpose of the given matrix
     * @param matrix
     * @return
     */
    static gf2Mat
    getTranspose(const gf2Mat &matrix) {
        assertMatrixPresent(matrix);
        gf2Mat transp(matrix.at(0).size());
        for (auto &i: transp) {
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
    static void rectMatrixMultiply(const gf2Mat &m1, const gf2Vec &vec, gf2Vec &result) {
        assertMatrixPresent(m1);
        assertVectorPresent(vec);
        if (m1.front().size() != vec.size() || m1.size() > result.capacity()) {
            throw QeccException("Cannot multiply, dimensions wrong");
        }
        for (std::size_t i = 0; i < m1.size(); i++) {
            const auto &row = m1.at(i);
            for (std::size_t k = 0; k < vec.size(); k++) {
                result.at(i) = result.at(i) ^ (row.at(k) && vec.at(k));
            }
        }
    }

    static void assertMatrixPresent(const gf2Mat &matrix) {
        if (matrix.empty() || matrix.at(0).empty()) {
            throw QeccException("Matrix is empty");
        }
    }

    static void assertVectorPresent(const gf2Vec &vector) {
        if (vector.empty()) {
            throw QeccException("Vector is empty");
        }
    }

    static void swapRows(gf2Mat &matrix, const std::size_t row1, const std::size_t row2) {
        for (std::size_t col = 0; col < matrix.at(0).size(); col++) {
            std::vector<bool>::swap(matrix.at(row1).at(col), matrix.at(row2).at(col));
        }
    }

    [[maybe_unused]] static void printGF2matrix(const gf2Mat &matrix) {
        std::cout << getStringFrom(matrix);
    }

    static void printGF2vector(const gf2Vec &vector) {
        std::cout << getStringFrom(vector);
    }

    static std::string getStringFrom(const gf2Mat &matrix) {
        if (matrix.empty()) {
            return "[]";
        }
        const auto &nrows = matrix.size();
        const auto &ncols = matrix.at(0).size();
        std::stringstream s;
        s << nrows << "x" << ncols << "matrix [\n";
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
            s << '\n';
        }
        s << "]";
        return s.str();
    }

    static std::string getStringFrom(const gf2Vec &vector) {
        if (vector.empty()) {
            return "[]";
        }
        const auto &nelems = vector.size();
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
        std::mt19937_64 gen(rd()); // NOLINT(cppcoreguidelines-init-variables)
        gf2Vec result = {};
        result.reserve(n);

        // Set up the weights, iid noise for each bit
        std::bernoulli_distribution d(physicalErrRate); // NOLINT(cppcoreguidelines-init-variables)
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
    static void computeResidualErr(const gf2Vec &error, gf2Vec &residual) {
        for (std::size_t j = 0; j < residual.size(); j++) {
            residual.at(j) = (residual.at(j) != error.at(j));
        }
    }

    static gf2Mat importGf2MatrixFromFile(const std::string &filepath) {
        std::string line;
        int word; // NOLINT(cppcoreguidelines-init-variables)
        std::ifstream inFile(filepath);
        gf2Mat result = {};

        if (!inFile) {
            throw QeccException("Cannot open file");
        }

        while (getline(inFile, line, '\n')) {
            gf2Vec tempVec = {};
            std::istringstream instream(line);
            while (instream >> word) {
                tempVec.push_back(static_cast<bool>(word));
            }
            result.emplace_back(tempVec);
        }
        return result;
    }

    [[maybe_unused]] static void
    printTimePerSampleRun(const std::map<std::string, std::size_t, std::less<>> &avgSampleRuns) {
        std::cout << "trial:timesum = {\n";
        for (const auto &[key, value]: avgSampleRuns) {
            std::cout << "  " << key << ":" << value << ",\n";
        }
        std::cout << "}\n";
    }

    [[maybe_unused]] static void
    readInFilePathsFromDirectory(const std::string &inPath, std::vector<std::string> &codePaths) {
        for (const auto &file: std::filesystem::directory_iterator(inPath)) {
            codePaths.emplace_back(file.path());
        }
    }

    static std::vector<std::vector<int>> toCsc(const std::vector<std::vector<bool>> &mat) {
        // convert redHz to int type and to csc format: matrix[col][row] = 1
        if (mat.empty()) {
            return std::vector<std::vector<int>>();
        }
        auto rows = mat.size();
        auto cols = mat.at(0).size();
        std::vector<std::vector<int>> result;
        for (auto i = 0; i < cols; i++) {
            std::vector<int> col = {};
            for (auto j = 0; j < rows; j++) {
                if (mat.at(j).at(i)) {
                    col.push_back(j);
                }
            }
            result.push_back(col);
        }
        return result;
    }
};
