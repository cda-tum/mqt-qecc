#include "Utils.hpp"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <gf2dense.hpp>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

bool Utils::isVectorInRowspace(const gf2Mat& inmat, const gf2Vec& vec) {
    assertMatrixPresent(inmat);
    assertVectorPresent(vec);
    if (std::none_of(vec.begin(), vec.end(), [](const bool val) { return val; })) { // all zeros vector trivial
        return true;
    }
    if (vec.size() != inmat.at(0).size()) {
        throw QeccException("Cannot check if in rowspace, dimensions of matrix and vector do not match");
    }
    const auto matrix    = getTranspose(inmat); // v is in rowspace of M <=> v is in col space of M^T
    auto       matrixCsc = Utils::toCsc(matrix);

    auto pluDecomp = ldpc::gf2dense::PluDecomposition(static_cast<int>(matrix.size()), static_cast<int>(matrix.at(0).size()), matrixCsc);
    pluDecomp.rref();

    std::vector<int> idxs{};
    for (size_t i = 0; i < vec.size(); i++) {
        if (vec.at(i)) {
            idxs.emplace_back(static_cast<int>(i));
        }
    }
    matrixCsc.emplace_back(idxs);
    auto pluExt = ldpc::gf2dense::PluDecomposition(static_cast<int>(matrix.size()), static_cast<int>(matrix.at(0).size() + 1), matrixCsc);
    pluExt.rref();

    return pluExt.matrix_rank == pluDecomp.matrix_rank;
}

gf2Mat Utils::getTranspose(const gf2Mat& matrix) {
    assertMatrixPresent(matrix);
    gf2Mat transp(matrix.at(0).size(), gf2Vec(matrix.size()));
    for (std::size_t i = 0; i < matrix.size(); i++) {
        const auto& row = matrix.at(i);
        for (std::size_t j = 0; j < row.size(); j++) {
            transp.at(j).at(i) = row.at(j);
        }
    }
    return transp;
}

void Utils::rectMatrixMultiply(const gf2Mat& m1, const gf2Vec& vec, gf2Vec& result) {
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

void Utils::assertMatrixPresent(const gf2Mat& matrix) {
    if (matrix.empty() || matrix.at(0).empty()) {
        throw QeccException("Matrix is empty");
    }
}

void Utils::assertVectorPresent(const gf2Vec& vector) {
    if (vector.empty()) {
        throw QeccException("Vector is empty");
    }
}

void Utils::swapRows(gf2Mat& matrix, const std::size_t row1, const std::size_t row2) {
    for (std::size_t col = 0; col < matrix.at(0).size(); col++) {
        std::vector<bool>::swap(matrix.at(row1).at(col), matrix.at(row2).at(col));
    }
}

void Utils::printGF2matrix(const gf2Mat& matrix) {
    std::cout << getStringFrom(matrix);
}

void Utils::printGF2vector(const gf2Vec& vector) {
    std::cout << getStringFrom(vector);
}

std::string Utils::getStringFrom(const gf2Mat& matrix) {
    if (matrix.empty()) {
        return "[]";
    }
    const auto&       nrows = matrix.size();
    const auto&       ncols = matrix.at(0).size();
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

std::string Utils::getStringFrom(const gf2Vec& vector) {
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

gf2Vec Utils::sampleErrorIidPauliNoise(const std::size_t n, const double physicalErrRate) {
    std::random_device rd;
    std::mt19937_64    gen(rd());
    gf2Vec             result{};
    result.reserve(n);

    // Set up the weights, iid noise for each bit
    std::bernoulli_distribution d(physicalErrRate);
    for (std::size_t i = 0; i < n; i++) {
        result.emplace_back(d(gen));
    }
    return result;
}

void Utils::computeResidualErr(const gf2Vec& error, gf2Vec& residual) {
    for (std::size_t j = 0; j < residual.size(); j++) {
        residual.at(j) = (residual.at(j) != error.at(j));
    }
}

gf2Mat Utils::importGf2MatrixFromFile(const std::string& filepath) {
    std::string   line;
    int           word{};
    std::ifstream inFile(filepath);
    gf2Mat        result{};

    if (!inFile) {
        throw QeccException("Cannot open file");
    }

    while (getline(inFile, line, '\n')) {
        gf2Vec             tempVec{};
        std::istringstream instream(line);
        while (instream >> word) {
            tempVec.emplace_back(static_cast<bool>(word));
        }
        result.emplace_back(tempVec);
    }
    return result;
}

void Utils::printTimePerSampleRun(const std::map<std::string, std::size_t, std::less<>>& avgSampleRuns) {
    std::cout << "trial:timesum = {\n";
    for (const auto& [key, value] : avgSampleRuns) {
        std::cout << "  " << key << ":" << value << ",\n";
    }
    std::cout << "}\n";
}

void Utils::readInFilePathsFromDirectory(const std::string& inPath, std::vector<std::string>& codePaths) {
    for (const auto& file : std::filesystem::directory_iterator(inPath)) {
        codePaths.emplace_back(file.path());
    }
}

std::vector<std::vector<int>> Utils::toCsc(const std::vector<std::vector<bool>>& mat) {
    // convert redHz to int type and to csc format: matrix[col][row] = 1
    if (mat.empty()) {
        return {};
    }
    auto                          rows = mat.size();
    auto                          cols = mat.at(0).size();
    std::vector<std::vector<int>> result;
    for (size_t i = 0; i < cols; i++) {
        std::vector<int> col = {};
        for (size_t j = 0; j < rows; j++) {
            if (mat.at(j).at(i)) {
                col.emplace_back(static_cast<int>(j));
            }
        }
        result.emplace_back(col);
    }
    return result;
}
