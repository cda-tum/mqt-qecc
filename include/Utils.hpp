#pragma once

#include "GF2.hpp"

#include <cassert>
#include <cstddef>
#include <functional>
#include <map>
#include <string>
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
    static bool isVectorInRowspace(const gf2Mat& inmat, const gf2Vec& vec);

    /**
     * Computes the transpose of the given matrix
     * @param matrix
     * @return
     */
    static gf2Mat
    getTranspose(const gf2Mat& matrix);

    /**
     * Computes matrix vector product and stores it in result vector
     * @param m1
     * @param vec
     * @param result
     */
    static void rectMatrixMultiply(const gf2Mat& m1, const gf2Vec& vec, gf2Vec& result);

    static void assertMatrixPresent(const gf2Mat& matrix);

    static void assertVectorPresent(const gf2Vec& vector);

    static void swapRows(gf2Mat& matrix, std::size_t row1, std::size_t row2);

    [[maybe_unused]] static void printGF2matrix(const gf2Mat& matrix);

    static void printGF2vector(const gf2Vec& vector);

    static std::string getStringFrom(const gf2Mat& matrix);

    static std::string getStringFrom(const gf2Vec& vector);

    /**
     * Returns a bitstring representing am n-qubit Pauli error (all Z or all X)
     * The qubits have iid error probabilities given by the parameter
     * @param n
     * @param physicalErrRate
     * @return
     */
    static gf2Vec sampleErrorIidPauliNoise(std::size_t n, double physicalErrRate);

    /**
     *
     * @param error bool vector representing error
     * @param residual estimate vector that contains residual error at end of function
     */
    static void computeResidualErr(const gf2Vec& error, gf2Vec& residual);

    static gf2Mat importGf2MatrixFromFile(const std::string& filepath);

    [[maybe_unused]] static void
    printTimePerSampleRun(const std::map<std::string, std::size_t, std::less<>>& avgSampleRuns);

    [[maybe_unused]] static void
    readInFilePathsFromDirectory(const std::string& inPath, std::vector<std::string>& codePaths);

    static CscMatrix toCsc(const gf2Mat& mat);
};
