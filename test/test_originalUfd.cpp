// to keep 0/1 in boolean areas without clang-tidy warnings:
// NOLINTBEGIN(readability-implicit-bool-conversion,modernize-use-bool-literals)

#include "Codes.hpp"
#include "UFDecoder.hpp"
#include "Utils.hpp"

#include <cstddef>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

class OriginalUFDtest : public testing::TestWithParam<std::vector<bool>> {};
class UniquelyCorrectableErrTestOriginal : public OriginalUFDtest {};
class InCorrectableErrTestOriginal : public OriginalUFDtest {};
class UpToStabCorrectableErrTestOriginal : public OriginalUFDtest {};

INSTANTIATE_TEST_SUITE_P(CorrectableSingleBitErrsSteane, UniquelyCorrectableErrTestOriginal,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 1, 0, 0, 0, 0}));

INSTANTIATE_TEST_SUITE_P(UptoStabCorrSteane, InCorrectableErrTestOriginal,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 1, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 1, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 0, 1, 0}));

INSTANTIATE_TEST_SUITE_P(UptoStabCorrSteane, UpToStabCorrectableErrTestOriginal,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 1},
                                 std::vector<bool>{1, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 1, 1, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 1}));

/**
 * Tests for unambiguous syndromes, estimates must be computed exactly
 */
TEST_P(UniquelyCorrectableErrTestOriginal, SteaneCodeDecodingTestEstim) {
    auto      code = SteaneXCode();
    UFDecoder decoder;
    decoder.setCode(code);
    std::cout << "code:\n"
              << code << '\n';
    const std::vector<bool> err = GetParam();

    auto syndr = code.getXSyndrome(err);
    std::cout << "syndrome: ";
    Utils::printGF2vector(syndr);
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "\nestiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    const gf2Vec sol = GetParam();

    std::cout << "\nEstim:\n";
    Utils::printGF2vector(estim);
    std::cout << "\nEstimIdx: ";
    Utils::printGF2vector(estim2);
    std::cout << "\nSol: ";
    Utils::printGF2vector(sol);
    std::cout << '\n';
    EXPECT_TRUE(sol == estim);
    EXPECT_TRUE(sol == estim2);
}

/**
 * Tests for ambiguous errors that cannot be corrected
 */
TEST_P(InCorrectableErrTestOriginal, SteaneCodeDecodingTestEstim) {
    auto      code = SteaneXCode();
    UFDecoder decoder;
    decoder.setCode(code);
    std::cout << "code:\n"
              << code << '\n';
    std::vector<bool> err = GetParam();

    auto syndr = code.getXSyndrome(err);
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "\nestiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    std::vector<bool> residualErr(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = (err[i] != estim[i]);
    }

    const gf2Vec sol = GetParam();

    std::cout << "Estim: \n";
    Utils::printGF2vector(estim);
    std::cout << "EstimIdx: \n";
    Utils::printGF2vector(estim2);
    std::cout << "Sol: \n";
    Utils::printGF2vector(sol);
    EXPECT_FALSE(sol == estim);
    EXPECT_FALSE(sol == estim2);
}

/**
 * Tests for errors that are correctable up to stabilizer
 */
TEST_P(UpToStabCorrectableErrTestOriginal, SteaneCodeDecodingTest) {
    auto      code = SteaneCode();
    UFDecoder decoder;
    decoder.setCode(code);
    std::cout << "code:\n"
              << code << '\n';
    std::vector<bool> err = GetParam();
    std::cout << "err:\n";
    Utils::printGF2vector(err);
    auto syndr = code.getXSyndrome(err);
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "\nestiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    std::cout << '\n';
    std::vector<bool> residualErr(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = (err[i] != estim[i]);
    }
    std::vector<bool> residualErr2(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr2.at(i) = (err[i] != estim2[i]);
    }

    EXPECT_TRUE(Utils::isVectorInRowspace(*code.gethX()->pcm, residualErr));
    EXPECT_TRUE(Utils::isVectorInRowspace(*code.gethX()->pcm, residualErr2));
}
// NOLINTEND(readability-implicit-bool-conversion,modernize-use-bool-literals)
