//
// Created by lucas on 13/06/22.
//

#include "Codes.hpp"
#include "OriginalUFD.hpp"

#include <gtest/gtest.h>

class OriginalUFDtest: public testing::TestWithParam<std::vector<bool>> {};
class UniquelyCorrectableErrTest_original: public OriginalUFDtest {};
class InCorrectableErrTest_original: public OriginalUFDtest {};
class UpToStabCorrectableErrTest_original: public OriginalUFDtest {};

INSTANTIATE_TEST_SUITE_P(CorrectableSingleBitErrs, UniquelyCorrectableErrTest_original,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 1, 0, 0, 0, 0}));

INSTANTIATE_TEST_SUITE_P(IncorrectableSingleBitErrs, InCorrectableErrTest_original,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 1, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 1, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 0, 1, 0}));

INSTANTIATE_TEST_SUITE_P(IncorrectableSingleBitErrs, UpToStabCorrectableErrTest_original,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 1},
                                 std::vector<bool>{1, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 1, 1, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 1}));
/**
 * Tests for unambigous syndromes, estimates must be computed exactly
 */
TEST_P(UniquelyCorrectableErrTest_original, SteaneCodeDecodingTestEstim) {
    auto code = new SteaneXCode();
    OriginalUFD decoder{code};
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = GetParam();

    auto syndr = code->getSyndrome(err);
    std::cout << "syndrome: ";
    Utils::printGF2vector(syndr);
    decoder.decode(syndr);
    auto   decodingResult = decoder.result;
    auto   estim          = decodingResult.estimBoolVector;
    auto   estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec estim2(err.size());
    std::cout << "estiIdxs: ";
    for (size_t i = 0; i < estimIdx.size(); i++) {
        estim2.at(estimIdx.at(i)) = true;
        std::cout << estimIdx.at(i);
    }
    std::cout << std::endl;
    gf2Vec sol = GetParam();

    std::cout << "Estim: " << std::endl;
    Utils::printGF2vector(estim);
    std::cout << "EstimIdx: " << std::endl;
    Utils::printGF2vector(estim2);
    std::cout << "Sol: " << std::endl;
    Utils::printGF2vector(sol);
    EXPECT_TRUE(sol == estim);
    EXPECT_TRUE(sol == estim2);
}

/**
 * Tests for ambigous errors that cannot be corrected
 */
TEST_P(InCorrectableErrTest_original, SteaneCodeDecodingTestEstim) {
    auto code = new SteaneXCode();
    OriginalUFD decoder{code};
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = GetParam();

    auto syndr = code->getSyndrome(err);
    decoder.decode(syndr);
    auto   decodingResult = decoder.result;
    auto   estim          = decodingResult.estimBoolVector;
    auto   estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec estim2(err.size());
    std::cout << "estiIdxs: ";
    for (size_t i = 0; i < estimIdx.size(); i++) {
        estim2.at(estimIdx.at(i)) = true;
        std::cout << estimIdx.at(i);
    }
    std::vector<bool> residualErr(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = err[i] ^ estim[i];
    }

    EXPECT_FALSE(Utils::isVectorInRowspace(code->Hz.pcm, residualErr));
}

/**
 * Tests for errors that are correctable up to stabilizer
 */
TEST_P(UpToStabCorrectableErrTest_original, SteaneCodeDecodingTest) {
    auto code = new SteaneXCode();
    OriginalUFD decoder{code};
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = GetParam();
    std::cout << "err :" << std::endl;
    Utils::printGF2vector(err);
    auto syndr = code->getSyndrome(err);
    decoder.decode(syndr);
    auto   decodingResult = decoder.result;
    auto   estim          = decodingResult.estimBoolVector;
    auto   estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec estim2(err.size());
    std::cout << "estiIdxs: ";
    for (size_t i = 0; i < estimIdx.size(); i++) {
        estim2.at(estimIdx.at(i)) = true;
        std::cout << estimIdx.at(i);
    }
    std::cout << std::endl;
    std::vector<bool> residualErr(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = err[i] ^ estim[i];
    }
    std::vector<bool> residualErr2(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr2.at(i) = err[i] ^ estim2[i];
    }

    EXPECT_TRUE(Utils::isVectorInRowspace(code->Hz.pcm, residualErr));
    EXPECT_TRUE(Utils::isVectorInRowspace(code->Hz.pcm, residualErr2));
}
