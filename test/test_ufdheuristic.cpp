//
// Created by lucas on 13/06/22.
//
// to keep 0/1 in boolean areas without clang-tidy warnings:
// NOLINTBEGIN(readability-implicit-bool-conversion,modernize-use-bool-literals)

#include "Codes.hpp"
#include "UFHeuristic.hpp"

#include <gtest/gtest.h>
class ImprovedUFDtestBase : public testing::TestWithParam<std::vector<bool>> {};
class UniquelyCorrectableErrTest : public ImprovedUFDtestBase {};
class IncorrectableErrTest : public ImprovedUFDtestBase {};
class UpToStabCorrectableErrTest : public ImprovedUFDtestBase {};
class CorrectableLargeToric : public ImprovedUFDtestBase {};

INSTANTIATE_TEST_SUITE_P(CorrectableSingleBitErrs, UniquelyCorrectableErrTest,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 1, 0, 0, 0, 0}));

INSTANTIATE_TEST_SUITE_P(IncorrectableSingleBitErrs, IncorrectableErrTest,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 1, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 1, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 0, 1, 0}));

INSTANTIATE_TEST_SUITE_P(UptoStabCorrectable, UpToStabCorrectableErrTest,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 1},
                                 std::vector<bool>{1, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 1, 1, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 1}));

INSTANTIATE_TEST_SUITE_P(CorrectableLargeToricTests, CorrectableLargeToric,
                         testing::Values(
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}));
/**
 * Tests for unambiguous syndromes, estimates must be computed exactly
 */
TEST_P(UniquelyCorrectableErrTest, SteaneCodeDecodingTestEstim) {
    auto        code = SteaneXCode();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "code: " << std::endl
              << code << std::endl;
    const std::vector<bool> err = GetParam();

    auto syndr = code.getXSyndrome(err);
    std::cout << "syndrome: " << Utils::getStringFrom(syndr) << std::endl;
    decoder.decode(syndr);
    auto const& decodingResult = decoder.result;
    auto const& estim          = decodingResult.estimBoolVector;
    auto const& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    std::cout << std::endl;
    const gf2Vec sol = GetParam();

    std::cout << "Estim: " << Utils::getStringFrom(estim) << std::endl;
    std::cout << "Sol: " << std::endl;
    Utils::printGF2vector(sol);
    EXPECT_TRUE(sol == estim);
    EXPECT_TRUE(sol == estim2);
}

/**
 * Tests for ambiguous errors that cannot be corrected
 */
TEST_P(IncorrectableErrTest, SteaneCodeDecodingTestEstim2) {
    auto        code = SteaneXCode();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err   = GetParam();
    auto              syndr = code.getXSyndrome(err);
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    std::cout << "estim: " << Utils::getStringFrom(estim) << std::endl;
    const auto& estimIdx = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    std::vector<bool> residualErr(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = (err[i] != estim[i]);
    }

    EXPECT_FALSE(Utils::isVectorInRowspace(*(code.gethZ()->pcm), residualErr));
}

/**
 * Tests for errors that are correctable up to stabilizer
 */
TEST_P(UpToStabCorrectableErrTest, SteaneCodeDecodingTest) {
    auto        code = SteaneXCode();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = GetParam();
    std::cout << "err :" << std::endl;
    Utils::printGF2vector(err);
    std::cout << std::endl;
    auto syndr = code.getXSyndrome(err);
    Utils::printGF2vector(syndr);
    std::cout << std::endl;
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
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
    std::cout << "estim: " << Utils::getStringFrom(estim) << std::endl;
    std::cout << "resid: " << Utils::getStringFrom(residualErr) << std::endl;
    EXPECT_TRUE(Utils::isVectorInRowspace(*code.gethZ()->pcm, residualErr));
    EXPECT_TRUE(Utils::isVectorInRowspace(*code.gethZ()->pcm, residualErr2));
}

/**
 * Tests for toric code one bit correctable errs
 */
TEST_F(ImprovedUFDtestBase, UniquelyCorrectableErrLargeToricCodeTest) {
    auto        code = ToricCode32();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "Adj lists code: " << std::endl
              << Utils::getStringFrom(*code.gethZ()->pcm) << std::endl;
    const std::vector<bool> err = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
    std::cout << "error: ";
    Utils::printGF2vector(err);
    std::cout << std::endl;
    auto syndr = code.getXSyndrome(err);
    std::cout << "syndrome: ";
    Utils::printGF2vector(syndr);
    std::cout << std::endl;
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx << "; ";
    }
    EXPECT_TRUE(estim == estim2);

    std::cout << std::endl;
    const gf2Vec& sol = err;
    std::cout << "Estim: " << Utils::getStringFrom(estim) << std::endl;
    std::cout << "Estim from Idx: " << Utils::getStringFrom(estim2) << std::endl;
    std::cout << "Sol: " << Utils::getStringFrom(sol) << std::endl;
    EXPECT_TRUE(sol == estim2);
}

/**
 * Tests for toric code one bit correctable errs
 */
TEST_P(CorrectableLargeToric, UniquelyCorrectableErrLargeToricCodeTest2) {
    auto        code = ToricCode32();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "Adj lists code: " << std::endl
              << Utils::getStringFrom(*code.gethZ()->pcm) << std::endl;
    const std::vector<bool> err = GetParam();
    std::cout << "error: ";
    Utils::printGF2vector(err);
    std::cout << std::endl;
    auto syndr = code.getXSyndrome(err);
    std::cout << "syndrome: ";
    Utils::printGF2vector(syndr);
    std::cout << std::endl;
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx << "; ";
    }
    EXPECT_TRUE(estim == estim2);

    std::cout << std::endl;
    const gf2Vec& sol = err;
    std::cout << "Estim: " << Utils::getStringFrom(estim) << std::endl;
    std::cout << "Estim from Idx: " << Utils::getStringFrom(estim2) << std::endl;
    std::cout << "Sol: " << Utils::getStringFrom(sol) << std::endl;
    EXPECT_TRUE(sol == estim2);
}
/**
 * Tests for errors that are correctable up to stabilizer
 */
TEST_F(ImprovedUFDtestBase, LargeCodeTest) {
    auto        code = HGPcode();
    UFHeuristic decoder;
    decoder.setCode(code);
    auto err  = gf2Vec(code.n);
    err.at(0) = 1;

    std::cout << "err :" << std::endl;
    Utils::printGF2vector(err);
    auto syndr = code.getXSyndrome(err);
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    std::cout << std::endl;
    std::vector<bool> residualErr(err.size());
    for (std::size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = (err.at(i) != estim.at(i));
    }
    std::vector<bool> residualErr2(err.size());
    for (std::size_t i = 0; i < err.size(); i++) {
        residualErr2.at(i) = (err.at(i) != estim2.at(i));
    }

    EXPECT_TRUE(Utils::isVectorInRowspace(*code.gethZ()->pcm, residualErr));
    EXPECT_TRUE(Utils::isVectorInRowspace(*code.gethZ()->pcm, residualErr2));
}
TEST_F(ImprovedUFDtestBase, BothErrsTest) {
    try {
        auto        code = SteaneCode();
        UFHeuristic decoder;
        decoder.setCode(code);
        std::vector<bool> err(code.n * 2);
        err.at(0)           = 1;
        err.at(code.getN()) = 1;

        std::cout << "err :" << std::endl;
        Utils::printGF2vector(err);
        auto syndr = code.getXSyndrome(err);
        decoder.decode(syndr);
        const auto& decodingResult = decoder.result;
        const auto& estim          = decodingResult.estimBoolVector;
        const auto& estimIdx       = decodingResult.estimNodeIdxVector;
        gf2Vec      estim2(err.size());
        std::cout << "estiIdxs: ";
        for (auto idx : estimIdx) {
            estim2.at(idx) = true;
            std::cout << idx;
        }
        std::cout << std::endl;
        std::vector<bool> residualErr(err.size());
        for (std::size_t i = 0; i < err.size(); i++) {
            residualErr.at(i) = (err.at(i) != estim.at(i));
        }
        std::vector<bool> residualErr2(err.size());
        for (std::size_t i = 0; i < err.size(); i++) {
            residualErr2.at(i) = (err.at(i) != estim2.at(i));
        }
        EXPECT_TRUE(code.isStabilizer(residualErr));
    } catch (QeccException& e) {
        std::cout << e.getMessage() << std::endl;
        EXPECT_TRUE(false);
    }
}
// NOLINTEND(readability-implicit-bool-conversion,modernize-use-bool-literals)
