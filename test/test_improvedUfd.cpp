//
// Created by lucas on 13/06/22.
//

#include "Codes.hpp"
#include "UFHeuristic.hpp"

#include <gtest/gtest.h>
class ImprovedUFDtestBase: public testing::TestWithParam<std::vector<bool>> {};
class UniquelyCorrectableErrTest: public ImprovedUFDtestBase {};
class UniquelyCorrectableErrToricCodeTest: public ImprovedUFDtestBase {}; // test cases wrong not uniquely determined
class IncorrectableErrToricCodeTest: public ImprovedUFDtestBase {};       // test cases not uniquely fixed
class IncorrectableErrTest: public ImprovedUFDtestBase {};
class UpToStabCorrectableErrTest: public ImprovedUFDtestBase {}; // first case not dec
class CorrectableLargeToric: public ImprovedUFDtestBase {};

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
                                 std::vector<bool>{0, 0, 0, 0, 0, 1, 0},
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 1}));

INSTANTIATE_TEST_SUITE_P(UptoStabCorrectable, UpToStabCorrectableErrTest,
                         testing::Values(
                                 std::vector<bool>{1, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 1, 1, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 1}));

INSTANTIATE_TEST_SUITE_P(CorrectableSingleBitErrsToric, UniquelyCorrectableErrToricCodeTest,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 0, 1},
                                 std::vector<bool>{0, 0, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 1, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 0, 1, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 0, 1}));

INSTANTIATE_TEST_SUITE_P(NotorrectableSingleBitErrsToric, IncorrectableErrToricCodeTest,
                         testing::Values(
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 0, 0}));

INSTANTIATE_TEST_SUITE_P(CorrectableLargeToricTests, CorrectableLargeToric,
                         testing::Values(
                                 // std::vector<bool>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}
                                 // std::vector<bool>{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 //  std::vector<bool>{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 //  std::vector<bool>{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 //  std::vector<bool>{1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                 //  std::vector<bool>{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                                 //  std::vector<bool>{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}
                                 ));
/**
 * Tests for unambigous syndromes, estimates must be computed exactly
 */
TEST_P(UniquelyCorrectableErrTest, SteaneCodeDecodingTestEstim) {
    auto        code = SteaneXCode();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = GetParam();

    auto syndr = code.getSyndrome(err);
    std::cout << "syndrome: " << Utils::getStringFrom(syndr) << std::endl;
    decoder.decode(syndr);
    auto const& decodingResult = decoder.result;
    auto const& estim          = decodingResult.estimBoolVector;
    auto const& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (unsigned long idx: estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    std::cout << std::endl;
    gf2Vec sol = GetParam();

    std::cout << "Estim: " << Utils::getStringFrom(estim) << std::endl;
    std::cout << "Sol: " << std::endl;
    Utils::printGF2vector(sol);
    EXPECT_TRUE(sol == estim);
    EXPECT_TRUE(sol == estim2);
}

/**
 * Tests for ambigous errors that cannot be corrected
 */
TEST_P(IncorrectableErrTest, SteaneCodeDecodingTestEstim2) {
    auto        code = SteaneXCode();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err   = GetParam();
    auto              syndr = code.getSyndrome(err);
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    std::cout << "estim: " << Utils::getStringFrom(estim) << std::endl;
    const auto& estimIdx = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (unsigned long idx: estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    std::vector<bool> residualErr(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = (err[i] != estim[i]);
    }

    EXPECT_FALSE(Utils::isVectorInRowspace(*(code.getHz()->pcm), residualErr));
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
    auto syndr = code.getSyndrome(err);
    Utils::printGF2vector(syndr);
    std::cout << std::endl;
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (unsigned long idx: estimIdx) {
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
    EXPECT_TRUE(Utils::isVectorInRowspace(*code.getHz()->pcm, residualErr));
    EXPECT_TRUE(Utils::isVectorInRowspace(*code.getHz()->pcm, residualErr2));
}

/**
 * Tests for toric code one bit correctable errs
 */
TEST_P(UniquelyCorrectableErrToricCodeTest, ToricCodeTest) {
    auto        code = ToricCode_8();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "Adj lists code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = GetParam();
    std::cout << "error: ";
    Utils::printGF2vector(err);
    std::cout << std::endl;
    auto syndr = code.getSyndrome(err);
    std::cout << "syndrome: ";
    Utils::printGF2vector(syndr);
    std::cout << std::endl;
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (unsigned long idx: estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx << "; ";
    }
    std::cout << std::endl;
    const gf2Vec& sol = err;

    std::cout << "Estim: " << Utils::getStringFrom(estim) << std::endl;
    std::cout << "Estim from Idx: " << Utils::getStringFrom(estim2) << std::endl;
    std::cout << "Sol: " << Utils::getStringFrom(sol) << std::endl;
    EXPECT_TRUE(sol == estim);
    EXPECT_TRUE(sol == estim2);
}
/**
 * Tests for toric code one bit not uniquely corr errs
 */
TEST_P(IncorrectableErrToricCodeTest, ToricCodeTest2) {
    auto        code = ToricCode_8();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "Adj lists code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = GetParam();
    std::cout << "error: ";
    Utils::printGF2vector(err);
    std::cout << std::endl;
    auto syndr = code.getSyndrome(err);
    std::cout << "syndrome: ";
    Utils::printGF2vector(syndr);
    std::cout << std::endl;
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (unsigned long idx: estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx << "; ";
    }
    std::cout << std::endl;
    gf2Vec            sol = err;
    std::vector<bool> residualErr(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = err[i] ^ estim[i];
    }

    std::cout << "Estim: " << Utils::getStringFrom(estim) << std::endl;
    std::cout << "Estim from Idx: " << Utils::getStringFrom(estim2) << std::endl;
    std::cout << "Sol: " << Utils::getStringFrom(sol) << std::endl;
    EXPECT_FALSE(Utils::isVectorInRowspace(*code.getHz()->pcm, residualErr));
}

/**
 * Tests for toric code one bit correctable errs
 */
TEST_F(ImprovedUFDtestBase, UniquelyCorrectableErrLargeToricCodeTest) {
    auto        code = ToricCode_32();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "Adj lists code: " << std::endl
              << Utils::getStringFrom(*code.getHz()->pcm) << std::endl;
    const std::vector<bool> err = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
    std::cout << "error: ";
    Utils::printGF2vector(err);
    std::cout << std::endl;
    auto syndr = code.getSyndrome(err);
    std::cout << "syndrome: ";
    Utils::printGF2vector(syndr);
    std::cout << std::endl;
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (unsigned long idx: estimIdx) {
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
    auto        code = ToricCode_32();
    UFHeuristic decoder;
    decoder.setCode(code);
    std::cout << "Adj lists code: " << std::endl
              << Utils::getStringFrom(*code.getHz()->pcm) << std::endl;
    std::vector<bool> err = GetParam();
    std::cout << "error: ";
    Utils::printGF2vector(err);
    std::cout << std::endl;
    auto syndr = code.getSyndrome(err);
    std::cout << "syndrome: ";
    Utils::printGF2vector(syndr);
    std::cout << std::endl;
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (unsigned long idx: estimIdx) {
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
    auto        err = gf2Vec(code.N);
    err.at(0)       = 1;

    std::cout << "err :" << std::endl;
    Utils::printGF2vector(err);
    auto syndr = code.getSyndrome(err);
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (unsigned long idx: estimIdx) {
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

    EXPECT_TRUE(Utils::isVectorInRowspace(*code.getHz()->pcm, residualErr));
    EXPECT_TRUE(Utils::isVectorInRowspace(*code.getHz()->pcm, residualErr2));
}
