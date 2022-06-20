//
// Created by luca on 13/06/22.
//

#include "Codes.hpp"
#include "ImprovedUFD.hpp"

#include <gtest/gtest.h>

class ImprovedUFDtest: public testing::TestWithParam<std::vector<bool>> {
protected:
    void setUp() {
    }
};
INSTANTIATE_TEST_SUITE_P(CorrectableErrs, ImprovedUFDtest,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 1, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 1, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 1, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 0, 1, 0},
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 1}));

TEST_F(ImprovedUFDtest, SteaneCodeDecodingTestEstim) {
    SteaneXCode code{};
    ImprovedUFD decoder{code};
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = {1, 0, 0, 0, 0, 0, 0};

    auto syndr = code.getSyndrome(err);
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
    gf2Vec sol = {1, 0, 0, 0, 0, 0, 0};

    std::cout << "Estim: " << std::endl;
    Utils::printGF2vector(estim);
    std::cout << "EstimIdx: " << std::endl;
    Utils::printGF2vector(estim2);
    std::cout << "Sol: " << std::endl;
    Utils::printGF2vector(sol);
    EXPECT_TRUE(sol == estim);
    EXPECT_TRUE(sol == estim2);
}

TEST_P(ImprovedUFDtest, SteaneCodeDecodingTestFull) {
    SteaneXCode code{};
    ImprovedUFD decoder{code};
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = GetParam();

    auto syndr = code.getSyndrome(err);
    decoder.decode(syndr);
    auto              decodingResult = decoder.result;
    auto              estim          = decodingResult.estimNodeIdxVector;
    std::vector<bool> estimate(code.getN());
    for (auto e: estim) {
        estimate.at(e) = true;
    }
    std::vector<bool> residualErr(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = err[i] ^ estimate[i];
    }

    std::cout << "err: ";
    Utils::printGF2vector(err);
    std::cout << "syndr: ";
    Utils::printGF2vector(syndr);
    std::cout << "est: ";
    Utils::printGF2vector(estimate);
    std::cout << "resid: ";
    Utils::printGF2vector(residualErr);
    auto succ = code.isVectorStabilizer(residualErr);
    std::cout << "is in rowspace: " << succ << std::endl;
    EXPECT_TRUE(succ);
}