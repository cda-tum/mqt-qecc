//
// Created by luca on 13/06/22.
//

#include "Codes.hpp"
#include "ImprovedUFD.hpp"

#include <gtest/gtest.h>

class ImprovedUFDtest: public testing::TestWithParam<std::string> {
protected:
    void setUp() {
    }
};

TEST(ImprovedUFDtest, SteaneCodeDecodingTest) {
    SteaneXCode code{};
    ImprovedUFD decoder{code};
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = {0, 0, 0, 0, 1, 0, 0};

    auto syndr = code.getSyndrome(err);
    decoder.decode(syndr);
    auto decodingResult = decoder.result;
    auto estim          = decodingResult.estimNodeIdxVector;

    EXPECT_TRUE(!estim.empty());
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