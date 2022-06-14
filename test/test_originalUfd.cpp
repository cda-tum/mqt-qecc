//
// Created by luca on 13/06/22.
//

#include "Codes.hpp"
#include "OriginalUFD.hpp"

#include <gtest/gtest.h>

class OriginalUFDtest: public testing::TestWithParam<std::string> {
protected:
    void setUp() {
    }
};

TEST(OriginalUFDtest, SteaneCodeDecodingTest) {
    SteaneXCode code{};
    OriginalUFD decoder{code};
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err   = {1, 0, 0, 0, 0, 0, 0};
    auto              syndr = code.getSyndrome(err);
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
    std::cout << "estim: ";
    Utils::printGF2vector(estimate);
    auto succ = code.isVectorStabilizer(residualErr);
    EXPECT_TRUE(succ);
}
