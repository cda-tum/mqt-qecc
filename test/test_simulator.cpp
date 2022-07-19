//
// Created by lucas on 21/04/2022.
//
#include "Codes.hpp"
#include "Decoder.hpp"
#include "DecodingRunInformation.hpp"
#include "DecodingSimulator.hpp"
#include "ImprovedUFD.hpp"
#include "OriginalUFD.hpp"

#include <bitset>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <locale>
#include <random>
using json = nlohmann::json;
class DecodingSimulatorTest: public testing::TestWithParam<std::string> {
};

TEST(DecodingSimulatorTest, TestRuntimeSim) {
    std::string rawOut = "./testRawFile", testOut = "./testStatFile";
    const double physicalErrRate = 0.01;
    std::size_t nrRuns = 2, nrSamples=2;
    const std::string codePath = "./resources/codes/inCodes";
    auto        code = SteaneXCode();
    try {
        ImprovedUFD decoder;
        decoder.setCode(code);
        DecodingSimulator::simulateAverageRuntime(rawOut, testOut, physicalErrRate, nrRuns, codePath, nrSamples);
    } catch (QeccException& e) {
        std::cerr << "Exception caught " << e.getMessage();
        EXPECT_TRUE(false);
    }
    EXPECT_TRUE(true);
}

TEST(DecodingSimulatorTest, TestPerformanceSim) {
    std::string rawOut = "./testRawFile", testOut = "./testStatFile";
    double      minErate = 0.01, maxErate = 0.03, stepSize = 0.01;
    std::size_t runsPerRate = 2, runsPerCode = 2;
    auto        code = SteaneXCode();
    try {
        ImprovedUFD decoder;
        decoder.setCode(code);
        DecodingSimulator::simulateWER(rawOut, testOut, minErate, maxErate, stepSize, runsPerRate, decoder);
    } catch (QeccException& e) {
        std::cerr << "Exception caught " << e.getMessage();
        EXPECT_TRUE(false);
    }
    EXPECT_TRUE(true);
}
