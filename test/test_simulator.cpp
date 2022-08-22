//
// Created by luca on 09/08/22.
//
#include "Codes.hpp"
#include "Decoder.hpp"
#include "DecodingRunInformation.hpp"
#include "DecodingSimulator.hpp"
#include "UFDecoder.hpp"

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
    std::size_t nrRuns = 1, nrSamples=1;
    const std::string codePath = "./resources/codes/inCodes";
    auto        code = SteaneXCode();
    try {
        UFDecoder decoder;
        decoder.setCode(code);
        DecodingSimulator::simulateAverageRuntime(rawOut, testOut, physicalErrRate, nrRuns, codePath, nrSamples, DecoderType::UF_HEURISTIC);
    } catch (QeccException& e) {
        std::cerr << "Exception caught " << e.getMessage();
        EXPECT_TRUE(false);
    }
    EXPECT_TRUE(true);
}

TEST(DecodingSimulatorTest, TestPerformanceSim) {
    std::string rawOut = "./testRawFile", testOut = "./testStatFile";
    double      minErate = 0.01, maxErate = 0.03, stepSize = 0.01;
    std::size_t runsPerRate = 2;
    auto        code = SteaneXCode();
    try {
        DecodingSimulator::simulateWER(rawOut, testOut, minErate, maxErate,runsPerRate,code , stepSize, DecoderType::UF_DECODER);
    } catch (QeccException& e) {
        std::cerr << "Exception caught " << e.getMessage();
        EXPECT_TRUE(false);
    }
    EXPECT_TRUE(true);
}
