#include "Codes.hpp"
#include "DecodingSimulator.hpp"
#include "QeccException.hpp"
#include "UFDecoder.hpp"

#include <cstddef>
#include <gtest/gtest.h>
#include <iostream>
#include <string>

class DecodingSimulatorTest : public testing::TestWithParam<std::string> {
};

TEST(DecodingSimulatorTest, TestRuntimeSim) {
    const std::string rawOut          = "./testRawFile";
    const std::string testOut         = "./testStatFile";
    const double      physicalErrRate = 0.01;
    const std::size_t nrRuns          = 1;
    const std::size_t nrSamples       = 1;
    const std::string codePath        = "./resources/codes/inCodes";
    auto              code            = SteaneXCode();
    try {
        UFDecoder decoder;
        decoder.setCode(code);
        DecodingSimulator::simulateAverageRuntime(rawOut, testOut, physicalErrRate, nrRuns, codePath, nrSamples, DecoderType::UfHeuristic);
    } catch (QeccException& e) {
        std::cerr << "Exception caught " << e.getMessage();
        EXPECT_TRUE(false);
    }
    EXPECT_TRUE(true);
}

TEST(DecodingSimulatorTest, TestPerformanceSim) {
    const std::string rawOut      = "./testRawFile";
    const std::string testOut     = "./testStatFile";
    const double      minErate    = 0.01;
    const double      maxErate    = 0.03;
    const double      stepSize    = 0.01;
    const std::size_t runsPerRate = 2;
    auto              code        = SteaneCode();
    try {
        DecodingSimulator::simulateWER(rawOut, testOut, minErate, maxErate, runsPerRate, code, stepSize, DecoderType::UfDecoder);
    } catch (QeccException& e) {
        std::cerr << "Exception caught " << e.getMessage();
        EXPECT_TRUE(false);
    }
    EXPECT_TRUE(true);
}
