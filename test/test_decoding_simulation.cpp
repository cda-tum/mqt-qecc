//
// Created by lucas on 21/04/2022.
//
#include "Codes.hpp"
#include "Decoder.hpp"
#include "DecodingRunInformation.hpp"
#include "DecodingSimulator.hpp"
#include "ImprovedUFD.hpp"

#include <bitset>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <locale>
#include <random>
using json = nlohmann::json;

class UnionFindSimulation: public testing::TestWithParam<std::string> {
protected:
    void setUp() {
    }
};

TEST(UnionFindSimulation, TestSimulator) {
    DecodingSimulator simulator;
    std::string       rawOut = "./testRawFile", testOut = "./testStatFile";
    double            minErate = 0.01, maxErate = 0.03, stepSize = 0.01;
    std::size_t       runsPerRate = 2, runsPerCode = 2;
    SteaneXCode       code;
    ImprovedUFD       decoder(code);
    simulator.simulateWER(rawOut, testOut, minErate, maxErate, stepSize, runsPerRate, decoder);
    std::vector<double> erRates = {minErate, maxErate};
    std::vector<Code>   codes   = {code};
    simulator.simulateRuntime(rawOut, testOut, erRates, runsPerCode, codes);
    EXPECT_TRUE(true);
}

/**
 * Simulate WER for growing number of physical err rate
 * Can also be used for threshold simulations
 */
TEST(UnionFindSimulation, EmpiricalEvaluationDecodingPerformance) {
    std::string outFilePath  = "./out";
    std::string dataFilePath = "./raw";
    std::cout << outFilePath;

    auto               t  = std::time(nullptr);
    auto               tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto          timestamp = oss.str();
    std::ofstream decodingResOutput(outFilePath + timestamp + ".json");
    std::ofstream rawDataOutput(dataFilePath + timestamp + ".json");
    // Basic Parameter setup
    const double                  normalizationConstant = 10'000.0; //
    double                        physicalErrRate       = 1.0 / normalizationConstant;
    double                        stepSize              = 1.0 / normalizationConstant;
    const double                  maxPhysicalErrRate    = 0.5;
    const std::size_t             nrOfRuns              = std::floor(maxPhysicalErrRate / physicalErrRate);
    std::size_t                   nrOfRunsPerRate       = 4096; // todo how deep to go?
    std::size_t                   nrOfFailedRuns        = 0U;
    double                        blockErrRate          = 0.0;
    double                        wordErrRate           = 0.0;
    std::size_t                   K                     = 0.0;
    std::map<std::string, double> wordErrRatePerPhysicalErrRate;
    decodingResOutput << "{ \"runs\" : [ ";
    for (std::size_t i = 0; i < nrOfRuns && physicalErrRate <= maxPhysicalErrRate; i++) {
        nrOfFailedRuns = 0U;
        blockErrRate   = 0.0;
        wordErrRate    = 0.0;
        decodingResOutput << R"({ "run": { "physicalErrRate":)" << physicalErrRate << ", \"data\": [ ";

        for (size_t j = 0; j < nrOfRunsPerRate; j++) {
            SteaneXCode code     = SteaneXCode();
            K = code.getK();
            ImprovedUFD decoder  = ImprovedUFD(code);
            auto        error    = Utils::sampleErrorIidPauliNoise(code.getN(), physicalErrRate);
            auto        syndrome = code.getSyndrome(error);
            decoder.decode(syndrome);
            auto              decodingResult = decoder.result;
            std::vector<bool> residualErr    = decodingResult.estimBoolVector;
            Utils::computeResidualErr(error, residualErr);
            auto                   success = code.isVectorStabilizer(residualErr);
            DecodingRunInformation stats;
            stats.result = decodingResult;

            if (success) {
                stats.status = SUCCESS;
                Utils::printGF2vector(residualErr);
            } else {
                stats.status = FAILURE;
                nrOfFailedRuns++;
            }
            decodingResOutput << decodingResult.to_json().dump(2U);
            if (j != nrOfRunsPerRate - 1) {
                decodingResOutput << ", ";
            }
        }
        //compute word error rate WER
        blockErrRate = (double)nrOfFailedRuns / (double)nrOfRunsPerRate;
        wordErrRate  = blockErrRate / (double)K;                                                            // rate of codewords re decoder does not give correct answer (fails or introduces logical operator)
        wordErrRatePerPhysicalErrRate.insert(std::make_pair(std::to_string(physicalErrRate), wordErrRate)); // to string for json parsing
        // only for json output
        if (i != nrOfRuns - 1) {
            decodingResOutput << "]}},";
        } else {
            decodingResOutput << "]}}";
        }
        physicalErrRate += stepSize;
    }
    decodingResOutput << "}";
    json dataj = wordErrRatePerPhysicalErrRate;
    rawDataOutput << dataj.dump(2U);
    decodingResOutput.close();
    rawDataOutput.close();
}

/**
 * Simulate average runtime for codes with growing nr of N for several physical err rates (err rates should only increase slope of curve)
 */
TEST(UnionFindSimulation, EmpiricalEvaluationDecoderRuntime) {
    std::string outFilePath  = "/home/luca/Documents/uf-simulations/runtime/out";
    std::string dataFilePath = "/home/luca/Documents/uf-simulations/runtime/data";
    std::cout << outFilePath;

    auto               t  = std::time(nullptr);
    auto               tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto          timestamp = oss.str();
    std::ofstream decodingResOutput(outFilePath + timestamp + ".json");
    std::ofstream rawDataOutput(dataFilePath + timestamp + ".json");
    // Basic Parameter setup
    const double                  physErrRates[]     = {0.01, 0.02, 0.03, 0.05, 0.1};
    std::size_t                   avgDecodingTimeAcc = 0U;
    const std::size_t             nrOfTrials         = 1'0;
    double                        avgDecTime         = 0.0;
    std::map<std::string, double> avgDecodingTimePerSize;
    std::vector<Code>             codes;
    std::string                   inPath = "/home/luca/Documents/codeRepos/qunionfind/examples/test/";
    for (const auto& file: std::filesystem::directory_iterator(inPath)) {
        codes.emplace_back(Code(ParityCheckMatrix(Utils::importGf2MatrixFromFile(file.path()))));
    }

    for (auto physErrRate: physErrRates) {
        for (const auto& code: codes) {
            avgDecodingTimeAcc = 0U;
            std::cout << "Code: " << Utils::getStringFrom(code.Hz.pcm) << std::endl;
            for (std::size_t i = 0; i < nrOfTrials; i++) {
                auto        c        = Code(code.Hz); // construct new for each trial
                ImprovedUFD decoder  = ImprovedUFD(c);
                auto        error    = Utils::sampleErrorIidPauliNoise(c.getN(), physErrRate);
                auto        syndrome = c.getSyndrome(error);

                std::cout << "error: " << Utils::getStringFrom(error) << ", rate: " << physErrRate << std::endl;
                std::cout << "starting decoding" << std::endl;
                decoder.decode(syndrome);
                std::cout << "decoding done " << std::endl;

                auto decodingResult = decoder.result;

                std::cout << "runtime: " << decodingResult.decodingTime << std::endl
                          << std::endl;

                avgDecodingTimeAcc = avgDecodingTimeAcc + decodingResult.decodingTime;
                decodingResOutput << decodingResult.to_json();
            }
            avgDecTime = (double)avgDecodingTimeAcc / (double)nrOfTrials;
            avgDecodingTimePerSize.insert(std::make_pair<>(std::to_string(code.getN()), avgDecTime));
        }
    }
    json dataj = avgDecodingTimePerSize;
    rawDataOutput << dataj.dump(2U);
}
