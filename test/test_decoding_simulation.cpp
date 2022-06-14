//
// Created by lucas on 21/04/2022.
//
#include "Codes.hpp"
#include "Decoder.hpp"
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

/**
 * Simulate WER for growing number of physical err rate
 * Can also be used for threshold simulations
 */
TEST(UnionFindSimulation, EmpiricalEvaluationDecodingPerformance) {
    std::string outFilePath  = "/home/luca/Documents/uf-simulations/decoding/out";
    std::string dataFilePath = "/home/luca/Documents/uf-simulations/decoding/data";
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
    const size_t                  nrOfRuns              = std::floor(maxPhysicalErrRate / physicalErrRate);
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
            SteaneXCode code{};
            K = code.getK();
            ImprovedUFD decoder{code};
            auto        error    = Utils::sampleErrorIidPauliNoise(code.getN(), physicalErrRate);
            auto        syndrome = code.getSyndrome(error);
            decoder.decode(syndrome);
            auto              decodingResult = decoder.result;
            std::vector<bool> residualErr    = decodingResult.estimBoolVector;
            Utils::computeResidualErr(error, residualErr);
            auto success = code.isVectorStabilizer(residualErr);

            if (success) {
                decodingResult.status = SUCCESS;
                Utils::printGF2vector(residualErr);
            } else {
                decodingResult.status = FAILURE;
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

    const double                  physErrRates[]     = {0.0001, 0.0002, 0.0005, 0.001, 0.05};
    std::size_t                   avgDecodingTimeAcc = 0U;
    const std::size_t             nrOfTrials         = 1'000'000;
    double                        avgDecTime         = 0.0;
    std::map<std::string, double> avgDecodingTimePerSize;

    for (auto physErrRate: physErrRates) {
        Code codes[] = {};
        for (const auto& code: codes) {
            avgDecodingTimeAcc = 0U;
            for (size_t i = 0; i < nrOfTrials; i++) {
                auto        c = Code(code.Hz); // construct new for each trial
                ImprovedUFD decoder{c};
                auto        error    = Utils::sampleErrorIidPauliNoise(c.getN(), physErrRate);
                auto        syndrome = c.getSyndrome(error);
                decoder.decode(syndrome);
                auto decodingResult = decoder.result;
                avgDecodingTimeAcc  = avgDecodingTimeAcc + decodingResult.decodingTime;
            }
            avgDecTime = (double)avgDecodingTimeAcc / (double)nrOfTrials;
            avgDecodingTimePerSize.insert(std::make_pair<>(std::to_string(code.getN()), avgDecTime));
        }
    }
    json dataj = avgDecodingTimePerSize;
    rawDataOutput << dataj.dump(2U);
}
