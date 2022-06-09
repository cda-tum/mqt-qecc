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

std::vector<bool> dummySampler(const std::size_t n) {
    std::vector<bool>               result(n);
    std::random_device              rd;
    std::mt19937                    gen(rd());
    std::uniform_int_distribution<> distr(0, n);

    result.at(0) = true;

    return result;
}

TEST(UnionFindSimulation, SteaneCodeDecodingTest) {
    SteaneXCode code{};
    ImprovedUFD decoder{code};
    std::cout << "code: " << std::endl
              << code << std::endl;
    auto err = dummySampler(code.getN());
    //std::cout << "error: " << err << std::endl;
    auto syndr = code.getSyndrome(err);
    //std::cout << "syndr: " << syndr << std::endl;
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

    //std::cout << "estim: " << estimate << std::endl;
    //std::cout << "r = e'+e " << residualErr << std::endl;
    auto succ = code.checkStabilizer(residualErr);
    EXPECT_TRUE(succ);

    if (succ) {
        std::cout << "Decoding successful, found estimate up to stabilizer: " << std::endl;
        //std::cout << estimate << std::endl;
        std::cout << "Elapsed time: " << decodingResult.decodingTime << "ms" << std::endl;
    } else {
        decodingResult.status = FAILURE;
        std::cout << "Decoding not successful, introduced logical opertor" << std::endl;
    }
}

std::vector<bool> sampleErrorIidPauliNoise(const std::size_t n, const double physicalErrRate) {
    std::random_device rd;
    std::mt19937       gen(rd());
    std::vector<bool>  result;

    // Setup the weights, iid noise for each bit
    std::discrete_distribution<> d({1 - physicalErrRate, physicalErrRate});
    for (std::size_t i = 0; i < n; i++) {
        result.emplace_back(d(gen));
    }
    return result;
}
/**
 *
 * @param error bool vector representing error
 * @param residual estimate vector that contains residual error at end of function
 */
void computeResidualErr(const std::vector<bool>& error, std::vector<bool>& residual) {
    for (std::size_t j = 0; j < residual.size(); j++) {
        residual.at(j) = residual.at(j) ^ error.at(j);
    }
}


// main simulation for empirical evaluation study
TEST(UnionFindSimulation, EmpiricalEvaluation) {
    std::string        outFilePath  = "/home/luca/Documents/uf-simulations/testrun/out";
    std::string        dataFilePath = "/home/luca/Documents/uf-simulations/testrun/data";
    std::cout << outFilePath;

    auto               t            = std::time(nullptr);
    auto               tm           = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto          timestamp = oss.str();
    std::ofstream decodingResOutput(outFilePath + timestamp + ".json");
    std::ofstream rawDataOutput(dataFilePath + timestamp + ".json");
    const double                  normalizationConstant = 10'000.0; //
    double                        physicalErrRate       = 1.0 / normalizationConstant;
    double                        stepSize              = 10.0 / normalizationConstant;
    const double                  maxPhysicalErrRate    = 0.5;
    const size_t                  nrOfRuns              = std::floor(maxPhysicalErrRate/physicalErrRate);
    std::size_t                   nrOfRunsPerRate       = 32;
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
            auto       error           = sampleErrorIidPauliNoise(code.getN(), physicalErrRate);
            auto syndrome = code.getSyndrome(error);
            decoder.decode(syndrome);
            auto              decodingResult = decoder.result;
            std::vector<bool> residualErr    = decodingResult.estimBoolVector;
            computeResidualErr(error, residualErr);
            auto success = code.checkStabilizer(residualErr);

            std::cout << "error: " ;
            Utils::printGF2vector(error);
            if (success) {
                decodingResult.status = SUCCESS;
                std::cout << "Decoding successful: " << std::endl;
                Utils::printGF2vector(residualErr);
                std::cout << "Elapsed time: " << decodingResult.decodingTime << "ms" << std::endl;
            } else {
                decodingResult.status = FAILURE;
                nrOfFailedRuns++;
                std::cout << "Decoding failure" << std::endl;
            }
            decodingResOutput << decodingResult.to_json().dump(2U);
            if (j != nrOfRunsPerRate - 1) {
                decodingResOutput << ", ";
            }
        }
        //compute word error rate WER
        blockErrRate = (double)nrOfFailedRuns / (double)nrOfRunsPerRate;
        wordErrRate  = blockErrRate / (double)K; // rate of codewords re decoder does not give correct answer (fails or introduces logical operator)
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
