//
// Created by lucas on 21/04/2022.
//
#include "Codes.hpp"
#include "Decoder.hpp"

#include <bitset>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <locale>
#include <random>

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

    //result.at(distr(gen)) = true;
    result.at(0) = true;

    return result;
}

TEST(SteaneCodeTest, SteaneCodeDecoding) {
    SteaneXCode code{};
    Decoder     decoder{code};
    std::cout << "code: " << std::endl
              << code << std::endl;
    auto err = dummySampler(code.getN());
    std::cout << "error: " << err << std::endl;
    auto syndr = code.getSyndrome(err);
    std::cout << "syndr: " << syndr << std::endl;
    std::set<std::shared_ptr<TreeNode>> syndrComponents;
    for (size_t i = 0; i < syndr.size(); i++) {
        if (syndr.at(i)) {
            auto snode     = code.tannerGraph.adjListNodes.at(i + code.getN()).at(0);
            snode->isCheck = true;
            snode->checkVertices.insert(snode->vertexIdx);
            syndrComponents.insert(snode);
        }
    }
    std::cout << "s component " << syndrComponents << std::endl;

    decoder.decode(syndrComponents);
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

    std::cout << "estim: " << estimate << std::endl;
    std::cout << "r = e'+e " << residualErr << std::endl;
    auto succ = code.checkStabilizer(residualErr);
    EXPECT_TRUE(succ);

    if (succ) {
        std::cout << "Decoding successful, found estimate up to stabilizer: " << std::endl;
        std::cout << estimate << std::endl;
        std::cout << "Elapsed time: " << decodingResult.decodingTime << "ms" << std::endl;
    } else {
        decodingResult.status = FAILURE;
        std::cout << "Decoding not successful, introduced logical opertor" << std::endl;
    }
}
std::vector<bool> sampleError(const std::size_t n, const double physicalErrRate) {
    std::random_device rd;
    std::mt19937       gen(rd());
    std::vector<bool>  result;

    // Setup the weights, iid for each qubit
    std::discrete_distribution<> d({1 - physicalErrRate, physicalErrRate});
    for (std::size_t i = 0; i < n; i++) {
        result.emplace_back(d(gen));
    }
    return result;
}
/**
 *
 * @param error bool vector representing error
 * @param residual estimate vector that contains residual error at end
 */
void computeResidualErr(const std::vector<bool>& error, std::vector<bool>& residual) {
    for (std::size_t j = 0; j < residual.size(); j++) {
        residual.at(j) = residual.at(j) ^ error.at(j);
    }
}
/**
 * returns list of tree node representations for syndrome
 * @param code
 * @param syndrome
 * @return
 */
std::set<std::shared_ptr<TreeNode>> computeInitTreeComponents(Code& code, const std::vector<bool>& syndrome) {
    std::set<std::shared_ptr<TreeNode>> result;
    for (size_t j = 0; j < syndrome.size(); j++) {
        if (syndrome.at(j)) {
            auto syndrNode     = code.tannerGraph.adjListNodes.at(j + code.getN()).at(0);
            syndrNode->isCheck = true;
            syndrNode->checkVertices.insert(syndrNode->vertexIdx);
            result.insert(syndrNode);
        }
    }
}

TEST(LogicalErrorRateSimulation, SteaneCodeDecoding) {
    std::string        outFilePath = "path";
    std::string        dataFilePath = "path";
    auto               t           = std::time(nullptr);
    auto               tm          = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto          timestamp = oss.str();
    std::ofstream            decodingResOutput(outFilePath + timestamp + ".json");
    std::ofstream            rawDataOutput(dataFilePath + timestamp + ".json");
    double        physicalErrRate = 0.0001;
    double maxPhysicalErrRate = 0.1;
    double stepSize = 0.0001;
    std::size_t   nrOfRunsPerRate = 10;
    int nrOfFailedRuns = 0;
    double blockErrRate = 0.0;
    double wordErrRate = 0.0;
    double K = 0.0;
    std::map<double, double> wordErrRateData;
    decodingResOutput << "{ \"runs\" : [ ";

    for (; physicalErrRate < maxPhysicalErrRate; physicalErrRate += stepSize) {
        nrOfFailedRuns = 0;
        blockErrRate = wordErrRate = 0.0;
        decodingResOutput << "{ \"run\": { \"physicalErrRate\":" << physicalErrRate << ", \"data\": [ ";
        for (size_t i = 0; i < nrOfRunsPerRate; i++) {
            SteaneXCode code{};
            K = code.getK();
            Decoder     decoder{code};
            auto        error           = sampleError(code.getN(), physicalErrRate);
            auto        syndrComponents = computeInitTreeComponents(code, code.getSyndrome(error));
            decoder.decode(syndrComponents);
            auto              decodingResult = decoder.result;
            std::vector<bool> residualErr    = decodingResult.estimBoolVector;
            computeResidualErr(error, residualErr);
            auto success = code.checkStabilizer(residualErr);

            if (success) {
                decodingResult.status = SUCCESS;
                std::cout << "Decoding successful, found residualErr up to stabilizer: " << std::endl;
                std::cout << residualErr << std::endl;
                std::cout << "Elapsed time: " << decodingResult.decodingTime << "ms" << std::endl;
            } else {
                decodingResult.status = FAILURE;
                nrOfFailedRuns++;
                std::cout << "Decoding not successful, introduced logical opertor" << std::endl;
            }
            decodingResOutput << decodingResult.to_json().dump(2U);
            if(i != nrOfRunsPerRate-1) {
                decodingResOutput << ", ";
            }
        }
        blockErrRate = nrOfFailedRuns/nrOfRunsPerRate;
        wordErrRate = blockErrRate/K;
        wordErrRateData.insert(std::make_pair(physicalErrRate, wordErrRate));
        if(physicalErrRate != maxPhysicalErrRate-stepSize) {
            decodingResOutput << "]}},";
        }else{
            decodingResOutput << "]}}";
        }
    }
    decodingResOutput << "}";
    json dataj = wordErrRate;
    rawDataOutput << dataj.dump(2U);
    decodingResOutput.close();
    rawDataOutput.close();
}
