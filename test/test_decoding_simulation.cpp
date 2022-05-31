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

std::vector<bool> sampleXError(const std::size_t n) {
    std::vector<bool>               result(n);
    std::random_device              rd;
    std::mt19937                    gen(rd());
    std::uniform_int_distribution<> distr(0, 1);

    for (std::size_t i = 0; i < n; i++) {
        result.emplace_back(distr(gen));
    }

    return result;
}
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
/*
TEST(BenchmarkSimulation, SteaneCodeDecoding) {
    std::string        outFilePath = "path";
    auto               t           = std::time(nullptr);
    auto               tm          = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto          timestamp = oss.str();
    std::ofstream outfile(outFilePath + "EC-" + timestamp + ".json");
    outfile << "{ \"runs\" : [";

    SteaneXCode                         code{};
    Decoder                             decoder{code};
    auto                                error = sampleXError(code.getN());
    auto                                syndr = code.getSyndrome(error);
    std::set<std::shared_ptr<TreeNode>> syndrComponents;
    for (size_t i = 0; i < syndr.size(); i++) {
        if (syndr.at(i)) {
            auto snode     = code.tannerGraph.adjListNodes.at(i + code.getN()).at(0);
            snode->isCheck = true;
            snode->checkVertices.insert(snode->vertexIdx);
            syndrComponents.insert(snode);
        }
    }
    decoder.decode(syndrComponents);
    auto              decodingResult = decoder.result;
    auto              estim          = decodingResult.estimNodeIdxVector;
    std::vector<bool> residualErr    = decodingResult.estimBoolVector;

    for (std::size_t i = 0; i < residualErr.size(); i++) {
        residualErr.at(i) = residualErr.at(i) ^ error.at(i);
    }
    auto success = code.checkStabilizer(residualErr);
    if (success) {
        std::cout << "Decoding successful, found residualErr up to stabilizer: " << std::endl;
        std::cout << residualErr << std::endl;
        std::cout << "Elapsed time: " << decodingResult.decodingTime << "ms" << std::endl;
    } else {
        decodingResult.status = FAILURE;
        std::cout << "Decoding not successful, introduced logical opertor" << std::endl;
    }
    outfile << decodingResult.to_json().dump(2U);
}
*/