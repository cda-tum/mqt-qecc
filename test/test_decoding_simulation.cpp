//
// Created by lucas on 21/04/2022.
//
#include "Codes.hpp"
#include "Decoder.hpp"

#include <filesystem>
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

TEST(UnionFindSimulation, SteaneCodeDecoding) {
    SteaneXCode code{};
    Decoder    decoder{code};
    std::cout << "code: " << std::endl<< code << std::endl;
    auto err = dummySampler(code.getN());
    std::cout << "error: "  << err << std::endl;
    auto syndr = code.getSyndrome(err);
    std::cout << "syndr: " << syndr << std::endl;
    std::set<std::shared_ptr<TreeNode>> syndrComponents;
    for (size_t i = 0; i < syndr.size(); i++) {
        if(syndr.at(i)){
            auto snode = code.tannerGraph.adjListNodes.at(i+code.getN()).at(0);
            snode->isCheck = true;
            snode->checkVertices.insert(snode->vertexIdx);
            syndrComponents.insert(snode);
        }
    }
    std::cout << "s component " << syndrComponents << std::endl;
    auto estim = decoder.decode(syndrComponents);
    if(estim.empty()){
        std::cout << "Decoding failure" << std::endl;
    }
    for (auto & x: estim) {
        std::cout << "estim: " << x << std::endl;
    }
    assert(!estim.empty());
}
