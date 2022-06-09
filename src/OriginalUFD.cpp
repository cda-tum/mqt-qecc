//
// Created by lucas on 21/04/2022.
//

#include "Decoder.hpp"
#include "OriginalUFD.hpp"
#include "TreeNode.hpp"

#include <cassert>
#include <chrono>
#include <queue>
#include <random>
#include <set>

void OriginalUFD::decode(std::vector<bool>& syndrome) {
    std::chrono::steady_clock::time_point decodingTimeBegin = std::chrono::steady_clock::now();
    std::vector<std::size_t>              res;
    auto components = syndrome;

   // grow

    std::chrono::steady_clock::time_point decodingTimeEnd = std::chrono::steady_clock::now();
    result.decodingTime                                   = std::chrono::duration_cast<std::chrono::milliseconds>(decodingTimeBegin - decodingTimeEnd).count();
    result.estimNodeIdxVector                             = res;
    result.estimBoolVector                                = std::vector<bool>(code.getN());
    for (unsigned long re: res) {
        result.estimBoolVector.at(re) = true;
    }
}

bool OriginalUFD::containsInvalidComponents(const std::set<std::shared_ptr<TreeNode>>& components){
    auto it = components.begin();
    while(it != components.end()){
        if(!isValidComponent(*it)){
            return true;
        }
        it++;
    }
    return false;
}

void OriginalUFD::standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                    std::map<std::size_t, bool>& presentMap, const std::set<std::shared_ptr<TreeNode>>& components) {
    for (auto& component: components) {
        presentMap.insert(std::make_pair(component->vertexIdx, true));
        assert(component->parent == nullptr); // at this point we can assume that component represents root of the component
        auto bndryNodes = component->boundaryVertices;

        for (const auto& bndryNode: bndryNodes) {
            auto nbrs = code.tannerGraph.getNeighboursIdx(bndryNode);
            for (auto& nbr: nbrs) {
                fusionEdges.emplace_back(std::pair(bndryNode, nbr));
            }
        }
    }
}

bool OriginalUFD::isValidComponent(const std::shared_ptr<TreeNode>& component) {
    return false; // todo
}
