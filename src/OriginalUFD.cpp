//
// Created by lucas on 21/04/2022.
//

#include "OriginalUFD.hpp"

#include "Decoder.hpp"
#include "TreeNode.hpp"

#include <cassert>
#include <chrono>
#include <queue>
#include <random>
#include <set>

/**
 * Original implementation of the generalized UF decoder for QLDPC codes using Gaussian elimination
 * @param syndrome
 */
void OriginalUFD::decode(std::vector<bool>& syndrome) {
    std::chrono::steady_clock::time_point decodingTimeBegin = std::chrono::steady_clock::now();
    std::vector<std::set<std::size_t>>    components;
    if (!syndrome.empty()) {
        // Init a single component for each syndrome vertex
        for (size_t i = 0; i < syndrome.size(); i++) {
            std::set<std::size_t> comp{};
            if (syndrome.at(i)) {
                comp.insert(code.getN() + i);
            }

            if (!comp.empty()) {
                components.emplace_back(comp);
            }
        }

        // grow all components (including valid ones) by 1
        while (containsInvalidComponents(components, syndrome)) {
            auto                               currCompIt = components.begin();
            std::vector<std::set<std::size_t>> neibrsToAdd{};
            std::set<std::size_t>              compNbrs;

            while (currCompIt != components.end()) {
                for (auto node: *currCompIt) {
                    auto nbrs = code.tannerGraph.getNeighboursIdx(node);
                    compNbrs.insert(nbrs.begin(), nbrs.end());
                }
                neibrsToAdd.emplace_back(compNbrs);
                currCompIt++;
            }
            for (size_t i = 0; i < components.size(); i++) {
                auto nbrs = neibrsToAdd.at(i);
                components.at(i).insert(nbrs.begin(), nbrs.end());
            }
        }
    }
    std::vector<std::set<std::size_t>> corrEstimates;
    for (auto& comp: components) {
        std::set<std::size_t> compEstimate = getEstimateForComponent(comp, syndrome);
        corrEstimates.emplace_back(compEstimate);
    }
    std::set<std::size_t> tmp;
    for (auto& estim: corrEstimates) {
        tmp.insert(estim.begin(), estim.end());
    }
    std::vector<std::size_t> res(tmp.begin(), tmp.end());

    std::chrono::steady_clock::time_point decodingTimeEnd = std::chrono::steady_clock::now();
    result.decodingTime                                   = std::chrono::duration_cast<std::chrono::milliseconds>(decodingTimeBegin - decodingTimeEnd).count();
    result.estimNodeIdxVector                             = res;
    result.estimBoolVector                                = std::vector<bool>(code.getN());
    for (unsigned long re: res) {
        result.estimBoolVector.at(re) = true;
    }
}

/**
 * Checks if there is a component in the list that is not valid
 * @param components
 * @param syndrome
 * @return
 */
bool OriginalUFD::containsInvalidComponents(std::vector<std::set<std::size_t>>& components, const std::vector<bool>& syndrome) {
    auto it = components.begin();
    while (it != components.end()) {
        if (!isValidComponent(*it, syndrome)) {
            return true;
        }
        it++;
    }
    return false;
}

/**
 * Checks if a component is valid
 * A component is valid if there is a set of (bit) nodes in its interior whose syndrome is equal to the given syndrome
 * @param component
 * @param syndrome
 * @return
 */
bool OriginalUFD::isValidComponent(std::set<std::size_t>& component, const std::vector<bool>& syndrome) {
    std::vector<bool> syndr;
    for (bool i: syndrome) {
        if (i) {
            syndr.push_back(true);
        } else {
            syndr.push_back(false);
        }
    }
    auto sol = Utils::solveSystem(code.Hz.pcm, syndr);
    if (sol.empty()) {
        return false;
    } else {
        std::vector<std::size_t> estimNodes;
        for (size_t i = 0; i < sol.size(); i++) {
            if (sol.at(i)) {
                estimNodes.emplace_back(i);
            }
        }
        auto interior = computeInteriorBitNodes(component);
        return std::includes(interior.begin(), interior.end(), estimNodes.begin(), estimNodes.end());
    }
}

/**
 * Computes a set of nodes s.t. for each n in the list, all neighbours of n are in the component
 * @param component
 * @return
 */
std::vector<std::size_t> OriginalUFD::computeInteriorBitNodes(std::set<std::size_t>& component) {
    std::vector<std::size_t> res;

    auto cIt = component.begin();
    while (cIt != component.end()) {
        auto nbrs = code.tannerGraph.getNeighboursIdx(*cIt);
        if (std::includes(component.begin(), component.end(), nbrs.begin(), nbrs.end()) && *cIt <= code.getN()) {
            res.emplace_back(*cIt);
        }
        cIt++;
    }
    return res;
}

/**
 * Computes estimate vector x for a component and a syndrome
 * @param set
 * @param syndrome
 * @return
 */
std::set<std::size_t> OriginalUFD::getEstimateForComponent(std::set<std::size_t>& set, const std::vector<bool>& syndrome) {
    std::vector<bool>     syndr;
    std::set<std::size_t> res;
    for (bool i: syndrome) {
        if (i) {
            syndr.push_back(true);
        } else {
            syndr.push_back(false);
        }
    }
    auto estim = Utils::solveSystem(code.Hz.pcm, syndr);
    for (auto&& i: estim) {
        if (i) {
            res.insert(i);
        }
    }
    return res;
}
