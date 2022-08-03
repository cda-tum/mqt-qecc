//
// Created by lucas on 21/04/2022.
//

#include "OriginalUFD.hpp"

#include "Decoder.hpp"

#include <cassert>
#include <chrono>
#include <queue>
#include <random>
#include <set>

/**
 * Original implementation of the generalized decoder for QLDPC codes using Gaussian elimination
 * @param syndrome
 */
void OriginalUFD::decode(const std::vector<bool>& syndrome) {
    const auto                                   decodingTimeBegin = std::chrono::high_resolution_clock::now();
    std::unordered_set<std::size_t>              components; // used to store vertex indices in E set
    std::vector<std::unordered_set<std::size_t>> invalidComponents;
    std::unordered_set<std::size_t>              syndr; // vertex indices of syndrome nodes
    for (std::size_t i = 0; i < syndrome.size(); i++) {
        if (syndrome.at(i)) {
            syndr.insert(getCode()->getN() + i);
        }
    }

    if (!syndr.empty()) {
        // Set set of nodes equal to syndrome E = syndrome
        for (auto s: syndr) {
            components.insert(s);
        }


        while (containsInvalidComponents(components, syndr, invalidComponents) && components.size() < (this->getCode()->Hz->pcm->size()+getCode()->Hz->pcm->front().size())) {
            if (this->growth == GrowthVariant::ALL_COMPONENTS) {
                // // grow all components (including valid ones) by 1
                standardGrowth(components);
            } else if (this->growth == GrowthVariant::INVALID_COMPONENTS) {
                // not implemented yet
            } else if (this->growth == GrowthVariant::SINGLE_SMALLEST) {
                // grow only by neighbours of single smallest cluster
                singleClusterSmallestFirstGrowth(components);
            } else if (this->growth == GrowthVariant::SINGLE_RANDOM) {
                // grow only by neighbours of single random cluster
                singleClusterRandomFirstGrowth(components);
            } else if (this->growth == GrowthVariant::SINGLE_QUBIT_RANDOM) {
                // grow only by neighbours of single qubit
                singleQubitRandomFirstGrowth(components);
            } else {
                throw std::invalid_argument("Unsupported growth variant");
            }
        }
    }

    std::vector<std::set<std::size_t>> estims;
    auto ccomps = getConnectedComps(components);
    for (const auto& comp: ccomps) {
        auto compEstimate = getEstimateForComponent(comp, syndr);
        estims.emplace_back(compEstimate.begin(), compEstimate.end());
    }
    std::set<std::size_t> tmp;
    for (auto& estim: estims) {
        tmp.insert(estim.begin(), estim.end());
    }
    std::vector<std::size_t> res(tmp.begin(), tmp.end());

    const auto decodingTimeEnd = std::chrono::high_resolution_clock::now();
    result.decodingTime        = std::chrono::duration_cast<std::chrono::milliseconds>(decodingTimeEnd - decodingTimeBegin).count();
    result.estimBoolVector     = std::vector<bool>(getCode()->getN());
    for (unsigned long re: res) {
        result.estimBoolVector.at(re) = true;
    }
    result.estimNodeIdxVector = std::move(res);
}

/**
 * Checks if there is a component in the list that is not valid
 * @param nodeSet
 * @param syndrome
 * @return
 */
bool OriginalUFD::containsInvalidComponents(const std::unordered_set<std::size_t>& nodeSet, const std::unordered_set<std::size_t>& syndrome,
                                            std::vector<std::unordered_set<std::size_t>>& invalidComps) const {
    auto ccomps = getConnectedComps(nodeSet);
    return std::any_of(ccomps.begin(), ccomps.end(), [&](const auto& comp) {
        bool res = isValidComponent(comp, syndrome);
        if (!res) {
            invalidComps.emplace_back(comp.begin(), comp.end());
        }
        return !res;
    });
}

/**
 * Checks if a component is valid
 * A component is valid if there is a set of (bit) nodes in its interior whose syndrome is equal to the given syndrome
 * @param nodeSet
 * @param syndrome
 * @return
 */
bool OriginalUFD::isValidComponent(const std::unordered_set<std::size_t>& nodeSet, const std::unordered_set<std::size_t>& syndrome) const {
    return !getEstimateForComponent(nodeSet, syndrome).empty();
}

/**
 * Computes a set of nodes s.t. for each n in the list, all neighbours of n are in the component
 * @param nodeSet
 * @return
 */
std::vector<std::size_t> OriginalUFD::computeInteriorBitNodes(const std::unordered_set<std::size_t>& nodeSet) const {
    std::vector<std::size_t> res;

    for (const auto idx: nodeSet) {
        const auto& nbrs = getCode()->Hz->getNbrs(idx);
        if (std::includes(nodeSet.begin(), nodeSet.end(), nbrs.begin(), nbrs.end()) && idx < getCode()->getN()) {
            res.emplace_back(idx);
        }
    }
    return res;
}

/**
 * Computes estimate vector x for a component and a syndrome. This is done by considering all vertices in Tanner Graph
 * that are in the Interior of the given node set and additionally the neighbours of the bit vertices in the interior.
 * Then, using Gaussian elimination, it is checked whether a solution for the local cluster that is consitent with the syndrome
 * can be found. If so, this local estimate is returned.
 * @param nodeSet
 * @param syndrome
 * @return
 */
std::unordered_set<std::size_t> OriginalUFD::getEstimateForComponent(const std::unordered_set<std::size_t>& nodeSet,
                                                                     const std::unordered_set<std::size_t>& syndrome) const {
    std::unordered_set<std::size_t> res{};

    auto                            intNodes = computeInteriorBitNodes(nodeSet);
    if (intNodes.empty()) {
        return std::unordered_set<std::size_t>{};
    }
    gf2Mat            redHz;
    std::size_t       idxCnt = 0;
    gf2Vec            redSyndr(idxCnt);
    std::vector<bool> used(this->getCode()->Hz->pcm->size());

    for (const auto it: nodeSet) {
        if (it >= getCode()->getN()) { // is a check node
            if (!used.at(it - getCode()->getN())) {
                redHz.emplace_back(this->getCode()->Hz->pcm->at(it - getCode()->getN()));
                used.at(it - getCode()->getN()) = true;
                if (syndrome.contains(it - getCode()->getN())) {
                    redSyndr.emplace_back(1); // If the check node is in the syndrome we need to satisfy check=1
                } else {
                    redSyndr.emplace_back(0);
                }
            }
        } else { // is a bit node
            const auto nbrs = this->getCode()->Hz->getNbrs(it);
            for (auto n: nbrs) { // add neighbouring checks (these are maybe not in the interior but to stay consistent with the syndrome we need to include these in the check)
                if (!used.at(n - getCode()->getN())) {
                    redHz.emplace_back(this->getCode()->Hz->pcm->at(n - getCode()->getN()));
                    if (syndrome.contains(n)) {
                        redSyndr.emplace_back(1);
                    } else {
                        redSyndr.emplace_back(0);
                    }
                    used.at(n - getCode()->getN()) = true;
                }
            }
        }
    }
    auto                            estim = Utils::solveSystem(redHz, redSyndr); // solves the system redHz*x=redSyndr by x to see if a solution can be found
    std::unordered_set<std::size_t> estIdx;
    for (std::size_t i = 0; i < estim.size(); i++) {
        if (estim.at(i)) {
            res.insert(i);
        }
    }
    return res;
}

/**
 * Grows the node set by the neighbours of ALL clusters
 * @param comps
 */
void OriginalUFD::standardGrowth(std::unordered_set<std::size_t>& comps) {
    for (auto currCompIt = comps.begin(); currCompIt != comps.end(); currCompIt++) {
        const auto nbrs = getCode()->Hz->getNbrs(*currCompIt);
        for (auto n: nbrs) {
            comps.insert(n);
        }
    }
}
/**
 * Grows the node set by the neighbours of the single smallest cluster
 * @param nodeSet
 */
void OriginalUFD::singleClusterSmallestFirstGrowth(std::unordered_set<std::size_t>& nodeSet) {
    auto                            ccomps = getConnectedComps(nodeSet);
    std::unordered_set<std::size_t> compNbrs;
    std::unordered_set<std::size_t> smallestComponent;
    std::size_t                     smallestSize;
    for (const auto& cId: ccomps) {
        if (cId.size() < smallestSize) {
            smallestComponent = cId;
            smallestSize      = cId.size();
        }
    }

    for (auto node: smallestComponent) {
        const auto& nbrs = getCode()->Hz->getNbrs(node);
        nodeSet.insert(nbrs.begin(), nbrs.end());
    }
}

/**
 * Grows the node set by the neighbours of a single random cluster
 * @param nodeSet
 */
void OriginalUFD::singleClusterRandomFirstGrowth(std::unordered_set<std::size_t>& nodeSet) {
    auto                            ccomps = getConnectedComps(nodeSet);
    std::unordered_set<std::size_t> chosenComponent;
    std::random_device              rd;
    std::mt19937                    gen(rd());
    std::uniform_int_distribution   d(static_cast<std::size_t>(0U), ccomps.size());
    std::size_t                     chosenIdx = d(gen);
    auto                            it        = ccomps.begin();
    std::advance(it, chosenIdx);
    chosenComponent = *it;

    for (auto node: chosenComponent) {
        const auto& nbrs = getCode()->Hz->getNbrs(node);
        nodeSet.insert(nbrs.begin(), nbrs.end());
    }
}

/**
 * Reset temporaily computed data
 */
void OriginalUFD::reset() {
    this->result = {};
    this->growth = GrowthVariant::ALL_COMPONENTS;
}

/**
 * Grows the node set by the neighbours of a single random qubit
 * @param comps
 */
void OriginalUFD::singleQubitRandomFirstGrowth(std::unordered_set<std::size_t>& comps) {
    auto                            ccomps = getConnectedComps(comps);
    std::unordered_set<std::size_t> compNbrs;
    std::unordered_set<std::size_t> chosenComponent;
    std::random_device              rd;
    std::mt19937                    gen(rd());
    std::uniform_int_distribution   d(static_cast<std::size_t>(0U), ccomps.size());
    std::size_t                     chosenIdx = d(gen);
    auto                            it        = ccomps.begin();
    std::advance(it, chosenIdx);
    chosenComponent = *it;

    const auto& nbrs = getCode()->Hz->getNbrs(*chosenComponent.begin());
    comps.insert(nbrs.begin(), nbrs.end());
}
/**
 * Given a set of nodes (the set of all nodes considered by the algorithm in the Tanner graph), compute the connected components in the Tanner graph
 * @param nodes
 * @return
 */
std::vector<std::unordered_set<std::size_t>> OriginalUFD::getConnectedComps(const std::unordered_set<std::size_t>& nodes) const {
    std::unordered_set<std::size_t>              visited;
    std::vector<std::unordered_set<std::size_t>> result;

    for (auto c: nodes) {
        if (!visited.contains(c)) {
            visited.insert(c);
            std::unordered_set<std::size_t> ccomp;

            std::queue<std::size_t> stack;
            stack.push(c);
            while (!stack.empty()) { // use DFS-like algorithm to compute connected component containing node 'c'
                auto curr = stack.back();
                stack.pop();
                if (!ccomp.contains(curr)) {
                    ccomp.insert(curr);
                    auto nbrs = getCode()->Hz->getNbrs(curr);
                    for (auto n: nbrs) {
                        if (!ccomp.contains(n) && nodes.contains(n)) {
                            stack.push(n);
                            visited.insert(n);
                        }
                    }
                }
            }
            result.emplace_back(ccomp);
        }
    }
    return result;
}
