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
 * Original implementation of the generalized UF decoder for QLDPC codes using Gaussian elimination
 * @param syndrome
 */
void OriginalUFD::decode(const std::vector<bool>& syndrome) {
    const auto                                   decodingTimeBegin = std::chrono::high_resolution_clock::now();
    std::unordered_set<std::size_t>              components;
    std::vector<std::unordered_set<std::size_t>> invalidComponents;
    std::unordered_set<std::size_t>              syndr;
    for (std::size_t i = 0; i < syndrome.size(); i++) {
        if (syndrome.at(i)) {
            syndr.insert(getCode()->getN() + i);
        }
    }

    if (!syndr.empty()) {
        // Init a single component for each syndrome vertex
        for (auto s: syndr) {
            components.insert(s);
        }

        // grow all components (including valid ones) by 1
        while (containsInvalidComponents(components, syndr, invalidComponents) && components.size() < (this->getCode()->Hz->pcm->size()+getCode()->Hz->pcm->front().size())) {
            if (this->growth == GrowthVariant::ALL_COMPONENTS) {
                // to grow all components (including valid ones)
                standardGrowth(components);
            } else if (this->growth == GrowthVariant::INVALID_COMPONENTS) {
                //
            } else if (this->growth == GrowthVariant::SINGLE_SMALLEST) {
                singleClusterSmallestFirstGrowth(components);
            } else if (this->growth == GrowthVariant::SINGLE_RANDOM) {
                singleClusterRandomFirstGrowth(components);
                throw std::invalid_argument("Unsupported growth variant");
            } else if (this->growth == GrowthVariant::SINGLE_QUBIT_RANDOM) {
                //singleQubitRandomFirstGrowth(components, neibrsToAdd);
                throw std::invalid_argument("Unsupported growth variant");
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
 * @param components
 * @param syndrome
 * @return
 */
bool OriginalUFD::containsInvalidComponents(const std::unordered_set<std::size_t>& components, const std::unordered_set<std::size_t>& syndrome,
                                            std::vector<std::unordered_set<std::size_t>>& invalidComps) const {
    auto ccomps = getConnectedComps(components);
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
 * @param component
 * @param syndrome
 * @return
 */
bool OriginalUFD::isValidComponent(const std::unordered_set<std::size_t>& component, const std::unordered_set<std::size_t>& syndrome) const {
    return !getEstimateForComponent(component, syndrome).empty();
}

/**
 * Computes a set of nodes s.t. for each n in the list, all neighbours of n are in the component
 * @param component
 * @return
 */
std::vector<std::size_t> OriginalUFD::computeInteriorBitNodes(const std::unordered_set<std::size_t>& component) const {
    std::vector<std::size_t> res;

    for (const auto idx: component) {
        const auto& nbrs = getCode()->Hz->getNbrs(idx);
        if (std::includes(component.begin(), component.end(), nbrs.begin(), nbrs.end()) && idx < getCode()->getN()) {
            res.emplace_back(idx);
        }
    }
    return res;
}

/**
 * Computes estimate vector x for a component and a syndrome
 * @param component
 * @param syndrome
 * @return
 */
std::unordered_set<std::size_t> OriginalUFD::getEstimateForComponent(const std::unordered_set<std::size_t>& component,
                                                                     const std::unordered_set<std::size_t>& syndrome) const {
    std::unordered_set<std::size_t> res{};

    auto                            intNodes = computeInteriorBitNodes(component);
    if (intNodes.empty()) {
        return std::unordered_set<std::size_t>{};
    }
    gf2Mat            redHz;
    std::size_t       idxCnt = 0;
    gf2Vec            redSyndr(idxCnt);
    std::vector<bool> used(this->getCode()->Hz->pcm->size());

    for (const auto it: component) {
        if (it >= getCode()->getN()) { // is a check node
            if (!used.at(it - getCode()->getN())) {
                redHz.emplace_back(this->getCode()->Hz->pcm->at(it - getCode()->getN()));
                used.at(it - getCode()->getN()) = true;
                if (syndrome.contains(it - getCode()->getN())) {
                    redSyndr.emplace_back(1);
                } else {
                    redSyndr.emplace_back(0);
                }
            }
        } else { // is a bit node
            const auto nbrs = this->getCode()->Hz->getNbrs(it);
            for (auto n: nbrs) { // add neighbouring checks
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
    auto                            estim = Utils::solveSystem(redHz, redSyndr);
    std::unordered_set<std::size_t> estIdx;
    for (std::size_t i = 0; i < estim.size(); i++) {
        if (estim.at(i)) {
            res.insert(i);
        }
    }
    return res;
}

void OriginalUFD::standardGrowth(std::unordered_set<std::size_t>& comps) {
    for (auto currCompIt = comps.begin(); currCompIt != comps.end(); currCompIt++) {
        const auto nbrs = getCode()->Hz->getNbrs(*currCompIt);
        for (auto n: nbrs) {
            comps.insert(n);
        }
    }
}

void OriginalUFD::singleClusterSmallestFirstGrowth(std::unordered_set<std::size_t>& comps) {
    auto                            ccomps = getConnectedComps(comps);
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
        comps.insert(nbrs.begin(), nbrs.end());
    }
}

void OriginalUFD::singleClusterRandomFirstGrowth(std::unordered_set<std::size_t>& comps) {
    auto                            ccomps = getConnectedComps(comps);
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
        comps.insert(nbrs.begin(), nbrs.end());
    }
}

void OriginalUFD::reset() {
    this->result = {};
    this->growth = GrowthVariant::ALL_COMPONENTS;
}
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
std::vector<std::unordered_set<std::size_t>> OriginalUFD::getConnectedComps(const std::unordered_set<std::size_t>& nodes) const {
    std::unordered_set<std::size_t>              visited;
    std::vector<std::unordered_set<std::size_t>> result;

    for (auto c: nodes) {
        if (!visited.contains(c)) {
            visited.insert(c);
            std::unordered_set<std::size_t> ccomp;

            std::queue<std::size_t> stack;
            stack.push(c);
            while (!stack.empty()) {
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
