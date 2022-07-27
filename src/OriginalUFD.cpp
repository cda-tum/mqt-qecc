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
    const auto                         decodingTimeBegin = std::chrono::high_resolution_clock::now();
    std::vector<std::set<std::size_t>> components;
    std::vector<std::set<std::size_t>> invalidComponents;
    if (!syndrome.empty()) {
        // Init a single component for each syndrome vertex
        for (std::size_t i = 0U; i < syndrome.size(); i++) {
            if (syndrome.at(i)) {
                components.emplace_back(std::set{getCode()->getN() + i});
            }
        }

        // grow all components (including valid ones) by 1
        while (containsInvalidComponents(components, syndrome, invalidComponents)) {
            std::vector<std::set<std::size_t>> neibrsToAdd{};
            if (this->growth == GrowthVariant::ALL_COMPONENTS) {
                // to grow all components (including valid ones)
                standardGrowth(components, neibrsToAdd);
            } else if (this->growth == GrowthVariant::INVALID_COMPONENTS) {
                standardGrowth(invalidComponents, neibrsToAdd);
            } else if (this->growth == GrowthVariant::SINGLE_SMALLEST) {
                singleClusterSmallestFirstGrowth(components, neibrsToAdd);
            } else if (this->growth == GrowthVariant::SINGLE_RANDOM) {
                singleClusterRandomFirstGrowth(components, neibrsToAdd);
            } else if (this->growth == GrowthVariant::SINGLE_QUBIT_RANDOM) {
                singleQubitRandomFirstGrowth(components, neibrsToAdd);
            } else {
                throw std::invalid_argument("Unsupported growth variant");
            }

            for (std::size_t i = 0; i < components.size(); i++) {
                const auto& nbrs = neibrsToAdd.at(i);
                components.at(i).insert(nbrs.begin(), nbrs.end());
            }
            neibrsToAdd.clear();
        }
    }

    std::vector<std::set<std::size_t>> corrEstimates;
    for (const auto& comp: components) {
        std::set<std::size_t> compEstimate = getEstimateForComponent(comp, syndrome);
        corrEstimates.emplace_back(compEstimate);
    }
    std::set<std::size_t> tmp;
    for (auto& estim: corrEstimates) {
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
bool OriginalUFD::containsInvalidComponents(const std::vector<std::set<std::size_t>>& components, const std::vector<bool>& syndrome,
                                            std::vector<std::set<std::size_t>>& invalidComps) const {
    return std::any_of(components.begin(), components.end(), [&](const auto& comp) {
        auto res = isValidComponent(comp, syndrome);
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
bool OriginalUFD::isValidComponent(const std::set<std::size_t>& component, const std::vector<bool>& syndrome) const {
    return !getEstimateForComponent(component, syndrome).empty();
}

/**
 * Computes a set of nodes s.t. for each n in the list, all neighbours of n are in the component
 * @param component
 * @return
 */
std::vector<std::size_t> OriginalUFD::computeInteriorBitNodes(const std::set<std::size_t>& component) const {
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
std::set<std::size_t> OriginalUFD::getEstimateForComponent(const std::set<std::size_t>& component, const std::vector<bool>& syndrome) const {
    std::set<std::size_t> res;
    auto                  intNodes = computeInteriorBitNodes(component);
    if (intNodes.empty()) {
        return std::set<std::size_t>{};
    }
    if (auto estim = Utils::solveSystem(*getCode()->Hz->pcm, syndrome); estim.empty()) {
        return res;
    } else {
        std::set<std::size_t> estIdx;
        for (std::size_t i = 0; i < estim.size(); i++) {
            if (estim.at(i)) {
                estIdx.insert(i);
            }
        }
        if (std::includes(intNodes.begin(), intNodes.end(), estIdx.begin(), estIdx.end())) {
            return estIdx;
        }
    }

    return res;
}

void OriginalUFD::standardGrowth(const std::vector<std::set<std::size_t>>& comps, std::vector<std::set<std::size_t>>& neibrsToAdd) {
    std::set<std::size_t> compNbrs;

    auto currCompIt = comps.begin();
    while (currCompIt != comps.end()) {
        for (auto node: *currCompIt) {
            const auto& nbrs = getCode()->Hz->getNbrs(node);
            compNbrs.insert(nbrs.begin(), nbrs.end());
        }
        neibrsToAdd.emplace_back(compNbrs);
        currCompIt++;
    }
}

void OriginalUFD::singleClusterSmallestFirstGrowth(const std::vector<std::set<std::size_t>>& comps, std::vector<std::set<std::size_t>>& neibrsToAdd) {
    std::set<std::size_t> compNbrs;
    std::set<std::size_t> smallestComponent;
    std::size_t           smallestSize = SIZE_MAX;
    for (const auto& cId: comps) {
        if (cId.size() < smallestSize) {
            smallestComponent = cId;
        }
    }

    for (auto node: smallestComponent) {
        const auto& nbrs = getCode()->Hz->getNbrs(node);
        compNbrs.insert(nbrs.begin(), nbrs.end());
    }
    neibrsToAdd.emplace_back(compNbrs);
}

void OriginalUFD::singleClusterRandomFirstGrowth(const std::vector<std::set<std::size_t>>& comps, std::vector<std::set<std::size_t>>& neibrsToAdd) {
    std::set<std::size_t>         compNbrs;
    std::set<std::size_t>         chosenComponent;
    std::random_device            rd;
    std::mt19937                  gen(rd());
    std::uniform_int_distribution d(static_cast<std::size_t>(0U), comps.size());
    std::size_t                   chosenIdx = d(gen);
    auto                          it        = comps.begin();
    std::advance(it, chosenIdx);
    chosenComponent = *it;

    for (auto node: chosenComponent) {
        const auto& nbrs = getCode()->Hz->getNbrs(node);
        compNbrs.insert(nbrs.begin(), nbrs.end());
    }
    neibrsToAdd.emplace_back(compNbrs);
}

void OriginalUFD::reset() {
    this->result = {};
    this->growth = GrowthVariant::ALL_COMPONENTS;
}
void OriginalUFD::singleQubitRandomFirstGrowth(const std::vector<std::set<std::size_t>>& comps, std::vector<std::set<std::size_t>>& neibrsToAdd) {
    std::set<std::size_t>         compNbrs;
    std::set<std::size_t>         chosenComponent;
    std::random_device            rd;
    std::mt19937                  gen(rd());
    std::uniform_int_distribution d(static_cast<std::size_t>(0U), comps.size());
    std::size_t                   chosenIdx = d(gen);
    auto                          it        = comps.begin();
    std::advance(it, chosenIdx);
    chosenComponent = *it;

    const auto& nbrs = getCode()->Hz->getNbrs(*chosenComponent.begin());
    compNbrs.insert(nbrs.begin(), nbrs.end());

    neibrsToAdd.emplace_back(compNbrs);
}
