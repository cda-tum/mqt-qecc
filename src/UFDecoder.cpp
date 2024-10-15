#include "UFDecoder.hpp"

#include "Decoder.hpp"
#include "GF2.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <queue>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

/**
 * Original implementation of the generalized decoder for QLDPC codes using Gaussian elimination
 * @param syndrome
 */
void UFDecoder::decode(const gf2Vec& syndrome) {
    if (syndrome.size() > this->getCode()->gethZ()->pcm->size()) {
        std::vector<bool> xSyndr;
        std::vector<bool> zSyndr;
        auto              mid = syndrome.begin() + (static_cast<std::int64_t>(syndrome.size()) / 2U);
        std::move(syndrome.begin(), mid, std::back_inserter(xSyndr));
        std::move(mid, syndrome.end(), std::back_inserter(zSyndr));
        doDecode(xSyndr, this->getCode()->gethZ());
        auto xres = this->result;
        this->reset();
        doDecode(zSyndr, this->getCode()->gethX());
        this->result.decodingTime += xres.decodingTime;
        std::move(xres.estimBoolVector.begin(), xres.estimBoolVector.end(),
                  std::back_inserter(this->result.estimBoolVector));
        std::move(xres.estimNodeIdxVector.begin(), xres.estimNodeIdxVector.end(),
                  std::back_inserter(this->result.estimNodeIdxVector));
    } else {
        this->doDecode(syndrome, getCode()->gethZ()); // X errs per default if single sided
    }
}

void UFDecoder::doDecode(const std::vector<bool>& syndrome, const std::unique_ptr<ParityCheckMatrix>& pcm) {
    const auto                         decodingTimeBegin = std::chrono::high_resolution_clock::now();
    std::set<std::size_t>              components; // used to store vertex indices in E set
    std::vector<std::set<std::size_t>> invalidComponents;
    std::set<std::size_t>              syndr; // vertex indices of syndrome nodes
    for (std::size_t i = 0; i < syndrome.size(); i++) {
        if (syndrome.at(i)) {
            syndr.insert(getCode()->getN() + i);
        }
    }

    if (!syndr.empty()) {
        // Set set of nodes equal to syndrome E = syndrome
        for (auto s : syndr) {
            components.insert(s);
        }

        while (containsInvalidComponents(components, syndr, invalidComponents, pcm) &&
               components.size() < (pcm->pcm->size() + pcm->pcm->front().size())) {
            if (this->growth == GrowthVariant::AllComponents) {
                // // grow all components (including valid ones) by 1
                standardGrowth(components);
            } else if (this->growth == GrowthVariant::InvalidComponents) {
                // not implemented yet
            } else if (this->growth == GrowthVariant::SingleSmallest) {
                // grow only by neighbours of single smallest cluster
                singleClusterSmallestFirstGrowth(components);
            } else if (this->growth == GrowthVariant::SingleRandom) {
                // grow only by neighbours of single random cluster
                singleClusterRandomFirstGrowth(components);
            } else if (this->growth == GrowthVariant::SingleQubitRandom) {
                // grow only by neighbours of single qubit
                singleQubitRandomFirstGrowth(components);
            } else {
                throw std::invalid_argument("Unsupported growth variant");
            }
        }
    }

    std::vector<std::set<std::size_t>> estims;
    auto                               ccomps = getConnectedComps(components);
    for (const auto& comp : ccomps) {
        auto compEstimate = getEstimateForComponent(comp, syndr, pcm);
        estims.emplace_back(compEstimate.begin(), compEstimate.end());
    }
    std::set<std::size_t> tmp;
    for (auto& estim : estims) {
        tmp.insert(estim.begin(), estim.end());
    }
    std::vector<std::size_t> res(tmp.begin(), tmp.end());

    const auto decodingTimeEnd = std::chrono::high_resolution_clock::now();
    result.decodingTime        = static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                                                           decodingTimeEnd - decodingTimeBegin)
                                                                  .count());
    result.estimBoolVector     = std::vector<bool>(getCode()->getN());
    for (auto re : res) {
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
bool UFDecoder::containsInvalidComponents(const std::set<std::size_t>&              nodeSet,
                                          const std::set<std::size_t>&              syndrome,
                                          std::vector<std::set<std::size_t>>&       invalidComps,
                                          const std::unique_ptr<ParityCheckMatrix>& pcm) const {
    auto ccomps = getConnectedComps(nodeSet);
    return std::any_of(ccomps.begin(), ccomps.end(), [&](const auto& comp) {
        bool const res = isValidComponent(comp, syndrome, pcm);
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
bool UFDecoder::isValidComponent(const std::set<std::size_t>&              nodeSet,
                                 const std::set<std::size_t>&              syndrome,
                                 const std::unique_ptr<ParityCheckMatrix>& pcm) const {
    return !getEstimateForComponent(nodeSet, syndrome, pcm).empty();
}

/**
 * Computes a set of nodes s.t. for each n in the list, all neighbours of n are in the component
 * @param nodeSet
 * @return
 */
std::vector<std::size_t> UFDecoder::computeInteriorBitNodes(const std::set<std::size_t>& nodeSet) const {
    std::vector<std::size_t> res;

    for (const auto idx : nodeSet) {
        const auto& nbrs = getCode()->gethZ()->getNbrs(idx);
        if (idx < getCode()->getN() && std::includes(nodeSet.begin(), nodeSet.end(), nbrs.begin(), nbrs.end())) {
            res.emplace_back(idx);
        }
    }
    return res;
}

/**
 * Computes estimate vector x for a component and a syndrome. This is done by considering all vertices in Tanner Graph
 * that are in the Interior of the given node set and additionally the neighbours of the bit vertices in the interior.
 * Then, using Gaussian elimination, it is checked whether a solution for the local cluster that is consistent with the syndrome
 * can be found. If so, this local estimate is returned.
 * @param nodeSet
 * @param syndrome
 * @return
 */
std::set<std::size_t> UFDecoder::getEstimateForComponent(const std::set<std::size_t>&              nodeSet,
                                                         const std::set<std::size_t>&              syndrome,
                                                         const std::unique_ptr<ParityCheckMatrix>& pcm) const {
    if (computeInteriorBitNodes(nodeSet).empty()) {
        return {};
    }

    gf2Mat            redHz;
    gf2Vec            redSyndr(0);
    std::vector<bool> used(pcm->pcm->size());

    for (const auto it : nodeSet) {
        if (it >= getCode()->getN()) { // is a check node
            if (!used.at(it - getCode()->getN())) {
                redHz.emplace_back(pcm->pcm->at(it - getCode()->getN()));
                used.at(it - getCode()->getN()) = true;
                if (syndrome.find(it - getCode()->getN()) != syndrome.end()) {
                    redSyndr.emplace_back(1); // If the check node is in the syndrome we need to satisfy check=1
                } else {
                    redSyndr.emplace_back(0);
                }
            }
        } else { // is a bit node
            const auto nbrs = this->getCode()->gethZ()->getNbrs(it);
            for (auto n : nbrs) { // add neighbouring checks (these are maybe not in the interior but to stay consistent with the syndrome we need to include these in the check)
                if (!used.at(n - getCode()->getN())) {
                    redHz.emplace_back(pcm->pcm->at(n - getCode()->getN()));
                    if (syndrome.find(n) != syndrome.end()) {
                        redSyndr.emplace_back(1);
                    } else {
                        redSyndr.emplace_back(0);
                    }
                    used.at(n - getCode()->getN()) = true;
                }
            }
        }
    }
    auto                  redHzCsc = Utils::toCsc(redHz);
    std::vector<uint64_t> redSyndInt(redSyndr.size());
    for (std::size_t i = 0; i < redSyndr.size(); i++) {
        redSyndInt.at(i) = redSyndr.at(i) ? 1 : 0;
    }
    auto pluDec = PluDecomposition(redHz.size(), redHz.at(0).size(), redHzCsc);
    auto estim  = pluDec.luSolve(redSyndInt); // solves the system redHz*x=redSyndr by x to see if a solution can be found

    std::set<std::size_t> res;
    for (std::size_t i = 0; i < estim.size(); i++) {
        if (estim.at(i) != 0U) {
            res.emplace(i);
        }
    }
    return res;
}

/**
 * Grows the node set by the neighbours of ALL clusters
 * @param comps
 */
void UFDecoder::standardGrowth(std::set<std::size_t>& comps) {
    for (auto currCompIt = comps.begin(); currCompIt != comps.end(); currCompIt++) {
        const auto nbrs = getCode()->gethZ()->getNbrs(*currCompIt);
        for (auto n : nbrs) {
            comps.insert(n);
        }
    }
}

/**
 * Grows the node set by the neighbours of the single smallest cluster
 * @param nodeSet
 */
void UFDecoder::singleClusterSmallestFirstGrowth(std::set<std::size_t>& nodeSet) {
    auto                  ccomps = getConnectedComps(nodeSet);
    std::set<std::size_t> smallestComponent;
    std::size_t           smallestSize = SIZE_MAX;
    for (const auto& cId : ccomps) {
        if (cId.size() < smallestSize) {
            smallestComponent = cId;
            smallestSize      = cId.size();
        }
    }

    for (auto node : smallestComponent) {
        const auto& nbrs = getCode()->gethZ()->getNbrs(node);
        nodeSet.insert(nbrs.begin(), nbrs.end());
    }
}

/**
 * Grows the node set by the neighbours of a single random cluster
 * @param nodeSet
 */
void UFDecoder::singleClusterRandomFirstGrowth(std::set<std::size_t>& nodeSet) {
    auto                          ccomps = getConnectedComps(nodeSet);
    std::set<std::size_t>         chosenComponent;
    std::random_device            rd;
    std::mt19937                  gen(rd());
    std::uniform_int_distribution d(static_cast<std::size_t>(0U), ccomps.size() - 1);
    const std::size_t             chosenIdx = d(gen);
    auto                          it        = ccomps.begin();
    std::advance(it, chosenIdx);
    chosenComponent = *it;

    for (auto node : chosenComponent) {
        const auto& nbrs = getCode()->gethZ()->getNbrs(node);
        nodeSet.insert(nbrs.begin(), nbrs.end());
    }
}

/**
 * Reset temporarily computed data
 */
void UFDecoder::reset() {
    this->result = {};
    this->growth = GrowthVariant::AllComponents;
}

/**
 * Grows the node set by the neighbours of a single random qubit
 * @param comps
 */
void UFDecoder::singleQubitRandomFirstGrowth(std::set<std::size_t>& comps) {
    auto                          ccomps = getConnectedComps(comps);
    std::set<std::size_t>         chosenComponent;
    std::random_device            rd;
    std::mt19937                  gen(rd());
    std::uniform_int_distribution d(static_cast<std::size_t>(0U), ccomps.size());
    const std::size_t             chosenIdx = d(gen);
    auto                          it        = ccomps.begin();
    std::advance(it, chosenIdx);
    chosenComponent = *it;

    const auto& nbrs = getCode()->gethZ()->getNbrs(*chosenComponent.begin());
    comps.insert(nbrs.begin(), nbrs.end());
}

/**
 * Given a set of nodes (the set of all nodes considered by the algorithm in the Tanner graph), compute the connected components in the Tanner graph
 * @param nodes
 * @return
 */
std::vector<std::set<std::size_t>>
UFDecoder::getConnectedComps(const std::set<std::size_t>& nodes) const {
    std::set<std::size_t>              visited;
    std::vector<std::set<std::size_t>> res;

    for (auto c : nodes) {
        if (visited.find(c) == visited.end()) {
            visited.insert(c);
            std::set<std::size_t> ccomp;

            std::queue<std::size_t> stack;
            stack.push(c);
            while (!stack.empty()) { // use DFS-like algorithm to compute connected component containing node 'c'
                auto curr = stack.back();
                stack.pop();
                if (ccomp.find(curr) == ccomp.end()) {
                    ccomp.insert(curr);
                    auto nbrs = getCode()->gethZ()->getNbrs(curr);
                    for (auto n : nbrs) {
                        if (ccomp.find(n) == ccomp.end() && nodes.find(n) != nodes.end()) {
                            stack.push(n);
                            visited.insert(n);
                        }
                    }
                }
            }
            res.emplace_back(ccomp);
        }
    }
    return res;
}
