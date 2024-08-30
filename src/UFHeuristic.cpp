#include "UFHeuristic.hpp"

#include "Decoder.hpp"
#include "QeccException.hpp"
#include "TreeNode.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <queue>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
/**
 * returns list of tree node (in UF data structure) representations for syndrome
 * @param code
 * @param syndrome
 * @return
 */
std::unordered_set<std::size_t> UFHeuristic::computeInitTreeComponents(const gf2Vec& syndrome) {
    std::unordered_set<std::size_t> res{};
    for (std::size_t i = 0; i < syndrome.size(); i++) {
        if (syndrome.at(i)) {
            const auto idx       = i + getCode()->getN();
            auto       syndrNode = std::make_unique<TreeNode>(idx);
            syndrNode->isCheck   = true;
            syndrNode->checkVertices.emplace_back(syndrNode->vertexIdx);
            nodeMap.try_emplace(syndrNode->vertexIdx, std::move(syndrNode));
            res.insert(idx);
        }
    }
    return res;
}

/**
 * Main part of the heuristic. Uses Union-Find datastructure for efficient cluster growth and validtiy check
 * @param syndrome
 */
void UFHeuristic::decode(const gf2Vec& syndrome) {
    if (syndrome.size() > this->getCode()->gethZ()->pcm->size()) {
        std::vector<bool> xSyndr;
        std::vector<bool> zSyndr;
        auto              mid = syndrome.begin() + (static_cast<std::int64_t>(std::size(syndrome)) / 2U);
        std::move(syndrome.begin(), mid, std::back_inserter(xSyndr));
        std::move(mid, syndrome.end(), std::back_inserter(zSyndr));
        doDecoding(xSyndr, this->getCode()->gethZ());
        auto xres = this->result;
        this->reset();
        doDecoding(zSyndr, this->getCode()->gethZ());
        this->result.decodingTime += xres.decodingTime;
        std::move(xres.estimBoolVector.begin(), xres.estimBoolVector.end(), std::back_inserter(this->result.estimBoolVector));
        std::move(xres.estimNodeIdxVector.begin(), xres.estimNodeIdxVector.end(), std::back_inserter(this->result.estimNodeIdxVector));
    } else {
        this->doDecoding(syndrome, getCode()->gethZ()); // X errs per default if single sided
    }
}
/**
 * Main part of the heuristic. Uses Union-Find datastructure for efficient cluster growth and validtiy check
 * @param syndrome
 */
void UFHeuristic::doDecoding(const gf2Vec& syndrome, const std::unique_ptr<ParityCheckMatrix>& pcm) {
    auto                     decodingTimeBegin = std::chrono::high_resolution_clock::now();
    std::vector<std::size_t> res;
    if (!syndrome.empty() && !std::all_of(syndrome.begin(), syndrome.end(), [](bool val) { return !val; })) {
        auto                            syndrComponents   = computeInitTreeComponents(syndrome);
        auto                            invalidComponents = syndrComponents;
        std::unordered_set<std::size_t> erasure;
        while (!invalidComponents.empty() && invalidComponents.size() < (this->getCode()->gethZ()->pcm->size() + getCode()->gethZ()->pcm->front().size())) {
            // Step 1 growth
            std::vector<std::pair<std::size_t, std::size_t>> fusionEdges;
            std::unordered_map<std::size_t, bool>            presentMap; // for step 4

            if (this->growth == GrowthVariant::AllComponents) {
                // to grow all components (including valid ones)
                for (auto e : erasure) {
                    invalidComponents.insert(e);
                }
                standardGrowth(fusionEdges, presentMap, invalidComponents, pcm);
            } else if (this->growth == GrowthVariant::InvalidComponents) {
                standardGrowth(fusionEdges, presentMap, invalidComponents, pcm);
            } else if (this->growth == GrowthVariant::SingleSmallest) {
                singleClusterSmallestFirstGrowth(fusionEdges, presentMap, invalidComponents, pcm);
            } else if (this->growth == GrowthVariant::SingleRandom) {
                singleClusterRandomFirstGrowth(fusionEdges, presentMap, invalidComponents, pcm);
            } else {
                throw std::invalid_argument("Unsupported growth variant");
            }
            // Fuse clusters that grew together
            auto eIt = fusionEdges.begin();
            while (eIt != fusionEdges.end()) {
                auto* n1    = getNodeFromIdx(eIt->first);
                auto* n2    = getNodeFromIdx(eIt->second);
                auto* root1 = TreeNode::Find(n1);
                auto* root2 = TreeNode::Find(n2);
                // compares vertexIdx only
                if (root1->vertexIdx == root2->vertexIdx) {
                    eIt = fusionEdges.erase(eIt); // passes eIt to erase
                } else {
                    auto s1 = root1->clusterSize; // sizes before union needed
                    auto s2 = root2->clusterSize;
                    TreeNode::Union(root1, root2);

                    if (s1 <= s2) { // Step 3 fusion of boundary lists
                        for (const auto& boundaryVertex : root1->boundaryVertices) {
                            root2->boundaryVertices.insert(boundaryVertex);
                        }
                        root1->boundaryVertices.clear();
                    } else {
                        for (const auto& boundaryVertex : root2->boundaryVertices) {
                            root1->boundaryVertices.insert(boundaryVertex);
                        }
                        root2->boundaryVertices.clear();
                    }
                    ++eIt;
                }
            }
            // Replace nodes in list by their roots avoiding duplicates
            std::vector<std::size_t> toAdd;
            auto                     idxIt = invalidComponents.begin();
            while (idxIt != invalidComponents.end()) {
                auto*       elem = getNodeFromIdx(*idxIt);
                const auto& root = TreeNode::Find(elem);
                if (elem->vertexIdx != root->vertexIdx && presentMap.find(root->vertexIdx) != presentMap.end()) {
                    // root already in component list, no replacement necessary
                    idxIt = invalidComponents.erase(idxIt);
                } else {
                    // root of component not yet in list, replace node by its root in components
                    idxIt = invalidComponents.erase(idxIt);
                    toAdd.emplace_back(root->vertexIdx);
                }
            }
            for (auto c : toAdd) {
                invalidComponents.insert(c);
            }

            // Update Boundary Lists: remove vertices that are not in boundary anymore
            for (const auto& compId : invalidComponents) {
                const auto& compNode = getNodeFromIdx(compId);
                auto        iter     = compNode->boundaryVertices.begin();
                while (iter != compNode->boundaryVertices.end()) {
                    const auto& nbrs     = pcm->getNbrs(*iter);
                    auto*       currNode = getNodeFromIdx(*iter);
                    const auto& currRoot = TreeNode::Find(currNode);
                    for (const auto& nbr : nbrs) {
                        auto* node = getNodeFromIdx(nbr);
                        if (auto* const nbrRoot = TreeNode::Find(node); currRoot->vertexIdx != nbrRoot->vertexIdx) {
                            // if we find one neighbour that is not in the same component the currNode is in the boundary
                            iter++;
                            break;
                        }
                        if (nbr == nbrs.back()) {
                            // if we have checked all neighbours and found none in another component remove from boundary list
                            // toRemoveList.emplace_back(*iter);
                            iter = compNode->boundaryVertices.erase(iter);
                            break;
                        }
                    }
                }
            }
            extractValidComponents(invalidComponents, erasure, pcm);
        }
        res = erasureDecoder(erasure, syndrComponents, pcm);
    }
    auto decodingTimeEnd   = std::chrono::high_resolution_clock::now();
    result                 = DecodingResult();
    result.decodingTime    = static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::milliseconds>(decodingTimeEnd - decodingTimeBegin).count());
    result.estimBoolVector = gf2Vec(getCode()->getN());
    for (const auto& re : res) {
        result.estimBoolVector.at(re) = true;
        result.estimNodeIdxVector.emplace_back(re);
    }
}

void UFHeuristic::standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                 std::unordered_map<std::size_t, bool>&            presentMap,
                                 const std::unordered_set<std::size_t>&            components,
                                 const std::unique_ptr<ParityCheckMatrix>&         pcm) {
    for (const auto& compId : components) {
        const auto& compNode = getNodeFromIdx(compId);
        presentMap.try_emplace(compNode->vertexIdx, true);
        const auto& bndryNodes = compNode->boundaryVertices;

        for (const auto& bndryNode : bndryNodes) {
            const auto& nbrs = pcm->getNbrs(bndryNode);
            for (const auto& nbr : nbrs) {
                fusionEdges.emplace_back(bndryNode, nbr);
            }
        }
    }
}

void UFHeuristic::singleClusterSmallestFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                   std::unordered_map<std::size_t, bool>& presentMap, const std::unordered_set<std::size_t>& components, const std::unique_ptr<ParityCheckMatrix>& pcm) {
    std::size_t smallestComponent = *components.begin();
    std::size_t smallestSize      = SIZE_MAX;
    for (const auto& cId : components) {
        const auto& comp = getNodeFromIdx(cId);
        if (comp->clusterSize < smallestSize) {
            smallestComponent = comp->vertexIdx;
            smallestSize      = comp->clusterSize;
        }
    }
    const auto& smallestC = getNodeFromIdx(smallestComponent);
    presentMap.try_emplace(smallestC->vertexIdx, true);
    const auto& bndryNodes = smallestC->boundaryVertices;

    for (const auto& bndryNode : bndryNodes) {
        const auto& nbrs = pcm->getNbrs(bndryNode);
        for (const auto& nbr : nbrs) {
            fusionEdges.emplace_back(bndryNode, nbr);
        }
    }
}

void UFHeuristic::singleClusterRandomFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                 std::unordered_map<std::size_t, bool>&            presentMap,
                                                 const std::unordered_set<std::size_t>&            components,
                                                 const std::unique_ptr<ParityCheckMatrix>&         pcm) {
    std::random_device rd;
    std::mt19937       gen(rd());
    if (components.size() > std::numeric_limits<int>::max()) {
        throw QeccException("cannot setup distribution, size too large for function");
    }
    std::uniform_int_distribution d(static_cast<std::size_t>(0U), components.size() - 1);
    const std::size_t             chosenIdx = d(gen);
    auto                          it        = components.begin();
    std::advance(it, chosenIdx);
    auto        chosenComponent = *it;
    const auto& chosenNode      = getNodeFromIdx(chosenComponent);

    presentMap.try_emplace(chosenNode->vertexIdx, true);
    const auto& bndryNodes = chosenNode->boundaryVertices;

    for (const auto& bndryNode : bndryNodes) {
        const auto& nbrs = pcm->getNbrs(bndryNode);
        for (const auto& nbr : nbrs) {
            fusionEdges.emplace_back(bndryNode, nbr);
        }
    }
}

/**
 * Computes interior of erasure with BFS algorithm then iterates over interior and removes neighbours of check vertices iteratively
 * @param erasure
 * @param syndrome
 * @return
 */
std::vector<std::size_t> UFHeuristic::erasureDecoder(std::unordered_set<std::size_t>& erasure, std::unordered_set<std::size_t>& syndr, const std::unique_ptr<ParityCheckMatrix>& pcm) {
    std::unordered_set<std::size_t>       syndrome(syndr.begin(), syndr.end());
    std::vector<std::vector<std::size_t>> erasureSet{};
    // compute interior of grown erasure components, that is nodes all of whose neighbours are also in the component
    for (const auto& currCompRootId : erasure) {
        std::vector<std::size_t> compErasure;
        const auto&              currCompRoot = getNodeFromIdx(currCompRootId);
        std::queue<std::size_t>  queue;

        // start traversal at component root
        queue.push(currCompRoot->vertexIdx);

        while (!queue.empty()) {
            const auto& currV = getNodeFromIdx(queue.front());
            queue.pop();
            if ((!currV->marked && (currCompRoot->boundaryVertices.find(currV->vertexIdx)) == currCompRoot->boundaryVertices.end()) || currV->isCheck) { // we need check nodes also if they are not in the "interior" or if there is only a restricted interior
                // add to interior by adding it to the list and marking it
                currV->marked = true;
                compErasure.emplace_back(currV->vertexIdx);
            }

            for (const auto& node : currV->children) {
                if ((!node->marked && (currCompRoot->boundaryVertices.find(node->vertexIdx)) == currCompRoot->boundaryVertices.end()) || node->isCheck) { // we need check nodes also if they are not in the "interior" or if there is only a restricted interior
                    // add to interior by adding it to the list and marking it
                    node->marked = true;
                    compErasure.emplace_back(node->vertexIdx);
                }
                queue.push(node->vertexIdx); // step into depth
            }
        }
        erasureSet.emplace_back(compErasure);
    }

    std::vector<std::size_t> resList;
    // go through nodes in erasure
    // if current node v is a bit node in Int, remove adjacent check nodes c_i and B(c_i,1)
    for (auto& component : erasureSet) {
        auto compNodeIt = component.begin();
        while (compNodeIt != component.end() && !syndrome.empty()) {
            const auto& currN = getNodeFromIdx(*compNodeIt++);
            if (!currN->isCheck && !currN->deleted) {
                resList.emplace_back(currN->vertexIdx); // add bit node to estimate
                // if we add a bit node we have to delete adjacent check nodes and their neighbours
                for (const auto& adjCheck : pcm->getNbrs(currN->vertexIdx)) {
                    const auto& adjCheckNode = getNodeFromIdx(adjCheck);
                    if (adjCheckNode->marked && !adjCheckNode->deleted) {
                        const auto& nNbrs = pcm->getNbrs(adjCheck);
                        auto        nnbr  = nNbrs.begin();
                        // remove bit nodes adjacent to neighbour check
                        while (nnbr != nNbrs.end()) {
                            const auto& nnbrNode = getNodeFromIdx(*nnbr++);
                            if (!nnbrNode->deleted && nnbrNode->marked) { // also removes currN
                                nnbrNode->deleted = true;
                            }
                        }
                        // remove check from syndrome and component
                        syndrome.erase(adjCheckNode->vertexIdx);
                        adjCheckNode->deleted = true;
                    }
                }
                // finally, remove node bitnode
                currN->deleted = true;
            }
        }
    }

    return resList;
}

/**
 * Add those components that are valid to the erasure
 * @param invalidComponents contains components to check validity for
 * @param validComponents contains valid components (including possible new ones at end of function)
 */
void UFHeuristic::extractValidComponents(std::unordered_set<std::size_t>& invalidComponents, std::unordered_set<std::size_t>& validComponents, const std::unique_ptr<ParityCheckMatrix>& pcm) {
    auto it = invalidComponents.begin();
    while (it != invalidComponents.end()) {
        if (isValidComponent(*it, pcm)) {
            validComponents.insert(*it);
            it = invalidComponents.erase(it);
        } else {
            it++;
        }
    }
}

// for each check node verify that there is no neighbour that is in the boundary of the component
// if there is no neighbour in the boundary for each check vertex the check is covered by a node in Int TODO prove this in paper
bool UFHeuristic::isValidComponent(const std::size_t& compId, const std::unique_ptr<ParityCheckMatrix>& pcm) {
    const auto& compNode = getNodeFromIdx(compId);
    gf2Vec      valid(compNode->checkVertices.size());
    std::size_t i = 0U;
    for (const auto& checkVertex : compNode->checkVertices) {
        const auto& nbrs = pcm->getNbrs(checkVertex);
        for (const auto& nbr : nbrs) {
            if (compNode->boundaryVertices.find(nbr) == compNode->boundaryVertices.end()) {
                valid.at(i) = true;
                break;
            }
        }
        i++;
    }
    return std::all_of(valid.begin(), valid.end(), [](bool v) { return v; });
}

// return raw ptr to leave ownership in list
TreeNode* UFHeuristic::getNodeFromIdx(const std::size_t idx) {
    if (auto nodeIt = nodeMap.find(idx); nodeIt != nodeMap.end()) {
        return nodeIt->second.get();
    }

    auto treeNode = std::make_unique<TreeNode>(idx);
    // determine if idx is a check
    if (idx >= getCode()->getN()) {
        treeNode->isCheck = true;
        treeNode->checkVertices.emplace_back(treeNode->vertexIdx);
    }
    const auto ins = nodeMap.try_emplace(idx, std::move(treeNode));
    return ins.first->second.get();
}

void UFHeuristic::reset() {
    for (auto& n : nodeMap) {
        for (auto& c : n.second->children) {
            c = nullptr;
        }
        n.second->children.clear();
        n.second->parent = nullptr;
        n.second.reset();
    }
    nodeMap.clear();
    this->result = {};
    this->growth = GrowthVariant::AllComponents;
    this->getCode()->gethZ()->nbrCache.clear();
    if (this->getCode()->gethZ()) {
        this->getCode()->gethZ()->nbrCache.clear();
    }
}
