//
// Created by lucas on 21/04/2022.
//

#include "ImprovedUFD.hpp"

#include "Decoder.hpp"
#include "TreeNode.hpp"

#include <cassert>
#include <chrono>
#include <queue>
#include <random>
/**
     * returns list of tree node (in UF data structure) representations for syndrome
     * @param code
     * @param syndrome
     * @return
    */
std::unordered_set<std::size_t> ImprovedUFD::computeInitTreeComponents(const gf2Vec& syndrome) {
    std::unordered_set<std::size_t> result;
    for (std::size_t i = 0; i < syndrome.size(); i++) {
        if (syndrome.at(i)) {
            auto syndrNode     = std::make_shared<TreeNode>(i + getCode()->getN());
            syndrNode->isCheck = true;
            syndrNode->checkVertices.emplace(syndrNode->vertexIdx);
            nodeMap.try_emplace(syndrNode->vertexIdx, syndrNode);
            result.emplace(syndrNode->vertexIdx);
        }
    }
    return result;
}

void ImprovedUFD::decode(gf2Vec& syndrome) {
    auto                            decodingTimeBegin = std::chrono::high_resolution_clock::now();
    std::unordered_set<std::size_t> res;

    if (!syndrome.empty() && !std::all_of(syndrome.begin(), syndrome.end(), [](bool val) { return !val; })) {
        auto                                   syndrComponents   = computeInitTreeComponents(syndrome);
        auto                                   invalidComponents = syndrComponents;
        std::vector<std::size_t> erasure;
        while (!invalidComponents.empty()) {
            for (size_t i = 0; i < invalidComponents.size(); i++) {
                // Step 1 growth
                std::vector<std::pair<std::size_t, std::size_t>> fusionEdges;
                std::map<std::size_t, bool>                      presentMap{}; // for step 4

                if (this->growth == GrowthVariant::ALL_COMPONENTS) {
                    // to grow all components (including valid ones)
                    invalidComponents.insert(erasure.begin(), erasure.end());
                    standardGrowth(fusionEdges, presentMap, invalidComponents);
                } else if (this->growth == GrowthVariant::INVALID_COMPONENTS) {
                    standardGrowth(fusionEdges, presentMap, invalidComponents);
                } else if (this->growth == GrowthVariant::SINGLE_SMALLEST) {
                    singleClusterSmallestFirstGrowth(fusionEdges, presentMap, invalidComponents);
                } else if (this->growth == GrowthVariant::SINGLE_RANDOM) {
                    singleClusterRandomFirstGrowth(fusionEdges, presentMap, invalidComponents);
                } else {
                    throw std::invalid_argument("Unsupported growth variant");
                }
                // Fuse clusters that grew together
                auto eIt = fusionEdges.begin();
                while (eIt != fusionEdges.end()) {
                    const auto& n1    = getNodeFromIdx(eIt->first).lock();
                    const auto& n2    = getNodeFromIdx(eIt->second).lock();
                    auto        root1 = TreeNode::Find(n1);
                    auto        root2 = TreeNode::Find(n2);

                    //compares vertexIdx only
                    if (*root1 == *root2) {
                        eIt = fusionEdges.erase(eIt); // passes eIt to erase
                    } else {
                        const std::size_t root1Size = root1->clusterSize;
                        const std::size_t root2Size = root2->clusterSize;
                        TreeNode::Union(root1, root2);

                        if (root1Size <= root2Size) { // Step 3 fusion of boundary lists
                            for (const auto& boundaryVertex: root1->boundaryVertices) {
                                root2->boundaryVertices.insert(boundaryVertex);
                            }
                            root1->boundaryVertices.clear();
                        } else {
                            for (const auto& boundaryVertex: root2->boundaryVertices) {
                                root1->boundaryVertices.insert(boundaryVertex);
                            }
                            root2->boundaryVertices.clear();
                        }
                        ++eIt;
                    }
                }

                // Replace nodes in list by their roots avoiding duplicates
                auto idxIt = invalidComponents.begin();
                while (idxIt != invalidComponents.end()) {
                    auto elem = getNodeFromIdx(*idxIt).lock();
                    auto root = TreeNode::Find(elem);
                    if (elem->vertexIdx != root->vertexIdx && presentMap.contains(root->vertexIdx)) {
                        // root already in component list, no replacement necessary
                        idxIt = invalidComponents.erase(idxIt);
                    } else {
                        // root of component not yet in list, replace node by its root in components;
                        idxIt = invalidComponents.erase(idxIt);
                        invalidComponents.emplace(root->vertexIdx);
                    }
                }

                // Update Boundary Lists: remove vertices that are not in boundary anymore
                for (auto& compId: invalidComponents) {
                    auto compNode = getNodeFromIdx(compId).lock();
                    auto iter = compNode->boundaryVertices.begin();
                    while (iter != compNode->boundaryVertices.end()) {
                        const auto nbrs     = getCode()->Hz.getNbrs(*iter);
                        const auto currNode = getNodeFromIdx(*iter).lock();
                        const auto currRoot = TreeNode::Find(currNode);
                        for (const auto& nbr: nbrs) {
                            const auto node = getNodeFromIdx(nbr).lock();
                            if (const auto nbrRoot = TreeNode::Find(node); currRoot->vertexIdx != nbrRoot->vertexIdx) {
                                // if we find one neighbour that is not in the same component the currNode is in the boundary
                                iter++;
                                break;
                            }
                            if (nbr == nbrs.back()) {
                                // if we have checked all neighbours and found none in another component remove from boundary list
                                //toRemoveList.emplace_back(*iter);
                                iter = compNode->boundaryVertices.erase(iter);
                                break;
                            }
                        }
                    }
                }
                extractValidComponents(invalidComponents, erasure);
            }
        }
        res = erasureDecoder(erasure, syndrComponents);
    }
    auto decodingTimeEnd   = std::chrono::high_resolution_clock::now();
    this->result           = DecodingResult();
    result.decodingTime    = std::chrono::duration_cast<std::chrono::milliseconds>(decodingTimeEnd - decodingTimeBegin).count();
    result.estimBoolVector = gf2Vec(getCode()->getN());
    for (auto re: res) {
        result.estimBoolVector.at(re) = true;
        result.estimNodeIdxVector.emplace_back(re);
    }
}

void ImprovedUFD::standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                 std::map<std::size_t, bool>& presentMap,
                                 const std::unordered_set<std::size_t>& components) {
    for (const auto& compId: components) {
        auto compNode = getNodeFromIdx(compId).lock();
        presentMap.try_emplace(compNode->vertexIdx, true);
        assert(compNode->parent == nullptr); // at this point we can assume that component represents root of the component
        const auto& bndryNodes = compNode->boundaryVertices;

        for (const auto& bndryNode: bndryNodes) {
            const auto nbrs = getCode()->Hz.getNbrs(bndryNode);
            for (const auto& nbr: nbrs) {
                fusionEdges.emplace_back(bndryNode, nbr);
            }
        }
    }
}

void ImprovedUFD::singleClusterSmallestFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                   std::map<std::size_t, bool>& presentMap, const std::unordered_set<std::size_t>& components) {
    std::shared_ptr<TreeNode> smallestComponent;
    std::size_t               smallestSize = SIZE_MAX;
    for (const auto& cId: components) {
        auto comp = getNodeFromIdx(cId).lock();
        if (comp->clusterSize < smallestSize) {
            smallestComponent = comp;
        }
    }
    presentMap.try_emplace(smallestComponent->vertexIdx, true);
    assert(smallestComponent->parent == nullptr);
    const auto& bndryNodes = smallestComponent->boundaryVertices;

    for (const auto& bndryNode: bndryNodes) {
        const auto nbrs = getCode()->Hz.getNbrs(bndryNode);
        for (const auto& nbr: nbrs) {
            fusionEdges.emplace_back(bndryNode, nbr);
        }
    }
    smallestComponent = {};
}

void ImprovedUFD::singleClusterRandomFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                 std::map<std::size_t, bool>& presentMap,
                                                 const std::unordered_set<std::size_t>& components) {
    std::size_t chosenComponent;
    std::random_device        rd;
    std::mt19937              gen(rd());
    gf2Vec                    result;
    if (components.size() > std::numeric_limits<int>::max()) {
        throw QeccException("cannot setup distribution, size too large for function");
    }
    std::uniform_int_distribution d(static_cast<std::size_t>(0U), components.size());
    std::size_t                   chosenIdx = d(gen);
    auto                          it        = components.begin();
    std::advance(it, chosenIdx);
    chosenComponent = *it;
    auto chosenNode = getNodeFromIdx(chosenComponent).lock();

    presentMap.try_emplace(chosenNode->vertexIdx, true);
    assert(chosenNode->parent == nullptr);
    const auto& bndryNodes = chosenNode->boundaryVertices;

    for (const auto& bndryNode: bndryNodes) {
        const auto nbrs = getCode()->Hz.getNbrs(bndryNode);
        for (const auto& nbr: nbrs) {
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
std::unordered_set<std::size_t> ImprovedUFD::erasureDecoder(std::vector<std::size_t>& erasure, std::unordered_set<std::size_t>& syndrome) {
    std::vector<std::unordered_set<std::size_t>>  erasureSet{};
    std::size_t                                   erasureSetIdx = 0;
    // compute interior of grown erasure components, that is nodes all of whose neighbours are also in the component
    for (auto& currCompRootId: erasure) {
        std::unordered_set<std::size_t> compErasure;
        auto currCompRoot = getNodeFromIdx(currCompRootId).lock();
        const auto                      boundaryVertices = currCompRoot->boundaryVertices; // boundary vertices are stored in root of each component
        if (currCompRoot->parent != nullptr) {
            currCompRoot = TreeNode::Find(currCompRoot);
        }
        assert(currCompRoot->parent == nullptr); // should be due to steps above otherwise we can just call Find here
        std::queue<std::shared_ptr<TreeNode>> queue;

        // start traversal at component root
        queue.push(currCompRoot);

        while (!queue.empty()) {
            auto currV = queue.front();
            queue.pop();
            if ((!currV->marked && !currCompRoot->boundaryVertices.contains(currV->vertexIdx)) || currV->isCheck) { // we need check nodes also if they are not in the "interior" or if there is only a restriced interior
                // add to interior by adding it to the list and marking it
                currV->marked = true;
                compErasure.insert(currV->vertexIdx);
            }
            std::vector<std::shared_ptr<TreeNode>> chldrn;
            for (auto& i: currV->children) {
                chldrn.emplace_back(i);
            }
            for (const auto& node: chldrn) {
                if ((!node->marked && !currCompRoot->boundaryVertices.contains(node->vertexIdx)) || node->isCheck) { // we need check nodes also if they are not in the "interior" or if there is only a restriced interior
                    // add to interior by adding it to the list and marking it
                    node->marked = true;
                    compErasure.insert(node->vertexIdx);
                }
                queue.push(node); // step into depth
            }
            chldrn.clear();
        }
        erasureSet.emplace_back(compErasure);
        erasureSetIdx++;
    }


    std::vector<std::unordered_set<std::size_t>> resList;
    // go through nodes in erasure
    // if current node v is a bit node in Int, remove adjacent check nodes c_i and B(c_i,1)
    for (auto& component: erasureSet) {
        std::unordered_set<std::size_t> xi;
        auto                            compNodeIt = component.begin();
        while (compNodeIt != component.end() && !syndrome.empty()) {
            auto currN = getNodeFromIdx(*compNodeIt++).lock();
            if (!currN->isCheck && !currN->deleted) {
                xi.insert(currN->vertexIdx); // add bit node to estimate
                // if we add a bit node we have to delete adjacent check nodes and their neighbours
                for (const auto& adjCheck: getCode()->Hz.getNbrs(currN->vertexIdx)) {
                    auto adjCheckNode = getNodeFromIdx(adjCheck).lock();
                    if (adjCheckNode->marked && !adjCheckNode->deleted) {
                        auto nNbrs = getCode()->Hz.getNbrs(adjCheck);
                        auto nnbr  = nNbrs.begin();
                        // remove bit nodes adjacent to neighbour check
                        while (nnbr != nNbrs.end()) {
                            auto nnbrNode = getNodeFromIdx(*nnbr++).lock();
                            if (!nnbrNode->deleted && nnbrNode->marked) { //also removes currN
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
        resList.emplace_back(xi);
    }

    std::unordered_set<std::size_t> res;
    for (const auto& i: resList) {
        for (auto& n: i) {
            res.insert(n);
        }
    }
    return res;
}

/**
 * Add those components that are valid to the erasure
 * @param invalidComponents containts components to check validity for
 * @param validComponents contains valid components (including possible new ones at end of function)
 */
void ImprovedUFD::extractValidComponents(std::unordered_set<std::size_t>& invalidComponents, std::vector<std::size_t>& validComponents) {
    auto it = invalidComponents.begin();
    while (it != invalidComponents.end()) {
        if (isValidComponent(*it)) {
            validComponents.emplace_back(*it);
            it = invalidComponents.erase(it);
        } else {
            it++;
        }
    }
}

// for each check node verify that there is no neighbour that is in the boundary of the component
// if there is no neighbour in the boundary for each check vertex the check is covered by a node in Int TODO prove this in paper
bool ImprovedUFD::isValidComponent(const std::size_t& compId) {
    auto compNode = getNodeFromIdx(compId).lock();
    gf2Vec      valid(compNode->checkVertices.size());
    std::size_t i = 0;
    for (const auto& checkVertex: compNode->checkVertices) {
        for (const auto nbrs = getCode()->Hz.getNbrs(checkVertex); const auto& nbr: nbrs) {
            if (!compNode->boundaryVertices.contains(nbr)) {
                valid.at(i) = true;
                break;
            }
        }
        i++;
    }
    return std::all_of(valid.begin(), valid.end(), [](bool i) { return i; });
}
// return weak_ptr to leave ownership in list
std::weak_ptr<TreeNode> ImprovedUFD::getNodeFromIdx(const std::size_t idx) {
    if (nodeMap.contains(idx)) {
        return nodeMap.at(idx);
    } else {
        auto res = std::make_shared<TreeNode>(idx);
        // determine if idx is a check
        if (idx >= getCode()->getN()) {
            res->isCheck = true;
            res->checkVertices.emplace(res->vertexIdx);
        }
        nodeMap.try_emplace(res->vertexIdx, res);
        return res;
    }
}

void ImprovedUFD::reset() {
    this->nodeMap.clear();
    this->result = {};
}
