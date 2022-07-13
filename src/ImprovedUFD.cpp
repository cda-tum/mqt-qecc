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
std::unordered_set<std::shared_ptr<TreeNode>> ImprovedUFD::computeInitTreeComponents(const gf2Vec& syndrome) {
    std::unordered_set<std::shared_ptr<TreeNode>> result{};
    for (std::size_t i = 0; i < syndrome.size(); i++) {
        if (syndrome.at(i)) {
            auto syndrNode     = std::make_shared<TreeNode>(TreeNode(i + code->getN()));
            syndrNode->isCheck = true;
            syndrNode->checkVertices.insert(syndrNode->vertexIdx);
            nodeMap.insert(std::make_pair(syndrNode->vertexIdx, syndrNode));
            result.insert(syndrNode);
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
        std::vector<std::shared_ptr<TreeNode>> erasure;
        while (!invalidComponents.empty()) {
            for (size_t i = 0; i < invalidComponents.size(); i++) {
                // Step 1 growth
                std::vector<std::pair<std::size_t, std::size_t>> fusionEdges;
                std::map<std::size_t, bool>                      presentMap{}; // for step 4

                if (!this->growth || this->growth == GrowthVariant::ALL_COMPONENTS) {
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
                    std::shared_ptr<TreeNode> n1    = getNodeFromIdx(eIt->first);
                    std::shared_ptr<TreeNode> n2    = getNodeFromIdx(eIt->second);
                    auto                      root1 = TreeNode::Find(n1);
                    auto                      root2 = TreeNode::Find(n2);

                    //compares vertexIdx only
                    if (*root1 == *root2) {
                        eIt = fusionEdges.erase(eIt); // passes eIt to erase
                    } else {
                        std::size_t root1Size = root1->clusterSize;
                        std::size_t root2Size = root2->clusterSize;
                        TreeNode::Union(root1, root2);

                        if (root1Size <= root2Size) { // Step 3 fusion of boundary lists
                            for (auto& boundaryVertex: root1->boundaryVertices) {
                                root2->boundaryVertices.insert(boundaryVertex);
                            }
                            root1->boundaryVertices.clear();
                        } else {
                            for (auto& boundaryVertex: root2->boundaryVertices) {
                                root1->boundaryVertices.insert(boundaryVertex);
                            }
                            root2->boundaryVertices.clear();
                        }
                        eIt++;
                    }
                }

                // Replace nodes in list by their roots avoiding duplicates
                auto it = invalidComponents.begin();
                while (it != invalidComponents.end()) {
                    auto elem = *it;
                    auto root = TreeNode::Find(elem);
                    if (elem->vertexIdx != root->vertexIdx && presentMap.contains(root->vertexIdx)) {
                        // root already in component list, no replacement necessary
                        it = invalidComponents.erase(it);
                    } else {
                        // root of component not yet in list, replace node by its root in components;
                        it = invalidComponents.erase(it);
                        invalidComponents.insert(root);
                    }
                }

                // Update Boundary Lists: remove vertices that are not in boundary anymore
                for (auto& component: invalidComponents) {
                    auto iter = component->boundaryVertices.begin();
                    while (iter != component->boundaryVertices.end()) {
                        auto nbrs     = code->Hz.getNbrs(*iter);
                        auto currNode = getNodeFromIdx(*iter);
                        auto currRoot = TreeNode::Find(currNode);
                        auto nbrIt    = nbrs.begin();
                        auto end      = nbrs.end();
                        while (nbrIt != end) {
                            auto nbr     = getNodeFromIdx(*nbrIt);
                            auto nbrRoot = TreeNode::Find(nbr);
                            if (currRoot->vertexIdx != nbrRoot->vertexIdx) {
                                // if we find one neighbour that is not in the same component the currNode is in the boundary
                                iter++;
                                break;
                            }
                            if (std::next(nbrIt) == end) {
                                // if we have checked all neighbours and found none in another component remove from boundary list
                                //toRemoveList.emplace_back(*iter);
                                iter = component->boundaryVertices.erase(iter);
                                break;
                            } else {
                                nbrIt++;
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
    result.estimBoolVector = gf2Vec(code->getN());
    for (auto re: res) {
        result.estimBoolVector.at(re) = true;
        result.estimNodeIdxVector.emplace_back(re);
    }
}

void ImprovedUFD::standardGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                 std::map<std::size_t, bool>& presentMap, const std::unordered_set<std::shared_ptr<TreeNode>>& components) {
    for (auto& component: components) {
        presentMap.insert(std::make_pair(component->vertexIdx, true));
        assert(component->parent == nullptr); // at this point we can assume that component represents root of the component
        auto bndryNodes = component->boundaryVertices;

        for (const auto& bndryNode: bndryNodes) {
            auto nbrs = code->Hz.getNbrs(bndryNode);
            for (auto& nbr: nbrs) {
                fusionEdges.emplace_back(std::pair(bndryNode, nbr));
            }
        }
    }
}

void ImprovedUFD::singleClusterSmallestFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                   std::map<std::size_t, bool>& presentMap, const std::unordered_set<std::shared_ptr<TreeNode>>& components) {
    std::shared_ptr<TreeNode> smallestComponent;
    std::size_t               smallestSize = SIZE_MAX;
    for (const auto& c: components) {
        if (c->clusterSize < smallestSize) {
            smallestComponent = c;
        }
    }
    presentMap.insert(std::make_pair(smallestComponent->vertexIdx, true));
    assert(smallestComponent->parent == nullptr);
    auto bndryNodes = smallestComponent->boundaryVertices;

    for (const auto& bndryNode: bndryNodes) {
        auto nbrs = code->Hz.getNbrs(bndryNode);
        for (auto& nbr: nbrs) {
            fusionEdges.emplace_back(std::pair(bndryNode, nbr));
        }
    }
}

void ImprovedUFD::singleClusterRandomFirstGrowth(std::vector<std::pair<std::size_t, std::size_t>>& fusionEdges,
                                                 std::map<std::size_t, bool>& presentMap, const std::unordered_set<std::shared_ptr<TreeNode>>& components) {
    std::shared_ptr<TreeNode> chosenComponent;
    std::random_device        rd;
    std::mt19937              gen(rd());
    gf2Vec                    result;
    if (components.size() > std::numeric_limits<int>::max()) {
        throw QeccException("cannot setup distribution, size too large for function");
    }
    std::uniform_int_distribution<> d(0U, components.size());
    std::size_t                     chosenIdx = d(gen);
    auto                            it        = components.begin();
    std::advance(it, chosenIdx);
    chosenComponent = *it;

    presentMap.insert(std::make_pair(chosenComponent->vertexIdx, true));
    assert(chosenComponent->parent == nullptr);
    auto bndryNodes = chosenComponent->boundaryVertices;

    for (const auto& bndryNode: bndryNodes) {
        auto nbrs = code->Hz.getNbrs(bndryNode);
        for (auto& nbr: nbrs) {
            fusionEdges.emplace_back(std::pair(bndryNode, nbr));
        }
    }
}

/**
 * Computes interior of erasure with BFS algorithm then iterates over interior and removes neighbours of check vertices iteratively
 * @param erasure
 * @param syndrome
 * @return
 */
std::unordered_set<std::size_t> ImprovedUFD::erasureDecoder(std::vector<std::shared_ptr<TreeNode>>& erasure, std::unordered_set<std::shared_ptr<TreeNode>>& syndrome) {
    std::unordered_set<std::shared_ptr<TreeNode>> interior;
    std::vector<std::unordered_set<std::size_t>>  erasureSet{};
    std::size_t                                   erasureSetIdx = 0;
    // compute interior of grown erasure components, that is nodes all of whose neighbours are also in the component
    for (auto& currCompRoot: erasure) {
        std::unordered_set<std::size_t> compErasure;
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
            for (auto& node: chldrn) {
                if ((!node->marked && !currCompRoot->boundaryVertices.contains(node->vertexIdx)) || node->isCheck) { // we need check nodes also if they are not in the "interior" or if there is only a restriced interior
                    // add to interior by adding it to the list and marking it
                    node->marked = true;
                    compErasure.insert(node->vertexIdx);
                }
                queue.push(node); // step into depth
            }
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
            std::cout << "not empty, next node" << std::endl;
            auto currN = getNodeFromIdx(*compNodeIt++);
            if (!currN->isCheck && !currN->deleted) {
                xi.insert(currN->vertexIdx); // add bit node to estimate
                // if we add a bit node we have to delete adjacent check nodes and their neighbours
                for (auto& adjCheck: code->Hz.getNbrs(currN->vertexIdx)) {
                    auto adjCheckNode = getNodeFromIdx(adjCheck);
                    if (adjCheckNode->marked && !adjCheckNode->deleted) {
                        auto nNbrs = code->Hz.getNbrs(adjCheck);
                        auto nnbr  = nNbrs.begin();
                        // remove bit nodes adjacent to neighbour check
                        while (nnbr != nNbrs.end()) {
                            auto nnbrNode = getNodeFromIdx(*nnbr++);
                            if (!nnbrNode->deleted && nnbrNode->marked) { //also removes currN
                                nnbrNode->deleted = true;
                            }
                        }
                        // remove check from syndrome and component
                        syndrome.erase(adjCheckNode);
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
    for (auto& i: resList) {
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
void ImprovedUFD::extractValidComponents(std::unordered_set<std::shared_ptr<TreeNode>>& invalidComponents, std::vector<std::shared_ptr<TreeNode>>& validComponents) {
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
bool ImprovedUFD::isValidComponent(const std::shared_ptr<TreeNode>& component) {
    gf2Vec      valid(component->checkVertices.size());
    std::size_t i = 0;
    for (const auto& checkVertex: component->checkVertices) {
        auto nbrs = code->Hz.getNbrs(checkVertex);

        for (auto& nbr: nbrs) {
            if (std::find(component->boundaryVertices.begin(), component->boundaryVertices.end(), nbr) == component->boundaryVertices.end()) {
                valid.at(i) = true;
                break;
            }
        }
        i++;
    }
    return std::all_of(valid.begin(), valid.end(), [](bool i) { return i; });
}
std::shared_ptr<TreeNode> ImprovedUFD::getNodeFromIdx(const std::size_t idx) {
    if (nodeMap.contains(idx)) {
        return nodeMap.at(idx);
    } else {
        auto res = std::make_shared<TreeNode>(TreeNode(idx));
        // determine if idx is a check
        if (idx >= code->getN()) {
            res->isCheck = true;
            res->checkVertices.insert(res->vertexIdx);
        }
        nodeMap.insert(std::make_pair(res->vertexIdx, res));
        return res;
    }
}
