//
// Created by luca on 26/04/2022.
//
#include "TreeNode.hpp"

#include <utility>
#include <vector>

#ifndef QUNIONFIND_CODE_HPP
    #define QUNIONFIND_CODE_HPP

struct ParityCheckMatrix {
    explicit ParityCheckMatrix(std::vector<std::vector<bool>> pcm):
        pcm(std::move(pcm)) {}
    std::vector<std::vector<bool>> pcm;
};

struct TannerGraph {
    std::vector<std::vector<bool>>                      adjMatrix;
    std::vector<std::vector<std::shared_ptr<TreeNode>>> adjListNodes;

    std::shared_ptr<TreeNode> getNodeForId(const std::size_t vertexId) {
        return adjListNodes.at(vertexId).at(0);
    }
    std::set<std::shared_ptr<TreeNode>> getNeighbours(const std::shared_ptr<TreeNode>& node) {
        std::set<std::shared_ptr<TreeNode>> result;
        for (size_t i = 1; i < adjListNodes.at(node->vertexIdx).size(); i++) { // at pos 0 is node itself
            result.insert(adjListNodes.at(node->vertexIdx).at(i));
        }
        return result;
    }
    std::set<std::shared_ptr<TreeNode>> getNeighbours(const std::size_t vertexId) {
        std::set<std::shared_ptr<TreeNode>> result;
        for (size_t i = 1; i < adjListNodes.at(vertexId).size(); i++) { // at pos 0 is node itself
            result.insert(adjListNodes.at(vertexId).at(i));
        }
        return result;
    }
    std::vector<std::size_t> getNeighboursIdx(const std::size_t& node) {
        std::vector<std::size_t> result;
        for (size_t i = 1; i < adjListNodes.at(node).size(); i++) { // at pos 0 is node itself
            result.emplace_back(adjListNodes.at(node).at(i)->vertexIdx);
        }
        return result;
    }

    friend std::ostream& operator<<(std::ostream& os, TannerGraph const& c) {
        for (const auto& i: c.adjMatrix) {
            os << "| ";
            for (bool j: i) {
                os << j << " ";
            }
            os << "|";
            os << std::endl;
        }
        return os;
    }
};
class Code {
public:
    TannerGraph       tannerGraph;
    ParityCheckMatrix Hx;

    explicit Code(ParityCheckMatrix hx):
        Hx(std::move(hx)) {
        auto                                                nrChecks = Hx.pcm.size();
        auto                                                nrData   = Hx.pcm.at(0).size();
        std::size_t                                         dim      = nrChecks + nrData;
        std::vector<std::vector<bool>>                      adjMatrBool(dim); // todo this contains bool values only not adjacency list check alg for errors
        std::vector<std::vector<std::shared_ptr<TreeNode>>> adjLstNodes(dim);
        std::map<std::size_t, std::shared_ptr<TreeNode>>    nodeMap;
        for (size_t i = 0; i < dim; i++) {
            std::shared_ptr<TreeNode> n = std::make_shared<TreeNode>(TreeNode(i));
            if (i >= nrData) {
                n->isCheck = true;
                n->checkVertices.insert(n->vertexIdx);
            }
            nodeMap.insert(std::make_pair(i, n));
        }
        for (size_t i = 0; i < dim; i++) {
            std::vector<bool>                      rowBool(dim);
            std::vector<std::shared_ptr<TreeNode>> nbrList;
            nbrList.emplace_back(nodeMap.at(i)); // adjacency list of node n contains n in first position
            if (i < dim - nrChecks) {
                for (size_t j = 0; j < nrChecks; j++) {
                    auto val               = Hx.pcm.at(j).at(i);
                    rowBool.at(nrData + j) = val;
                    if (val) {
                        nbrList.emplace_back(nodeMap.at(nrData + j));
                    }
                }
            } else {
                for (size_t j = 0; j < nrData; j++) {
                    auto val      = Hx.pcm.at(i - nrData).at(j);
                    rowBool.at(j) = val;
                    if (val) {
                        nbrList.emplace_back(nodeMap.at(j)); // insert index of vertex here to generate adjacency list containing the nodes
                    }
                }
            }
            adjMatrBool.at(i) = rowBool;
            adjLstNodes.at(i) = nbrList;
        }
        tannerGraph.adjMatrix    = adjMatrBool;
        tannerGraph.adjListNodes = adjLstNodes;
    }

    std::size_t getN() {
        if (!Hx.pcm.empty()) {
            return Hx.pcm.at(0).size();
        } else {
            return 0;
        }
    }
    static std::vector<std::vector<bool>> rectMatrixMultiply(const std::vector<std::vector<bool>>& m1, const std::vector<std::vector<bool>>& m2) {
        std::vector<std::vector<bool>> result(m1.size());
        for (std::size_t i = 0; i < m1.size(); i++) {
            result.at(i) = std::vector<bool>(m2.at(0).size());
            for (std::size_t j = 0; j < m2.at(0).size(); j++) {
                result[i][j] = false;

                for (std::size_t k = 0; k < m2.size(); k++) {
                    result[i][j] = result[i][j] + (m1[i][k] * m2[k][j]);
                }
            }
        }
        return result;
    }
    std::vector<bool> getSyndrome(const std::vector<bool>& err) {
        std::vector<std::vector<bool>> errMat(err.size());
        for (size_t i = 0; i < err.size(); i++) {
            errMat.at(i) = std::vector<bool>{err.at(i)}; //transpose
        }
        auto res = rectMatrixMultiply(Hx.pcm, errMat);
        if (!res.empty()) {
            std::vector<bool> rres(Hx.pcm.size());
            for (size_t i = 0; i < rres.size(); i++) {
                rres.at(i) = res.at(i).at(0); // transpose back
            }
            return rres;
        } else {
            return std::vector<bool>{};
        }
    }
    friend std::ostream& operator<<(std::ostream& os, Code const& c) {
        return os << c.tannerGraph;
    }
};

inline std::ostream& operator<<(std::ostream& os, const std::vector<bool>& v) {
    if (v.empty()) {
        os << "[]";
        return os;
    }
    os << "[";
    for (const auto& i: v) {
        os << i;
        os << ", ";
    }
    os << "\b\b";
    os << "]";
    return os;
}
#endif //QUNIONFIND_CODE_HPP
