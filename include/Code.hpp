//
// Created by luca on 26/04/2022.
//

#ifndef QUNIONFIND_CODE_HPP
#define QUNIONFIND_CODE_HPP
#include "TreeNode.hpp"
#include "Utils.hpp"

#include <utility>
#include <vector>

struct ParityCheckMatrix {
    explicit ParityCheckMatrix(gf2Mat pcm):
        pcm(std::move(pcm)) {}
    const gf2Mat pcm;
};

struct TannerGraph {
    gf2Mat                                              adjMatrix;
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
/**
 * Considers X errors with Z checks only.
 * Z errors are corrected completely analogously in a symmetric way for Z and X.
 */
class Code {
public:
    TannerGraph       tannerGraph;
    ParityCheckMatrix Hz;
    std::size_t       K = 0U;
    std::size_t       N = 0U;

    /*
     * Takes matrix Hz over GF(2) and constructs respective code for X errors with Z checks represented by Hz (Convention: Rows in first dim, columns in second)
     * Additionally, the adjacency matrix representation using the Union Find datastructure is construced for Hz
     */
    explicit Code(ParityCheckMatrix hz):
        Hz(std::move(hz)) {
        auto                                                nrChecks = Hz.pcm.size();
        auto                                                nrData   = Hz.pcm.at(0).size();
        std::size_t                                         dim      = nrChecks + nrData;
        gf2Mat                                              adjMatrBool(dim);
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
            gf2Vec                                 rowBool(dim);
            std::vector<std::shared_ptr<TreeNode>> nbrList;
            nbrList.emplace_back(nodeMap.at(i)); // adjacency list of node n contains n in first position by our convention
            if (i < dim - nrChecks) {
                for (size_t j = 0; j < nrChecks; j++) {
                    auto val               = Hz.pcm.at(j).at(i);
                    rowBool.at(nrData + j) = val;
                    if (val) {
                        nbrList.emplace_back(nodeMap.at(nrData + j));
                    }
                }
            } else {
                for (size_t j = 0; j < nrData; j++) {
                    auto val      = Hz.pcm.at(i - nrData).at(j);
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
        N                        = Hz.pcm.at(0).size();
    }

    [[nodiscard]] std::size_t getN() const {
        return N;
    }

    [[nodiscard]] std::size_t getK() const {
        return K;
    }

    [[nodiscard]] gf2Vec getSyndrome(const gf2Vec& err) const {
        gf2Mat errMat(err.size());
        for (size_t i = 0; i < err.size(); i++) {
            errMat.at(i) = gf2Vec{err.at(i)}; //transpose
        }
        auto res = Utils::rectMatrixMultiply(Hz.pcm, errMat);
        if (!res.empty()) {
            gf2Vec rres(Hz.pcm.size());
            for (size_t i = 0; i < rres.size(); i++) {
                rres.at(i) = res.at(i).at(0); // transpose back
            }
            return rres;
        } else {
            return gf2Vec{};
        }
    }

    [[nodiscard]] bool isVectorStabilizer(const gf2Vec& est) const {
        return Utils::isVectorInRowspace(Hz.pcm, est);
    }

    friend std::ostream& operator<<(std::ostream& os, Code const& c) {
        return os << c.tannerGraph;
    }
};
#endif //QUNIONFIND_CODE_HPP
