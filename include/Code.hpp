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
    std::size_t K =0;

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

    std::size_t getK(){
        return K;
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

    bool checkStabilizer(const std::vector<bool>& est) {
        return checkVectorInRowspace(Hx.pcm, est);
    }

    static void swapRows(std::vector<std::vector<bool>>& matrix, const std::size_t row1, const std::size_t row2) {
        for (std::size_t col = 0; col <= matrix.at(0).size(); col++) {
            std::swap(matrix.at(row1).at(col), matrix.at(row2).at(col));
        }
    }

    static void printGF2matrix(const std::vector<std::vector<bool>>& matrix) {
        for (const auto & i : matrix) {
            for (size_t j = 0; j < i.size(); j++) {
                std::cout << i.at(j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    static bool checkVectorInRowspace(std::vector<std::vector<bool>> M, std::vector<bool> vec) { //https://stackoverflow.com/questions/11483925/how-to-implementing-gaussian-elimination-for-binary-equations
        std::size_t                    nrCols = M.size();
        std::size_t                    nrRows = M.at(0).size();
        std::size_t                    row    = 0;
        std::vector<std::vector<bool>> tranp(nrRows);
        for (size_t i = 0; i < tranp.size(); i++) {
            tranp.at(i) = std::vector<bool>(nrCols);
        }
        for (size_t i = 0; i < M.size(); i++) {
            for (size_t j = 0; j < M.at(i).size(); j++) {
                tranp[j][i] = M[i][j];
            }
        }
        M = tranp;

        std::cout << "vector: " << std::endl;
        for (size_t i = 0; i < vec.size(); i++) {
            std::cout << vec.at(i) << " ";
        }
        printGF2matrix(M);

        for (std::size_t col = 0; col < nrCols && row < nrRows; col++, row++) {
            if (M[row][col] == 0) {
                for (std::size_t i = 0; i < nrRows; i++) {
                    if (M[i][col] != 0) {
                        swapRows(M, i, row);
                        std::swap(vec.at(i), vec.at(row));
                    }
                }
            }
            if (M[row][col] == 0) {
                return false;
            }
            for (std::size_t j = 0; j < nrRows; ++j) {
                if (j != col) {
                    if (M[j][col]) {
                        std::size_t k;
                        for (k = col; k < nrCols; ++k) {
                            M[j][k] = M[j][k] ^ M[col][k];
                        }
                        vec[k] = vec[k] ^ vec[col];
                    }
                }
            }
        }
        for (size_t i = 0; i < vec.size(); i++) {
            if (vec[i]) {
                for (size_t j = 0; j < M.at(i).size(); j++) {
                    if (M[i][j]) {
                        return false;
                    }
                }
            }
        }
        return true;
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
