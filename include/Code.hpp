//
// Created by luca on 26/04/2022.
//
#include <vector>
#include "TreeNode.hpp"

#ifndef QUNIONFIND_CODE_HPP
#define QUNIONFIND_CODE_HPP

#endif //QUNIONFIND_CODE_HPP
using nodeIdx = size_t;

struct Graph {
    std::vector<std::vector<bool>> adjMatrix;
    std::vector<std::vector<TreeNode>> adjMatrixNodes;

    std::vector<TreeNode> getNeighbours(const TreeNode &node) {
        std::vector<TreeNode> result;
        for (size_t i = 0; i < adjMatrix.at(node.vertexIdx).size(); i++) {
            if (adjMatrix.at(node.vertexIdx).at(i)) {
                result.emplace_back(adjMatrixNodes.at(node.vertexIdx).at(i));
            }
        }
        return result;
    }
    std::vector<std::size_t> getNeighboursIdx(const std::size_t &node) {
        std::vector<std::size_t> result;
        for (auto && i : adjMatrix.at(node)) {
            if (i) {
                result.emplace_back(i);
            }
        }
        return result;
    }
};

class Code {
public:
    Graph tannerGraph;
};