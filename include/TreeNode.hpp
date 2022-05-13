//
// Created by luca on 26/04/2022.
//

#ifndef QUNIONFIND_TREENODE_HPP
#define QUNIONFIND_TREENODE_HPP

#include <map>
#include <memory>
#include <set>
#include <vector>

class TreeNode {
public:
    TreeNode():
        TreeNode(-1) {}
    TreeNode(TreeNode&&) = default; // forces a move constructor anyway

    TreeNode(const TreeNode& other) {
        vertexIdx = other.vertexIdx;
        isCheck   = other.isCheck;
        parent    = std::make_unique<TreeNode>(*other.parent);
        for (const auto& i: other.children) {
            children.emplace_back(new TreeNode(*i));
        }
        clusterSize      = other.clusterSize;
        boundaryVertices = other.boundaryVertices;
        checkVertices    = other.checkVertices;
    }

    TreeNode(const std::size_t vertexIdx) { this->vertexIdx = vertexIdx; }

    TreeNode& operator=(const TreeNode& other) {
        if (this == &other) {
            return *this;
        }
        vertexIdx = other.vertexIdx;
        isCheck   = other.isCheck;
        parent    = std::make_unique<TreeNode>(*other.parent);
        for (const auto& i: other.children) {
            children.emplace_back(new TreeNode(*i));
        }
        clusterSize      = other.clusterSize;
        boundaryVertices = other.boundaryVertices;
        checkVertices    = other.checkVertices;

        return *this;
    }

    bool operator==(const TreeNode& other) const {
        return this->vertexIdx == other.vertexIdx;
    }
    bool operator<(const TreeNode& other) const {
        return this->vertexIdx < other.vertexIdx;
    }

public:
    std::size_t                            vertexIdx = 0;
    bool                                   isCheck   = false;
    std::unique_ptr<TreeNode>              parent{};
    std::vector<std::unique_ptr<TreeNode>> children{};

    size_t             clusterSize = 0;
    std::set<TreeNode> boundaryVertices{};
    std::set<TreeNode> checkVertices{};
};

/*
 * Recursive Find using path compression
 */
TreeNode Find(TreeNode& node) {
    if (!node.parent) {
        return node;
    } else {
        node.parent = std::make_unique<TreeNode>(Find(*node.parent));
        return *node.parent;
    }
}

/*
 * Merge two trees with given roots
 */
void Union(TreeNode& tree1, TreeNode& tree2) {
    auto root1 = Find(tree1);
    auto root2 = Find(tree2);

    if (root1 == root2) {
        return;
    }

    if (root1.clusterSize <= root2.clusterSize) {
        root1.parent = std::make_unique<TreeNode>(root2);
        root2.children.emplace_back(new TreeNode(root1));
    } else if (root1.clusterSize > root2.clusterSize) {
        root2.parent = std::make_unique<TreeNode>(root1);
        root1.children.emplace_back(new TreeNode(root2));
    }
}

#endif //QUNIONFIND_TREENODE_HPP
