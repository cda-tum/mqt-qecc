//
// Created by luca on 26/04/2022.
//

#ifndef QUNIONFIND_TREENODE_HPP
#define QUNIONFIND_TREENODE_HPP

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <vector>

class TreeNode {
public:
    std::size_t               vertexIdx = 0;
    bool                      isCheck   = false;
    std::shared_ptr<TreeNode> parent    = nullptr;
    std::vector<std::shared_ptr<TreeNode>>     children{};
    size_t                    clusterSize = 1;
    std::set<std::shared_ptr<TreeNode>>        boundaryVertices{};
    std::set<std::shared_ptr<TreeNode>>        checkVertices{};

    TreeNode():
        TreeNode(-1) {}
    TreeNode(TreeNode&&) = default; // forces a move constructor anyway

    TreeNode(const TreeNode& other) {
        vertexIdx = other.vertexIdx;
        isCheck   = other.isCheck;
        if (other.parent) {
            parent = other.parent;
        } else {
            parent = nullptr;
        }
        for (const auto& i: other.children) {
            children.emplace_back(i);
        }
        clusterSize      = other.clusterSize;
        boundaryVertices = other.boundaryVertices;
        checkVertices    = other.checkVertices;
    }

    explicit TreeNode(const std::size_t vertexIdx) { this->vertexIdx = vertexIdx; }
    /*
     * Recursive Find using path compression
     */
    static std::shared_ptr<TreeNode> Find(std::shared_ptr<TreeNode>& node) {
        std::cout << "in find" << std::endl;
        if (node->parent == nullptr) {
            return node;
        } else {
            node->parent = std::make_shared<TreeNode>(*Find(node->parent));
            return node->parent;
        }
    }

    /*
     * Merge two trees with given roots
     */
    static void Union(std::shared_ptr<TreeNode>& tree1, std::shared_ptr<TreeNode>& tree2) {
        std::cout << "in union" << std::endl;
        auto root1 = Find(tree1);
        auto root2 = Find(tree2);

        if (*root1 == *root2) {
            return;
        }

        if (root1->clusterSize <= root2->clusterSize) {
            root1->parent = root2;
            root2->children.emplace_back(root1);
            root2->clusterSize += root2->clusterSize;
        } else if (root1->clusterSize > root2->clusterSize) {
            root2->parent =root1;
            root1->children.emplace_back(root2);
            root2->clusterSize += root2->clusterSize;
        }
    }
    TreeNode& operator=(const TreeNode& other) {
        if (this == &other) {
            return *this;
        }
        vertexIdx = other.vertexIdx;
        isCheck   = other.isCheck;
        if (other.parent) {
            parent = other.parent;
        } else {
            parent = nullptr;
        }
        for (const auto& i: other.children) {
            children.emplace_back(i);
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
    friend std::ostream& operator<<(std::ostream& os, const std::set<TreeNode>& v) {
        if (v.empty()) {
            os << "[]";
            return os;
        }
        os << "[";
        for (const auto& i: v) {
            os << i.vertexIdx;
            os << ", ";
        }
        os << "\b\b";
        os << "]";
        return os;
    }
    friend std::ostream& operator<<(std::ostream& os, const std::vector<TreeNode>& v) {
        if (v.empty()) {
            os << "[]";
            return os;
        }
        os << "[";
        for (const auto& i: v) {
            os << i.vertexIdx;
            os << ", ";
        }
        os << "\b\b";
        os << "]";
        return os;
    }
    friend std::ostream& operator<<(std::ostream& os, const std::vector<std::shared_ptr<TreeNode>>& v) {
        if (v.empty()) {
            os << "[]";
            return os;
        }
        os << "[";
        for (const auto& i: v) {
            os << i->vertexIdx;
            os << ", ";
        }
        os << "\b\b";
        os << "]";
        return os;
    }
    friend std::ostream& operator<<(std::ostream& os, const std::set<std::shared_ptr<TreeNode>>& v) {
        if (v.empty()) {
            os << "[]";
            return os;
        }
        os << "[";
        for (const auto& i: v) {
            os << i->vertexIdx;
            os << ", ";
        }
        os << "\b\b";
        os << "]";
        return os;
    }
};

#endif //QUNIONFIND_TREENODE_HPP
