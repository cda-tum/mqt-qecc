#pragma once

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <unordered_set>
#include <vector>

/**
 * Union Find data structure as tree with additional information in root node of tree
 */
class TreeNode {
public:
    std::size_t                     vertexIdx = 0U;
    bool                            isCheck   = false;
    TreeNode*                       parent    = nullptr;
    std::vector<TreeNode*>          children;
    std::size_t                     clusterSize = 1U;
    std::unordered_set<std::size_t> boundaryVertices;
    std::vector<std::size_t>        checkVertices;
    bool                            marked  = false;
    bool                            deleted = false;

    TreeNode() : TreeNode(std::numeric_limits<std::size_t>::max()) {}

    explicit TreeNode(const std::size_t& idx) : vertexIdx(idx) {
        boundaryVertices.emplace(idx);
    }

    /*
     * find using path compression
     */
    static TreeNode* Find(TreeNode* node) { // NOLINT(readability-identifier-naming)
        auto* parent = node->parent;
        while (parent != nullptr && parent->parent != nullptr) {
            node   = parent;
            parent = parent->parent;
        }
        if (parent == nullptr) {
            return node;
        }
        return parent;
    }
    /*
     * Merge two trees with given roots
     */
    static void Union(TreeNode* tree1, TreeNode* tree2) { // NOLINT(readability-identifier-naming)
        auto* root1 = Find(tree1);
        auto* root2 = Find(tree2);

        if (root1->vertexIdx == root2->vertexIdx) {
            return;
        }
        if (root1->clusterSize <= root2->clusterSize) {
            addFirstToSecondTree(root1, root2);
        } else {
            addFirstToSecondTree(root2, root1);
        }
    }

    static void addFirstToSecondTree(TreeNode* first, TreeNode* second) {
        first->parent = second;
        second->children.emplace_back(first);
        second->clusterSize += first->clusterSize;
        std::move(first->checkVertices.begin(), first->checkVertices.end(), std::back_inserter(second->checkVertices));

        first->checkVertices.clear();
        if (first->isCheck) {
            first->checkVertices.emplace_back(first->vertexIdx);
        }
    }

    bool operator==(const TreeNode& other) const {
        return vertexIdx == other.vertexIdx;
    }
    auto operator<=(const TreeNode& other) const {
        return (vertexIdx <= other.vertexIdx);
    }

    auto operator>=(const TreeNode& other) const {
        return (vertexIdx >= other.vertexIdx);
    }

    auto operator<(const TreeNode& other) const {
        return (vertexIdx < other.vertexIdx);
    }

    auto operator>(const TreeNode& other) const {
        return (vertexIdx > other.vertexIdx);
    }

    friend std::ostream& operator<<(std::ostream& os, const TreeNode& v) {
        return os << "idx: " << v.vertexIdx << "parentIx: " << v.parent->vertexIdx << "check: " << v.isCheck << "children: " << v.children.size();
    }

    friend std::ostream& operator<<(std::ostream& os, const std::unordered_set<TreeNode>& v) {
        if (v.empty()) {
            os << "[]";
            return os;
        }
        os << "[";
        for (const auto& i : v) {
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
        for (const auto& i : v) {
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
        for (const auto& i : v) {
            os << i->vertexIdx;
            os << ", ";
        }
        os << "\b\b";
        os << "]";
        return os;
    }
    friend std::ostream& operator<<(std::ostream& os, const std::unordered_set<std::shared_ptr<TreeNode>>& v) {
        if (v.empty()) {
            os << "[]";
            return os;
        }
        os << "[";
        for (const auto& i : v) {
            os << i->vertexIdx;
            os << ", ";
        }
        os << "\b\b";
        os << "]";
        return os;
    }
};
