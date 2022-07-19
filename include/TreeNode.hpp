//
// Created by lucas on 26/04/2022.
//

#ifndef QUNIONFIND_TREENODE_HPP
#define QUNIONFIND_TREENODE_HPP

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <vector>

/**
 * Union Find data structure as tree with additional information in root node of tree
 */
class TreeNode {
public:
    std::size_t                            vertexIdx = 0U;
    bool                                   isCheck   = false;
    std::shared_ptr<TreeNode>              parent    = nullptr;
    std::vector<std::shared_ptr<TreeNode>> children{};
    size_t                                 clusterSize = 1U;
    std::set<std::size_t>                  boundaryVertices{};
    std::set<std::size_t>                  checkVertices{};
    // for interior calculation
    std::set<std::size_t> markedNeighbours{};
    bool                  marked  = false;
    bool                  deleted = false;

    TreeNode():
        TreeNode(-1) {}

    explicit TreeNode(const std::size_t vertexIdx):
        vertexIdx(vertexIdx) {
        boundaryVertices.emplace(vertexIdx);
    }

    /*
 * Recursive find using path compression
 */
    static std::shared_ptr<TreeNode> Find(std::shared_ptr<TreeNode>& node) {
        std::cout << "in find" << std::endl;
        while(node->parent != nullptr && node->parent->parent != nullptr){
            auto parent = node->parent;
            node->parent = node->parent->parent;
            node = parent;
        }
        if(node->parent == nullptr){
            return node;
        }else{
            return node->parent;
        }
    }
    /*
     * Merge two trees with given roots
     */
    static void Union(std::shared_ptr<TreeNode>& tree1, std::shared_ptr<TreeNode>& tree2) {
        std::cout << "in union " << std::endl;
        auto root1 = Find(tree1);
        auto root2 = Find(tree2);

        if (*root1 == *root2) {
            return;
        }

        if (root1->clusterSize <= root2->clusterSize) {
            addFirstToSecondTree(root1, root2);
        } else {
            addFirstToSecondTree(root2, root1);
        }
    }

    static void addFirstToSecondTree(std::shared_ptr<TreeNode>& first, const std::shared_ptr<TreeNode>& second) {
        std::cout << "first to 2" << std::endl;
        first->parent = second;
        second->children.emplace_back(first);
        second->clusterSize += first->clusterSize;
        for (const auto& cv: first->checkVertices) {
            second->checkVertices.insert(cv);
        }
        first->checkVertices.clear();
        if (first->isCheck) {
            first->checkVertices.emplace(first->vertexIdx);
        }
        std::cout << "end 1 to 2" << std::endl;
    }

    bool operator==(const TreeNode& other) const {
        return vertexIdx == other.vertexIdx;
    }
    auto operator<=>(const TreeNode& other) const {
        return (vertexIdx <=> other.vertexIdx);
    }

    friend std::ostream& operator<<(std::ostream& os, const TreeNode& v) {
        return os << "idx: " << v.vertexIdx << "parentIx: " << v.parent->vertexIdx << "check: " << v.isCheck << "children: " << v.children;
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
