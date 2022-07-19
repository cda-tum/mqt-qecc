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
    std::size_t                          vertexIdx = 0U;
    bool                                 isCheck   = false;
    std::weak_ptr<TreeNode>              parent;
    std::vector<std::weak_ptr<TreeNode>> children{};
    size_t                               clusterSize = 1U;
    std::set<std::size_t>                boundaryVertices{};
    std::set<std::size_t>                checkVertices{};
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
    static std::weak_ptr<TreeNode> Find(std::weak_ptr<TreeNode>& node) {
        //std::cout << "in find" << std::endl;
        std::shared_ptr<TreeNode> n = node.lock();
        std::shared_ptr<TreeNode> parent = node.lock()->parent.lock();
        while(parent != nullptr && parent->parent.lock() != nullptr){
            auto p = parent;
            parent = parent->parent.lock();
            node         = p;
        }
        if(parent == nullptr){
            return node;
        } else {
            return parent;
        }
    }
    /*
     * Merge two trees with given roots
     */
    static void Union(std::weak_ptr<TreeNode>& tree1, std::weak_ptr<TreeNode>& tree2) {
        auto root1 = Find(tree1);
        auto root2 = Find(tree2);
        auto r1 = root1.lock();
        auto r2 = root2.lock();

        if (*r1 == *r2) {
            return;
        }

        if (r1->clusterSize <= r2->clusterSize) {
            addFirstToSecondTree(root1, root2);
        } else {
            addFirstToSecondTree(root2, root1);
        }
    }

    static void addFirstToSecondTree(std::weak_ptr<TreeNode>& first, const std::weak_ptr<TreeNode>& second) {
        auto f = first.lock();
        auto s = second.lock();
        f->parent = s;
        s->children.emplace_back(first);
        s->clusterSize += f->clusterSize;
        for (const auto& cv: f->checkVertices) {
            s->checkVertices.insert(cv);
        }
        f->checkVertices.clear();
        if (f->isCheck) {
            f->checkVertices.emplace(f->vertexIdx);
        }
    }

    bool operator==(const TreeNode& other) const {
        return vertexIdx == other.vertexIdx;
    }
    auto operator<=>(const TreeNode& other) const {
        return (vertexIdx <=> other.vertexIdx);
    }

    friend std::ostream& operator<<(std::ostream& os, const TreeNode& v) {
        return os << "idx: " << v.vertexIdx << "parentIx: " << v.parent.lock()->vertexIdx << "check: " << v.isCheck << "children: " << v.children.size();
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
