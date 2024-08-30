#include "TreeNode.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <string>

class TreeNodeTest : public testing::TestWithParam<std::string> {};

TEST(TreeNodeTest, TestFindDisjoint) {
    auto n1 = std::make_unique<TreeNode>(0);
    auto n2 = std::make_unique<TreeNode>(1);

    auto* n1w = n1.get();
    auto* n2w = n2.get();

    EXPECT_FALSE(TreeNode::Find(n1w) == TreeNode::Find(n2w));
}

TEST(TreeNodeTest, TestFindSameTreeChild) {
    auto n1    = std::make_unique<TreeNode>(0);
    auto n2    = std::make_unique<TreeNode>(1);
    n2->parent = n1.get();

    auto* n1w = n1.get();
    auto* n2w = n2.get();

    EXPECT_TRUE(TreeNode::Find(n1w) == TreeNode::Find(n2w));
}

TEST(TreeNodeTest, TestFindSameTreeParent) {
    auto n1    = std::make_unique<TreeNode>(0);
    auto n2    = std::make_unique<TreeNode>(1);
    n1->parent = n2.get();
    auto* n1w  = n1.get();
    auto* n2w  = n2.get();

    EXPECT_TRUE(TreeNode::Find(n1w) == TreeNode::Find(n2w));
}

TEST(TreeNodeTest, TestFindThreeNodesCheckParentChild) {
    auto n1    = std::make_unique<TreeNode>(0);
    auto n2    = std::make_unique<TreeNode>(1);
    auto n3    = std::make_unique<TreeNode>(2);
    n3->parent = n2.get();
    n2->parent = n1.get();
    auto* n1w  = n1.get();
    auto* n2w  = n2.get();

    EXPECT_TRUE(TreeNode::Find(n1w) == TreeNode::Find(n2w));
}
TEST(TreeNodeTest, TestFindThreeNodesCheckParentGrandchild) {
    auto n1    = std::make_unique<TreeNode>(0);
    auto n2    = std::make_unique<TreeNode>(1);
    auto n3    = std::make_unique<TreeNode>(2);
    n3->parent = n2.get();
    n2->parent = n1.get();
    auto* n1w  = n1.get();
    auto* n3w  = n3.get();

    EXPECT_TRUE(TreeNode::Find(n1w) == TreeNode::Find(n3w));
}
TEST(TreeNodeTest, TestFindThreeNodesCheckChildGrandchild) {
    auto n1    = std::make_unique<TreeNode>(0);
    auto n2    = std::make_unique<TreeNode>(1);
    auto n3    = std::make_unique<TreeNode>(2);
    n3->parent = n2.get();
    n2->parent = n1.get();
    auto* n2w  = n1.get();
    auto* n3w  = n1.get();

    EXPECT_TRUE(TreeNode::Find(n2w) == TreeNode::Find(n3w));
}

TEST(TreeNodeTest, TestFindFourNodesCheckChildGrandchild) {
    auto n1    = std::make_unique<TreeNode>(0);
    auto n2    = std::make_unique<TreeNode>(1);
    auto n3    = std::make_unique<TreeNode>(2);
    auto n4    = std::make_unique<TreeNode>(2);
    n4->parent = n3.get();
    n3->parent = n2.get();
    n2->parent = n1.get();
    auto* n2w  = n2.get();
    auto* n4w  = n4.get();

    EXPECT_TRUE(TreeNode::Find(n2w) == TreeNode::Find(n4w));
}

TEST(TreeNodeTest, TestFindFourNodesCheckSecondChildGrandchild) {
    auto n1    = std::make_unique<TreeNode>(0);
    auto n2    = std::make_unique<TreeNode>(1);
    auto n3    = std::make_unique<TreeNode>(2);
    auto n4    = std::make_unique<TreeNode>(2);
    n4->parent = n3.get();
    n3->parent = n2.get();
    n2->parent = n1.get();
    auto* n3w  = n3.get();
    auto* n4w  = n4.get();

    EXPECT_TRUE(TreeNode::Find(n3w) == TreeNode::Find(n4w));
}

TEST(TreeNodeTest, TestFindTwoNodesCheckSecondChildGrandchild) {
    auto n1    = std::make_unique<TreeNode>(0);
    auto n2    = std::make_unique<TreeNode>(1);
    n2->parent = n1.get();
    auto* n1w  = n1.get();
    auto* n2w  = n2.get();
    EXPECT_TRUE(TreeNode::Find(n1w) == TreeNode::Find(n2w));
}

TEST(TreeNodeTest, TestFindMultiBranchCheckChildGrandchild) {
    auto n1 = std::make_unique<TreeNode>(0);
    auto n2 = std::make_unique<TreeNode>(1);
    auto n3 = std::make_unique<TreeNode>(2);
    auto n4 = std::make_unique<TreeNode>(3);

    n4->parent = n2.get();
    n3->parent = n2.get();
    n2->parent = n1.get();
    auto* n3w  = n3.get();
    auto* n4w  = n4.get();

    EXPECT_TRUE(TreeNode::Find(n4w) == TreeNode::Find(n3w));
}

TEST(TreeNodeTest, TestUnion) {
    auto  n1  = std::make_unique<TreeNode>(0);
    auto  n2  = std::make_unique<TreeNode>(1);
    auto* n1w = n1.get();
    auto* n2w = n2.get();

    TreeNode::Union(n1w, n2w);

    EXPECT_TRUE(TreeNode::Find(n1w) == TreeNode::Find(n2w));
}

TEST(TreeNodeTest, TestUnionThreeNodes) {
    auto  n1   = std::make_unique<TreeNode>(0);
    auto  n2   = std::make_unique<TreeNode>(1);
    auto  n3   = std::make_unique<TreeNode>(2);
    auto  n4   = std::make_unique<TreeNode>(3);
    auto* n3w  = n3.get();
    auto* n4w  = n4.get();
    n3->parent = n2.get();
    n2->parent = n1.get();

    TreeNode::Union(n3w, n4w);
    EXPECT_TRUE(TreeNode::Find(n3w) == TreeNode::Find(n4w));
}
