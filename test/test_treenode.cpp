//
// Created by lucas on 13/06/22.
//

#include "Codes.hpp"
#include "OriginalUFD.hpp"

#include <gtest/gtest.h>

class TreeNodeTest: public testing::TestWithParam<std::string> {
protected:
    void setUp() {
    }
};

TEST(TreeNodeTest, TestFindDisjoint) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);


    EXPECT_FALSE(TreeNode::Find(n1) == TreeNode::Find(n2));
}


TEST(TreeNodeTest, TestFindSameTreeChild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    n2->parent = n1;

    EXPECT_TRUE(TreeNode::Find(n1) == TreeNode::Find(n2));
}

TEST(TreeNodeTest, TestFindSameTreeParent) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    n1->parent = n2;

    EXPECT_TRUE(TreeNode::Find(n1) == TreeNode::Find(n2));
}

TEST(TreeNodeTest, TestFindThreeNodesCheckParentChild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    n3->parent = n2;
    n2->parent = n1;

    EXPECT_TRUE(TreeNode::Find(n1) == TreeNode::Find(n2));
}
TEST(TreeNodeTest, TestFindThreeNodesCheckParentGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    n3->parent = n2;
    n2->parent = n1;

    EXPECT_TRUE(TreeNode::Find(n1) == TreeNode::Find(n3));
}
TEST(TreeNodeTest, TestFindThreeNodesCheckChildGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    n3->parent = n2;
    n2->parent = n1;

    EXPECT_TRUE(TreeNode::Find(n2) == TreeNode::Find(n3));
}

TEST(TreeNodeTest, TestFindFourNodesCheckChildGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    auto n4 = std::make_shared<TreeNode>(2);
    n4->parent = n3;
    n3->parent = n2;
    n2->parent = n1;

    EXPECT_TRUE(TreeNode::Find(n2) == TreeNode::Find(n4));
}

TEST(TreeNodeTest, TestFindFourNodesCheckSecondChildGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    auto n4 = std::make_shared<TreeNode>(2);
    n4->parent = n3;
    n3->parent = n2;
    n2->parent = n1;

    EXPECT_TRUE(TreeNode::Find(n3) == TreeNode::Find(n4));
}

TEST(TreeNodeTest, TestFindTwoNodesCheckSecondChildGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    n2->parent = n1;

    EXPECT_TRUE(TreeNode::Find(n1) == TreeNode::Find(n2));
}

TEST(TreeNodeTest, TestFindMultiBranchCheckChildGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    auto n4 = std::make_shared<TreeNode>(3);

    n4->parent = n2;
    n3->parent = n2;
    n2->parent = n1;

    EXPECT_TRUE(TreeNode::Find(n4) == TreeNode::Find(n3));
}


TEST(TreeNodeTest, TestUnion) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);

    TreeNode::Union(n1,n2);

    EXPECT_TRUE(TreeNode::Find(n1) == TreeNode::Find(n2));
}

TEST(TreeNodeTest, TestUnionThreeNodes) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    auto n4 = std::make_shared<TreeNode>(3);

    n3->parent = n2;
    n2->parent = n1;

    TreeNode::Union(n3,n4);
    EXPECT_TRUE(TreeNode::Find(n3) == TreeNode::Find(n4));
}






