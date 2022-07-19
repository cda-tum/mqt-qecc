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

    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);


    EXPECT_FALSE(TreeNode::Find(n1w).lock() == TreeNode::Find(n2w).lock());
}


TEST(TreeNodeTest, TestFindSameTreeChild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    n2->parent = n1;
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);

    EXPECT_TRUE(TreeNode::Find(n1w).lock() == TreeNode::Find(n2w).lock());
}

TEST(TreeNodeTest, TestFindSameTreeParent) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    n1->parent = n2;
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);


    EXPECT_TRUE(TreeNode::Find(n1w).lock() == TreeNode::Find(n2w).lock());
}

TEST(TreeNodeTest, TestFindThreeNodesCheckParentChild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    n3->parent = n2;
    n2->parent = n1;
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);
    auto n3w = std::weak_ptr<TreeNode>(n3);


    EXPECT_TRUE(TreeNode::Find(n1w).lock() == TreeNode::Find(n2w).lock());
}
TEST(TreeNodeTest, TestFindThreeNodesCheckParentGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    n3->parent = n2;
    n2->parent = n1;
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);
    auto n3w = std::weak_ptr<TreeNode>(n3);

    EXPECT_TRUE(TreeNode::Find(n1w).lock() == TreeNode::Find(n3w).lock());
}
TEST(TreeNodeTest, TestFindThreeNodesCheckChildGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    n3->parent = n2;
    n2->parent = n1;
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n1);
    auto n3w = std::weak_ptr<TreeNode>(n1);


    EXPECT_TRUE(TreeNode::Find(n2w).lock() == TreeNode::Find(n3w).lock());
}

TEST(TreeNodeTest, TestFindFourNodesCheckChildGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    auto n4 = std::make_shared<TreeNode>(2);
    n4->parent = n3;
    n3->parent = n2;
    n2->parent = n1;
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);
    auto n3w = std::weak_ptr<TreeNode>(n3);
    auto n4w = std::weak_ptr<TreeNode>(n4);


    EXPECT_TRUE(TreeNode::Find(n2w).lock() == TreeNode::Find(n4w).lock());
}

TEST(TreeNodeTest, TestFindFourNodesCheckSecondChildGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    auto n4 = std::make_shared<TreeNode>(2);
    n4->parent = n3;
    n3->parent = n2;
    n2->parent = n1;
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);
    auto n3w = std::weak_ptr<TreeNode>(n3);
    auto n4w = std::weak_ptr<TreeNode>(n4);

    EXPECT_TRUE(TreeNode::Find(n3w).lock() == TreeNode::Find(n4w).lock());
}

TEST(TreeNodeTest, TestFindTwoNodesCheckSecondChildGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    n2->parent = n1;
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);
    EXPECT_TRUE(TreeNode::Find(n1w).lock() == TreeNode::Find(n2w).lock());
}

TEST(TreeNodeTest, TestFindMultiBranchCheckChildGrandchild) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    auto n4 = std::make_shared<TreeNode>(3);

    n4->parent = n2;
    n3->parent = n2;
    n2->parent = n1;
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);
    auto n3w = std::weak_ptr<TreeNode>(n3);
    auto n4w = std::weak_ptr<TreeNode>(n4);

    EXPECT_TRUE(TreeNode::Find(n4w).lock() == TreeNode::Find(n3w).lock());
}


TEST(TreeNodeTest, TestUnion) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);

    TreeNode::Union(n1w,n2w);

    EXPECT_TRUE(TreeNode::Find(n1w).lock() == TreeNode::Find(n2w).lock());
}

TEST(TreeNodeTest, TestUnionThreeNodes) {
    auto n1 = std::make_shared<TreeNode>(0);
    auto n2 = std::make_shared<TreeNode>(1);
    auto n3 = std::make_shared<TreeNode>(2);
    auto n4 = std::make_shared<TreeNode>(3);
    auto n1w = std::weak_ptr<TreeNode>(n1);
    auto n2w = std::weak_ptr<TreeNode>(n2);
    auto n3w = std::weak_ptr<TreeNode>(n3);
    auto n4w = std::weak_ptr<TreeNode>(n4);
    n3->parent = n2;
    n2->parent = n1;

    TreeNode::Union(n3w,n4w);
    EXPECT_TRUE(TreeNode::Find(n3w).lock() == TreeNode::Find(n4w).lock());
}






