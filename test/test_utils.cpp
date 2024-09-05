// to keep 0/1 in boolean areas without clang-tidy warnings:
// NOLINTBEGIN(readability-implicit-bool-conversion,modernize-use-bool-literals)

#include "QeccException.hpp"
#include "Utils.hpp"

#include <gtest/gtest.h>
#include <string>

class UtilsTest : public testing::TestWithParam<std::string> {};

TEST(UtilsTest, TestSwapRows) {
    gf2Mat       matrix = {{1, 1, 0, 1, 0, 0, 1},
                           {1, 0, 1, 0, 1, 0, 0},
                           {0, 1, 1, 0, 0, 1, 0}};
    const gf2Mat sol    = {{1, 1, 0, 1, 0, 0, 1},
                           {0, 1, 1, 0, 0, 1, 0},
                           {1, 0, 1, 0, 1, 0, 0}};

    Utils::swapRows(matrix, 1, 2);
    EXPECT_TRUE(sol == matrix);
}

TEST(UtilsTest, TestTranspose) {
    const gf2Mat matrix  = {{1, 0, 0, 1, 0, 1, 1},
                            {0, 1, 0, 1, 1, 0, 1},
                            {0, 0, 1, 0, 1, 1, 1}};
    const gf2Mat matrixT = {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1},
            {1, 1, 0},
            {0, 1, 1},
            {1, 0, 1},
            {1, 1, 1}};
    EXPECT_TRUE(matrixT == Utils::getTranspose(matrix));
}

TEST(UtilsTest, TestTranspose2) {
    const gf2Mat matrix  = {{1, 0, 0, 1, 0, 1, 1},
                            {0, 1, 0, 1, 1, 0, 1},
                            {0, 0, 1, 0, 1, 1, 1}};
    const gf2Mat matrixT = {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1},
            {1, 1, 0},
            {0, 1, 1},
            {1, 0, 1},
            {1, 1, 1}};
    EXPECT_TRUE(matrix == Utils::getTranspose(matrixT));
}

TEST(UtilsTest, GaussGF2testNotInRS) {
    const gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                           {0, 1, 0, 1, 1, 0, 1},
                           {0, 0, 1, 0, 1, 1, 1}};

    const gf2Vec vector = {0, 0, 0, 1, 0, 1, 1};

    auto isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_FALSE(isInRowSpace);
}

TEST(UtilsTest, GaussGF2testNotInRS2) {
    const gf2Mat matrix2 = {{1, 0, 0, 0},
                            {0, 1, 0, 0},
                            {0, 0, 1, 0}};
    const gf2Vec vector2 = {1, 1, 1, 1};

    auto isInRowSpace2 = Utils::isVectorInRowspace(matrix2, vector2);
    EXPECT_FALSE(isInRowSpace2);
}

TEST(UtilsTest, GaussGF2testInRS) {
    const gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                           {0, 1, 0, 1, 1, 0, 1},
                           {0, 0, 1, 0, 1, 1, 1}};

    const gf2Vec vector       = {1, 0, 0, 1, 0, 1, 1};
    auto         isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_TRUE(isInRowSpace);
}

TEST(UtilsTest, GaussGF2testInRS2) {
    const gf2Mat matrix2       = {{1, 0, 0, 0},
                                  {0, 1, 0, 0},
                                  {0, 0, 1, 0}};
    const gf2Vec vector2       = {1, 1, 1, 0};
    auto         isInRowSpace2 = Utils::isVectorInRowspace(matrix2, vector2);
    EXPECT_TRUE(isInRowSpace2);
}

TEST(UtilsTest, GaussGF2testNotInRSTrivial) {
    const gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                           {0, 1, 0, 1, 1, 0, 1},
                           {0, 0, 1, 0, 1, 1, 1}};

    const gf2Vec vector       = {1, 1, 1, 1, 1, 1, 1};
    auto         isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_FALSE(isInRowSpace);
}

TEST(UtilsTest, GaussGF2testInRSTrivial) {
    const gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                           {0, 1, 0, 1, 1, 0, 1},
                           {0, 0, 1, 0, 1, 1, 1}};

    const gf2Vec vector       = {0, 0, 0, 0, 0, 0, 0};
    auto         isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_TRUE(isInRowSpace);
}

TEST(UtilsTest, ImportCode) {
    const gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                           {0, 1, 0, 1, 1, 0, 1},
                           {0, 0, 1, 0, 1, 1, 1}};

    auto res = Utils::importGf2MatrixFromFile("./resources/codes/testCode.txt");
    EXPECT_TRUE(res == matrix);
}

TEST(UtilsTest, ImportCodeFileNotExists) {
    EXPECT_THROW(Utils::importGf2MatrixFromFile("code_that_does_not_exist.txt"), QeccException);
}

TEST(UtilsTest, MatrixMultiply) {
    const gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                           {0, 1, 0, 1, 1, 0, 1},
                           {0, 0, 1, 0, 1, 1, 1}};
    const gf2Vec vec    = {1, 0, 0, 0, 0, 0, 0};
    const gf2Vec sol    = {1, 0, 0};
    gf2Vec       syndr(matrix.size());
    Utils::rectMatrixMultiply(matrix, vec, syndr);
    EXPECT_TRUE(syndr == sol);
}

TEST(UtilsTest, MatrixMultiply2) {
    const gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                           {0, 1, 0, 1, 1, 0, 1},
                           {0, 0, 1, 0, 1, 1, 1}};
    const gf2Vec vec    = {1, 1, 0, 0, 0, 0, 0};
    const gf2Vec sol    = {1, 1, 0};
    gf2Vec       comp(matrix.size());
    Utils::rectMatrixMultiply(matrix, vec, comp);

    EXPECT_TRUE(comp == sol);
}

// NOLINTEND(readability-implicit-bool-conversion,modernize-use-bool-literals)
