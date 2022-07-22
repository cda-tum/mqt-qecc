//
// Created by lucas on 13/06/22.
//

#include "Codes.hpp"
#include "OriginalUFD.hpp"

#include <gtest/gtest.h>

class UtilsTest: public testing::TestWithParam<std::string> {
protected:
    void setUp() {
    }
};

TEST(UtilsTest, MatConversion) {
    auto ctxx = flint::nmodxx_ctx(2);

    gf2Mat matrix = {{1, 1, 0, 1, 0, 0, 1},
                     {1, 0, 1, 0, 1, 0, 0},
                     {0, 1, 1, 0, 0, 1, 0}};
    auto   sol    = flint::nmod_matxx(matrix.size(), matrix.at(0).size(), 2);
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix.at(0).size(); j++) {
            sol.at(i, j) = flint::nmodxx::red(matrix[i][j], ctxx);
        }
    }

    auto res = Utils::getFlintMatrix(matrix);
    EXPECT_TRUE(sol == res);
}

TEST(UtilsTest, MatConversionBack) {
    auto ctxx = flint::nmodxx_ctx(2);

    gf2Mat matrix = {{1, 1, 0, 1, 0, 0, 1},
                     {1, 0, 1, 0, 1, 0, 0},
                     {0, 1, 1, 0, 0, 1, 0}};
    auto   sol    = flint::nmod_matxx(matrix.size(), matrix.at(0).size(), 2);
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix.at(0).size(); j++) {
            sol.at(i, j) = flint::nmodxx::red(matrix[i][j], ctxx);
        }
    }
    auto res = Utils::getMatrixFromFlint(sol);
    EXPECT_TRUE(res == matrix);
}

TEST(UtilsTest, TestSwapRows) {
    gf2Mat matrix = {{1, 1, 0, 1, 0, 0, 1},
                     {1, 0, 1, 0, 1, 0, 0},
                     {0, 1, 1, 0, 0, 1, 0}};
    gf2Mat sol    = {
               {1, 1, 0, 1, 0, 0, 1},
               {0, 1, 1, 0, 0, 1, 0},
               {1, 0, 1, 0, 1, 0, 0}};

    Utils::swapRows(matrix, 1, 2);
    EXPECT_TRUE(sol == matrix);
}

TEST(UtilsTest, TestReduce) {
    gf2Mat matrix = {{1, 1, 0, 1, 0, 0, 1},
                     {1, 0, 1, 0, 1, 0, 0},
                     {0, 1, 1, 0, 0, 1, 0}};
    gf2Mat sol    = {
               {1, 0, 1, 0, 1, 0, 0},
               {0, 1, 1, 0, 0, 1, 0},
               {0, 0, 0, 1, 1, 1, 1}};
    auto res = Utils::gauss(matrix);
    EXPECT_TRUE(Utils::getFlintMatrix(sol) == res);
}

TEST(UtilsTest, TestTranspose) {
    gf2Mat matrix  = {{1, 0, 0, 1, 0, 1, 1},
                      {0, 1, 0, 1, 1, 0, 1},
                      {0, 0, 1, 0, 1, 1, 1}};
    gf2Mat matrixT = {
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
    gf2Mat matrix  = {{1, 0, 0, 1, 0, 1, 1},
                      {0, 1, 0, 1, 1, 0, 1},
                      {0, 0, 1, 0, 1, 1, 1}};
    gf2Mat matrixT = {
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
    gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                     {0, 1, 0, 1, 1, 0, 1},
                     {0, 0, 1, 0, 1, 1, 1}};

    gf2Vec vector = {0, 0, 0, 1, 0, 1, 1};

    auto isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_FALSE(isInRowSpace);
}

TEST(UtilsTest, GaussGF2testNotInRS2) {
    gf2Mat matrix2 = {{1, 0, 0, 0},
                      {0, 1, 0, 0},
                      {0, 0, 1, 0}};
    gf2Vec vector2 = {1, 1, 1, 1};

    auto isInRowSpace2 = Utils::isVectorInRowspace(matrix2, vector2);
    EXPECT_FALSE(isInRowSpace2);
}

TEST(UtilsTest, GaussGF2testInRS) {
    gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                     {0, 1, 0, 1, 1, 0, 1},
                     {0, 0, 1, 0, 1, 1, 1}};

    gf2Vec vector       = {1, 0, 0, 1, 0, 1, 1};
    auto   isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_TRUE(isInRowSpace);
}

TEST(UtilsTest, GaussGF2testInRS2) {
    gf2Mat matrix2       = {{1, 0, 0, 0},
                            {0, 1, 0, 0},
                            {0, 0, 1, 0}};
    gf2Vec vector2       = {1, 1, 1, 0};
    auto   isInRowSpace2 = Utils::isVectorInRowspace(matrix2, vector2);
    EXPECT_TRUE(isInRowSpace2);
}

TEST(UtilsTest, GaussGF2testNotInRSTrivial) {
    gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                     {0, 1, 0, 1, 1, 0, 1},
                     {0, 0, 1, 0, 1, 1, 1}};

    gf2Vec vector       = {1, 1, 1, 1, 1, 1, 1};
    auto   isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_FALSE(isInRowSpace);
}

TEST(UtilsTest, GaussGF2testInRSTrivial) {
    gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                     {0, 1, 0, 1, 1, 0, 1},
                     {0, 0, 1, 0, 1, 1, 1}};

    gf2Vec vector       = {0, 0, 0, 0, 0, 0, 0};
    auto   isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_TRUE(isInRowSpace);
}

TEST(UtilsTest, LinEqSolvTest) {
    gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                     {0, 1, 0, 1, 1, 0, 1},
                     {0, 0, 1, 0, 1, 1, 1}};

    gf2Vec vector   = {0, 1, 0};
    auto   res      = Utils::solveSystem(matrix, vector);
    gf2Vec solution = {0, 1, 0, 0, 0, 0, 0};
    std::cout << "res: ";
    Utils::printGF2vector(res);
    EXPECT_TRUE(res == solution);
}

TEST(UtilsTest, LinEqSolvTest2) {
    gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                     {0, 1, 0, 1, 1, 0, 1},
                     {0, 0, 1, 0, 1, 1, 1}};

    gf2Vec vector   = {0, 1, 1};
    auto   res      = Utils::solveSystem(matrix, vector);
    gf2Vec solution = {0, 1, 1, 0, 0, 0, 0};
    std::cout << "sol";
    Utils::printGF2vector(res);
    EXPECT_TRUE(res == solution);
}
TEST(UtilsTest, ImportCode) {
    gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                     {0, 1, 0, 1, 1, 0, 1},
                     {0, 0, 1, 0, 1, 1, 1}};

    auto res = Utils::importGf2MatrixFromFile("./resources/codes/testCode.txt");
    EXPECT_TRUE(res == matrix);
}

TEST(UtilsTest, GetFlintMatrix) {
    gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                     {0, 1, 0, 1, 1, 0, 1},
                     {0, 0, 1, 0, 1, 1, 1}};
    auto s = Utils::getFlintMatrix(matrix);
    print_pretty(s);
}

TEST(UtilsTest, MatrixMultiply) {
    gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                     {0, 1, 0, 1, 1, 0, 1},
                     {0, 0, 1, 0, 1, 1, 1}};
    gf2Vec vec = {1,0,0,0,0,0,0};
    gf2Vec sol = {1,0,0};
    EXPECT_TRUE(Utils::rectMatrixMultiply(matrix,vec) == sol);
}

TEST(UtilsTest, MatrixMultiply2) {
    gf2Mat matrix = {{1, 0, 0, 1, 0, 1, 1},
                     {0, 1, 0, 1, 1, 0, 1},
                     {0, 0, 1, 0, 1, 1, 1}};
    gf2Vec vec = {1,1,0,0,0,0,0};
    gf2Vec sol = {1,1,0};
    EXPECT_TRUE(Utils::rectMatrixMultiply(matrix,vec) == sol);
}