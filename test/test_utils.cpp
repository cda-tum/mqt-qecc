//
// Created by luca on 13/06/22.
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

    std::vector<std::vector<bool>> matrix = {{1, 1, 0, 1, 0, 0, 1},
                                             {1, 0, 1, 0, 1, 0, 0},
                                             {0, 1, 1, 0, 0, 1, 0}};
    auto                           sol    = flint::nmod_matxx(matrix.size(), matrix.at(0).size(), 2);
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

    std::vector<std::vector<bool>> matrix = {{1, 1, 0, 1, 0, 0, 1},
                                             {1, 0, 1, 0, 1, 0, 0},
                                             {0, 1, 1, 0, 0, 1, 0}};
    auto                           sol    = flint::nmod_matxx(matrix.size(), matrix.at(0).size(), 2);
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix.at(0).size(); j++) {
            sol.at(i, j) = flint::nmodxx::red(matrix[i][j], ctxx);
        }
    }
    auto res = Utils::getMatrixFromFlint(sol);
    EXPECT_TRUE(res == matrix);
}

TEST(UtilsTest, TestSwapRows) {
    std::vector<std::vector<bool>> matrix = {{1, 1, 0, 1, 0, 0, 1},
                                             {1, 0, 1, 0, 1, 0, 0},
                                             {0, 1, 1, 0, 0, 1, 0}};
    std::vector<std::vector<bool>> sol    = {
               {1, 1, 0, 1, 0, 0, 1},
               {0, 1, 1, 0, 0, 1, 0},
               {1, 0, 1, 0, 1, 0, 0}};

    Utils::swapRows(matrix, 1, 2);
    EXPECT_TRUE(sol == matrix);
}

TEST(UtilsTest, TestReduce) {
    std::vector<std::vector<bool>> matrix = {{1, 1, 0, 1, 0, 0, 1},
                                             {1, 0, 1, 0, 1, 0, 0},
                                             {0, 1, 1, 0, 0, 1, 0}};
    std::vector<std::vector<bool>> sol    = {
               {1, 0, 1, 0, 1, 0, 0},
               {0, 1, 1, 0, 0, 1, 0},
               {0, 0, 0, 1, 1, 1, 1}};
    auto res = Utils::gauss(matrix);
    std::cout << "res:" << std::endl;
    Utils::printGF2matrix(res);
    std::cout << "sol" << std::endl;
    Utils::printGF2matrix(sol);
    EXPECT_TRUE(sol == res);
}

TEST(UtilsTest, TestTranspose) {
    std::vector<std::vector<bool>> matrix  = {{1, 0, 0, 1, 0, 1, 1},
                                              {0, 1, 0, 1, 1, 0, 1},
                                              {0, 0, 1, 0, 1, 1, 1}};
    std::vector<std::vector<bool>> matrixT = {
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
    std::vector<std::vector<bool>> matrix  = {{1, 0, 0, 1, 0, 1, 1},
                                              {0, 1, 0, 1, 1, 0, 1},
                                              {0, 0, 1, 0, 1, 1, 1}};
    std::vector<std::vector<bool>> matrixT = {
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
    std::vector<std::vector<bool>> matrix = {{1, 0, 0, 1, 0, 1, 1},
                                             {0, 1, 0, 1, 1, 0, 1},
                                             {0, 0, 1, 0, 1, 1, 1}};

    std::vector<bool> vector = {0, 0, 0, 1, 0, 1, 1};

    auto isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_FALSE(isInRowSpace);
};

TEST(UtilsTest, GaussGF2testNotInRS2) {
    std::vector<std::vector<bool>> matrix2 = {{1, 0, 0, 0},
                                              {0, 1, 0, 0},
                                              {0, 0, 1, 0}};
    std::vector<bool>              vector2 = {1, 1, 1, 1};

    auto isInRowSpace2 = Utils::isVectorInRowspace(matrix2, vector2);
    EXPECT_FALSE(isInRowSpace2);
};

TEST(UtilsTest, GaussGF2testInRS) {
    std::vector<std::vector<bool>> matrix = {{1, 0, 0, 1, 0, 1, 1},
                                             {0, 1, 0, 1, 1, 0, 1},
                                             {0, 0, 1, 0, 1, 1, 1}};

    std::vector<bool> vector       = {1, 0, 0, 1, 0, 1, 1};
    auto              isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_TRUE(isInRowSpace);
};

TEST(UtilsTest, GaussGF2testInRS2) {
    std::vector<std::vector<bool>> matrix2       = {{1, 0, 0, 0},
                                                    {0, 1, 0, 0},
                                                    {0, 0, 1, 0}};
    std::vector<bool>              vector2       = {1, 1, 1, 0};
    auto                           isInRowSpace2 = Utils::isVectorInRowspace(matrix2, vector2);
    EXPECT_TRUE(isInRowSpace2);
};

TEST(UtilsTest, GaussGF2testNotInRSTrivial) {
    std::vector<std::vector<bool>> matrix = {{1, 0, 0, 1, 0, 1, 1},
                                             {0, 1, 0, 1, 1, 0, 1},
                                             {0, 0, 1, 0, 1, 1, 1}};

    std::vector<bool> vector       = {1, 1, 1, 1, 1, 1, 1};
    auto              isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_FALSE(isInRowSpace);
};

TEST(UtilsTest, GaussGF2testInRSTrivial) {
    std::vector<std::vector<bool>> matrix = {{1, 0, 0, 1, 0, 1, 1},
                                             {0, 1, 0, 1, 1, 0, 1},
                                             {0, 0, 1, 0, 1, 1, 1}};

    std::vector<bool> vector       = {0, 0, 0, 0, 0, 0, 0};
    auto              isInRowSpace = Utils::isVectorInRowspace(matrix, vector);
    EXPECT_TRUE(isInRowSpace);
};

TEST(UtilsTest, LinEqSolvTest) {
    std::vector<std::vector<bool>> matrix = {{1, 0, 0, 1, 0, 1, 1},
                                             {0, 1, 0, 1, 1, 0, 1},
                                             {0, 0, 1, 0, 1, 1, 1}};

    std::vector<bool> vector   = {0, 1, 0};
    auto              res      = Utils::solveSystem(matrix, vector);
    std::vector<bool> solution = {0, 1, 0, 0, 0, 0, 0};
    std::cout << "res: ";
    Utils::printGF2vector(res);
    EXPECT_TRUE(res == solution);
};

TEST(UtilsTest, LinEqSolvTest2) {
    std::vector<std::vector<bool>> matrix = {{1, 0, 0, 1, 0, 1, 1},
                                             {0, 1, 0, 1, 1, 0, 1},
                                             {0, 0, 1, 0, 1, 1, 1}};

    std::vector<bool> vector   = {0, 1, 1};
    auto              res      = Utils::solveSystem(matrix, vector);
    std::vector<bool> solution = {0, 1, 1, 0, 0, 0, 0};
    std::cout << "sol";
    Utils::printGF2vector(res);
    EXPECT_TRUE(res == solution);
};
