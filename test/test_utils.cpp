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
TEST(OriginalUFDtest, TestTranspose) {
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

TEST(OriginalUFDtest, TestTranspose2) {
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
    std::vector<std::vector<bool>> matrix  = {{1, 0, 0, 1, 0, 1, 1},
                                              {0, 1, 0, 1, 1, 0, 1},
                                              {0, 0, 1, 0, 1, 1, 1}};
    std::vector<std::vector<bool>> matrix2 = {{1, 0, 0, 0},
                                              {0, 1, 0, 0},
                                              {0, 0, 1, 0}};
    std::vector<bool>              vector  = {0, 0, 0, 1, 0, 1, 1};
    std::vector<bool>              vector2 = {1, 1, 1, 1};

    auto transp        = Utils::getTranspose(matrix);
    auto transp2       = Utils::getTranspose(matrix2);
    auto isInRowSpace2 = Utils::isVectorInRowspace(transp2, vector2);

    auto isInRowSpace = Utils::isVectorInRowspace(transp, vector);
    EXPECT_FALSE(isInRowSpace);
    EXPECT_FALSE(isInRowSpace2);
};

TEST(UtilsTest, GaussGF2testInRS) {
    std::vector<std::vector<bool>> matrix = {{1, 0, 0, 1, 0, 1, 1},
                                             {0, 1, 0, 1, 1, 0, 1},
                                             {0, 0, 1, 0, 1, 1, 1}};

    std::vector<std::vector<bool>> matrix2       = {{1, 0, 0, 0},
                                                    {0, 1, 0, 0},
                                                    {0, 0, 1, 0}};
    std::vector<bool>              vector        = {1, 0, 0, 1, 0, 1, 1};
    std::vector<bool>              vector2       = {1, 1, 1, 0};
    auto                           transp        = Utils::getTranspose(matrix);
    auto                           transp2       = Utils::getTranspose(matrix2);
    auto                           isInRowSpace  = Utils::isVectorInRowspace(transp, vector);
    auto                           isInRowSpace2 = Utils::isVectorInRowspace(transp2, vector2);
    EXPECT_TRUE(isInRowSpace);
    EXPECT_TRUE(isInRowSpace2);
};

TEST(UtilsTest, GaussGF2testNotInRSTrivial) {
    std::vector<std::vector<bool>> matrix = {{1, 0, 0, 1, 0, 1, 1},
                                             {0, 1, 0, 1, 1, 0, 1},
                                             {0, 0, 1, 0, 1, 1, 1}};

    std::vector<bool> vector       = {1, 1, 1, 1, 1, 1, 1};
    auto              transp       = Utils::getTranspose(matrix);
    auto              isInRowSpace = Utils::isVectorInRowspace(transp, vector);
    EXPECT_FALSE(isInRowSpace);
};

TEST(UtilsTest, GaussGF2testInRSTrivial) {
    std::vector<std::vector<bool>> matrix = {{1, 0, 0, 1, 0, 1, 1},
                                             {0, 1, 0, 1, 1, 0, 1},
                                             {0, 0, 1, 0, 1, 1, 1}};

    std::vector<bool> vector       = {0, 0, 0, 0, 0, 0, 0};
    auto              transp       = Utils::getTranspose(matrix);
    auto              isInRowSpace = Utils::isVectorInRowspace(transp, vector);
    EXPECT_TRUE(isInRowSpace);
};

TEST(UtilsTest, LinEqSolvTest) {
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