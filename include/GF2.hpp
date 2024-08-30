/// Originally taken from the LDPC library by quantumgizmos released under the MIT license.
/// https://github.com/quantumgizmos/ldpc_v2/blob/cb85ee8601ee3fe59482985ab20840525c452a98/src_cpp/gf2dense.hpp

#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

using CscMatrix = std::vector<std::vector<uint64_t>>;
using CsrMatrix = std::vector<std::vector<uint64_t>>;

/**
 * A class to represent the PLU decomposition.
 * That is for a matrix A, we have perm*A = lower*upper, where perm is a permutation matrix.
 */
class PluDecomposition {
private:
    CscMatrix cscMat; // csc_mat[column][row]

    /**
     * @brief Compute the reduced row-echelon form
     *
     * @details The algorithm operates column-wise. The original matrix is not modified.
     * Instead, rrCol is used to represent the column currently being eliminated and row operations,
     * e.g., swaps and addition are performed column by column.
     * First, all previous row operations are applied to the current column.
     * Then pivoting is performed for the current column and corresponding swaps are applied.
     */
    void rref();

    bool eliminateColumn(std::size_t colIdx);

    CsrMatrix                             lower;
    CsrMatrix                             upper;
    std::size_t                           matrixRank{};
    std::size_t                           rowCount{};
    std::size_t                           colCount{};
    std::vector<std::size_t>              rows;
    std::vector<std::size_t>              swapRows;
    std::vector<std::vector<std::size_t>> eliminationRows;
    std::vector<std::size_t>              pivotCols;

public:
    PluDecomposition(const std::size_t nrows, const std::size_t ncols, CscMatrix mat)
        : cscMat(std::move(mat)),
          rowCount(nrows),
          colCount(ncols) {
        rref();
    }

    PluDecomposition() = default;

    [[nodiscard]] std::size_t getMatrixRank() const {
        return matrixRank;
    }

    std::vector<uint64_t> luSolve(const std::vector<uint64_t>& y);
};
