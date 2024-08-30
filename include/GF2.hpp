/// Originally taken from the LDPC library by quantumgizmos released under the MIT license.
/// https://github.com/quantumgizmos/ldpc_v2/blob/cb85ee8601ee3fe59482985ab20840525c452a98/src_cpp/gf2dense.hpp

#pragma once

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

using CscMatrix = std::vector<std::vector<std::size_t>>;
using CsrMatrix = std::vector<std::vector<std::size_t>>;

/**
 * A class to represent the PLU decomposition.
 * That is for a matrix A, we have perm*A = lower*upper, where perm is a permutation matrix.
 */
class PluDecomposition {
private:
    CscMatrix            cscMat; // csc_mat[column][row]
    std::vector<uint8_t> yImageCheckVector;

public:
    CsrMatrix                             lower;
    CsrMatrix                             upper;
    CscMatrix                             perm;
    std::size_t                           matrixRank{};
    std::size_t                           rowCount{};
    std::size_t                           colCount{};
    std::vector<std::size_t>              rows;
    std::vector<std::size_t>              swapRows;
    std::vector<std::vector<std::size_t>> eliminationRows;
    std::vector<std::size_t>              pivotCols;
    std::vector<std::size_t>              notPivotCols;
    bool                                  luConstructed = false;

    PluDecomposition(const std::size_t nrows, const std::size_t ncols, std::vector<std::vector<std::size_t>>& mat)
        : cscMat(mat),
          rowCount(nrows),
          colCount(ncols) {
    }

    PluDecomposition() = default;

    ~PluDecomposition() =
            default;

    /**
     * Reset all internal information.
     */
    void reset() {
        matrixRank = 0;
        rows.clear();
        swapRows.clear();
        pivotCols.clear();
        notPivotCols.clear();
        yImageCheckVector.clear();

        for (auto& col : lower) {
            col.clear();
        }
        lower.clear();

        for (auto& col : eliminationRows) {
            col.clear();
        }
        eliminationRows.clear();

        for (auto& col : upper) {
            col.clear();
        }
        upper.clear();

        for (auto& col : perm) {
            col.clear();
        }
        perm.clear();

        luConstructed = false;
    }

    /**
     * Compute the reduced row-echelon form
     * The algorithm operates column, wise. The original matrix, csc_mat is not modified.
     * Instead, rr_col is used to represent the column currently eliminated and row operations, e.g., swaps
     * and addition is done column per column.
     * First, all previous row operations are applied to the current column.
     * Then pivoting is done for the current column and corresponding swaps are applied.
     * @param constructU
     */
    void rref(const bool constructL = true, const bool constructU = true) {
        reset();
        for (std::size_t i = 0; i < rowCount; i++) {
            rows.push_back(i);
        }

        auto maxRank = std::min(rowCount, colCount);

        if (constructL) {
            lower.resize(rowCount, std::vector<std::size_t>{});
        }

        for (std::size_t colIdx = 0; colIdx < colCount; colIdx++) {
            eliminateColumn(colIdx, constructL, constructU);
            if (matrixRank == maxRank) {
                break;
            }
        }

        if (constructL && constructU) {
            luConstructed = true;
        }
    }

    std::vector<std::size_t> luSolve(std::vector<std::size_t>& y) {
        /*
        Equation: Ax = y

        We use LU decomposition to arrange the above into the form:
        LU(Qx) = PAQ^T(Qx)=Py

        We can then solve for x using forward-backward substitution:
        1. Forward substitution: Solve Lb = Py for b
        2. Backward substitution: Solve UQx = b for x
        */
        if (y.size() != rowCount) {
            throw std::invalid_argument("Input parameter `y` is of the incorrect length for lu_solve.");
        }
        if (!luConstructed) {
            throw std::invalid_argument("LU decomposition has not been constructed. Please call rref() first.");
        }
        auto x = std::vector<std::size_t>(colCount, 0);
        auto b = std::vector<std::size_t>(matrixRank, 0);
        // First we solve Lb = y, where b = Ux
        // Solve Lb=y with forwarded substitution
        for (std::size_t rowIndex = 0; rowIndex < matrixRank; rowIndex++) {
            std::size_t rowSum = 0;
            for (auto colIndex : lower[rowIndex]) {
                rowSum ^= b[colIndex];
            }
            b[rowIndex] = rowSum ^ y[rows[rowIndex]];
        }
        // Solve Ux = b with backwards substitution
        for (int rowIndex = static_cast<int>(matrixRank) - 1; rowIndex >= 0; rowIndex--) {
            std::size_t rowSum = 0;
            for (std::size_t const colIndex : upper[static_cast<uint64_t>(rowIndex)]) {
                rowSum ^= x[colIndex];
            }
            x[pivotCols[static_cast<uint64_t>(rowIndex)]] = rowSum ^ b[static_cast<uint64_t>(rowIndex)];
        }
        return x;
    }

    bool eliminateColumn(std::size_t colIdx, const bool constructL, const bool constructU) {
        auto rrCol = std::vector<uint8_t>(rowCount, 0);

        for (auto rowIndex : cscMat[colIdx]) {
            rrCol[rowIndex] = 1;
        }
        // apply previous operations to current column
        for (std::size_t i = 0; i < matrixRank; i++) {
            std::swap(rrCol[i], rrCol[swapRows[i]]);
            if (rrCol[i] == 1) {
                // if row elem is one, do elimination for current column below the pivot
                // elimination operations to apply are stored in the `elimination_rows` attribute,
                for (std::size_t const rowIdx : eliminationRows[i]) {
                    rrCol[rowIdx] ^= 1;
                }
            }
        }
        bool pivotFound = false;
        for (auto i = matrixRank; i < rowCount; i++) {
            if (rrCol[i] == 1) {
                pivotFound = true;
                swapRows.push_back(i);
                pivotCols.push_back(colIdx);
                break;
            }
        }
        // if no pivot was found, we go to next column
        if (!pivotFound) {
            notPivotCols.push_back(colIdx);
            return false;
        }

        std::swap(rrCol[matrixRank], rrCol[swapRows[matrixRank]]);
        std::swap(rows[matrixRank], rows[swapRows[matrixRank]]);
        eliminationRows.emplace_back();

        if (constructL) {
            std::swap(lower[matrixRank], lower[swapRows[matrixRank]]);
            lower[matrixRank].push_back(matrixRank);
        }

        for (auto i = matrixRank + 1; i < rowCount; i++) {
            if (rrCol[i] == 1) {
                eliminationRows[matrixRank].push_back(i);
                if (constructL) {
                    lower[i].push_back(matrixRank);
                }
            }
        }

        if (constructU) {
            upper.emplace_back();
            for (std::size_t i = 0; i <= matrixRank; i++) {
                if (rrCol[i] == 1) {
                    upper[i].push_back(colIdx);
                }
            }
        }

        matrixRank++;
        return true;
    }
};
