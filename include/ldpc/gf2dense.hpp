#ifndef GF2DENSE_H
#define GF2DENSE_H

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ldpc::gf2dense {
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

    PluDecomposition(std::size_t rowCount, std::size_t colCount, std::vector<std::vector<std::size_t>>& cscMat)
        : cscMat(cscMat),
          rowCount(rowCount),
          colCount(colCount) {
    }

    PluDecomposition() = default;

    ~PluDecomposition() =
            default;

    /**
     * Reset all internal information.
     */
    void reset() {
        this->matrixRank = 0;
        this->rows.clear();
        this->swapRows.clear();
        this->pivotCols.clear();
        this->notPivotCols.clear();
        this->yImageCheckVector.clear();

        for (auto& col : this->lower) {
            col.clear();
        }
        this->lower.clear();

        for (auto& col : this->eliminationRows) {
            col.clear();
        }
        this->eliminationRows.clear();

        for (auto& col : this->upper) {
            col.clear();
        }
        this->upper.clear();

        for (auto& col : this->perm) {
            col.clear();
        }
        this->perm.clear();

        this->luConstructed = false;
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
        this->reset();
        for (std::size_t i = 0; i < this->rowCount; i++) {
            this->rows.push_back(i);
        }

        auto maxRank = std::min(this->rowCount, this->colCount);

        if (constructL) {
            this->lower.resize(this->rowCount, std::vector<std::size_t>{});
        }

        for (std::size_t colIdx = 0; colIdx < this->colCount; colIdx++) {
            this->eliminateColumn(colIdx, constructL, constructU);
            if (this->matrixRank == maxRank) {
                break;
            }
        }

        if (constructL && constructU) {
            this->luConstructed = true;
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
        if (y.size() != this->rowCount) {
            throw std::invalid_argument("Input parameter `y` is of the incorrect length for lu_solve.");
        }
        if (!this->luConstructed) {
            throw std::invalid_argument("LU decomposition has not been constructed. Please call rref() first.");
        }
        auto x = std::vector<std::size_t>(this->colCount, 0);
        auto b = std::vector<std::size_t>(this->matrixRank, 0);
        // First we solve Lb = y, where b = Ux
        // Solve Lb=y with forwarded substitution
        for (std::size_t rowIndex = 0; rowIndex < this->matrixRank; rowIndex++) {
            std::size_t rowSum = 0;
            for (auto colIndex : this->lower[rowIndex]) {
                rowSum ^= b[colIndex];
            }
            b[rowIndex] = rowSum ^ y[this->rows[rowIndex]];
        }
        // Solve Ux = b with backwards substitution
        for (int rowIndex = static_cast<int>(this->matrixRank) - 1; rowIndex >= 0; rowIndex--) {
            std::size_t rowSum = 0;
            for (std::size_t const colIndex : this->upper[static_cast<uint64_t>(rowIndex)]) {
                rowSum ^= x[colIndex];
            }
            x[this->pivotCols[static_cast<uint64_t>(rowIndex)]] = rowSum ^ b[static_cast<uint64_t>(rowIndex)];
        }
        return x;
    }

    bool eliminateColumn(std::size_t colIdx, const bool constructL, const bool constructU) {
        auto rrCol = std::vector<uint8_t>(this->rowCount, 0);

        for (auto rowIndex : this->cscMat[colIdx]) {
            rrCol[rowIndex] = 1;
        }
        // apply previous operations to current column
        for (std::size_t i = 0; i < this->matrixRank; i++) {
            std::swap(rrCol[i], rrCol[this->swapRows[i]]);
            if (rrCol[i] == 1) {
                // if row elem is one, do elimination for current column below the pivot
                // elimination operations to apply are stored in the `elimination_rows` attribute,
                for (std::size_t const rowIdx : this->eliminationRows[i]) {
                    rrCol[rowIdx] ^= 1;
                }
            }
        }
        bool pivotFound = false;
        for (auto i = this->matrixRank; i < this->rowCount; i++) {
            if (rrCol[i] == 1) {
                pivotFound = true;
                this->swapRows.push_back(i);
                this->pivotCols.push_back(colIdx);
                break;
            }
        }
        // if no pivot was found, we go to next column
        if (!pivotFound) {
            this->notPivotCols.push_back(colIdx);
            return false;
        }

        std::swap(rrCol[this->matrixRank], rrCol[this->swapRows[this->matrixRank]]);
        std::swap(this->rows[this->matrixRank], this->rows[this->swapRows[this->matrixRank]]);
        this->eliminationRows.emplace_back();

        if (constructL) {
            std::swap(this->lower[this->matrixRank], this->lower[this->swapRows[this->matrixRank]]);
            this->lower[this->matrixRank].push_back(this->matrixRank);
        }

        for (auto i = this->matrixRank + 1; i < this->rowCount; i++) {
            if (rrCol[i] == 1) {
                this->eliminationRows[this->matrixRank].push_back(i);
                if (constructL) {
                    this->lower[i].push_back(this->matrixRank);
                }
            }
        }

        if (constructU) {
            this->upper.emplace_back();
            for (std::size_t i = 0; i <= this->matrixRank; i++) {
                if (rrCol[i] == 1) {
                    this->upper[i].push_back(colIdx);
                }
            }
        }

        this->matrixRank++;
        return true;
    }
};
} // namespace ldpc::gf2dense

#endif
