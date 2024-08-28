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
    CscMatrix            csc_mat; // csc_mat[column][row]
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
    bool                                  LU_constructed = false;

    PluDecomposition(std::size_t rowCount, std::size_t col_count, std::vector<std::vector<std::size_t>>& cscMat)
        : csc_mat(cscMat),
          rowCount(rowCount),
          colCount(col_count) {
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

        this->LU_constructed = false;
    }

    /**
     * Compute the reduced row-echelon form
     * The algorithm operates column, wise. The original matrix, csc_mat is not modified.
     * Instead, rr_col is used to represent the column currently eliminated and row operations, e.g., swaps
     * and addition is done column per column.
     * First, all previous row operations are applied to the current column.
     * Then pivoting is done for the current column and corresponding swaps are applied.
     * @param construct_U
     */
    void rref(const bool construct_L = true, const bool construct_U = true) {
        this->reset();
        for (std::size_t i = 0; i < this->rowCount; i++) {
            this->rows.push_back(i);
        }

        auto max_rank = std::min(this->rowCount, this->colCount);

        if (construct_L) {
            this->lower.resize(this->rowCount, std::vector<std::size_t>{});
        }

        for (std::size_t col_idx = 0; col_idx < this->colCount; col_idx++) {
            this->eliminate_column(col_idx, construct_L, construct_U);
            if (this->matrixRank == max_rank) {
                break;
            }
        }

        if (construct_L && construct_U) {
            this->LU_constructed = true;
        }
    }

    bool eliminate_column(std::size_t col_idx, const bool construct_L = true, const bool construct_U = true) {
        auto rr_col = std::vector<uint8_t>(this->rowCount, 0);

        for (auto row_index : this->csc_mat[col_idx]) {
            rr_col[row_index] = 1;
        }
        // apply previous operations to current column
        for (std::size_t i = 0; i < this->matrixRank; i++) {
            std::swap(rr_col[i], rr_col[this->swapRows[i]]);
            if (rr_col[i] == 1) {
                // if row elem is one, do elimination for current column below the pivot
                // elimination operations to apply are stored in the `elimination_rows` attribute,
                for (std::size_t const row_idx : this->eliminationRows[i]) {
                    rr_col[row_idx] ^= 1;
                }
            }
        }
        bool PIVOT_FOUND = false;
        for (auto i = this->matrixRank; i < this->rowCount; i++) {
            if (rr_col[i] == 1) {
                PIVOT_FOUND = true;
                this->swapRows.push_back(i);
                this->pivotCols.push_back(col_idx);
                break;
            }
        }
        // if no pivot was found, we go to next column
        if (!PIVOT_FOUND) {
            this->notPivotCols.push_back(col_idx);
            return false;
        }

        std::swap(rr_col[this->matrixRank], rr_col[this->swapRows[this->matrixRank]]);
        std::swap(this->rows[this->matrixRank], this->rows[this->swapRows[this->matrixRank]]);
        this->eliminationRows.emplace_back();

        if (construct_L) {
            std::swap(this->lower[this->matrixRank], this->lower[this->swapRows[this->matrixRank]]);
            this->lower[this->matrixRank].push_back(this->matrixRank);
        }

        for (auto i = this->matrixRank + 1; i < this->rowCount; i++) {
            if (rr_col[i] == 1) {
                this->eliminationRows[this->matrixRank].push_back(i);
                if (construct_L) {
                    this->lower[i].push_back(this->matrixRank);
                }
            }
        }

        if (construct_U) {
            this->upper.emplace_back();
            for (std::size_t i = 0; i <= this->matrixRank; i++) {
                if (rr_col[i] == 1) {
                    this->upper[i].push_back(col_idx);
                }
            }
        }

        this->matrixRank++;
        return true;
    }

    std::vector<std::size_t> lu_solve(std::vector<std::size_t>& y) {
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
        if (!this->LU_constructed) {
            throw std::invalid_argument("LU decomposition has not been constructed. Please call rref() first.");
        }
        auto x = std::vector<std::size_t>(this->colCount, 0);
        auto b = std::vector<std::size_t>(this->matrixRank, 0);
        // First we solve Lb = y, where b = Ux
        // Solve Lb=y with forwared substitution
        for (std::size_t row_index = 0; row_index < this->matrixRank; row_index++) {
            std::size_t row_sum = 0;
            for (auto col_index : this->lower[row_index]) {
                row_sum ^= b[col_index];
            }
            b[row_index] = row_sum ^ y[this->rows[row_index]];
        }
        // Solve Ux = b with backwards substitution
        for (int row_index = static_cast<int>(this->matrixRank) - 1; row_index >= 0; row_index--) {
            std::size_t row_sum = 0;
            for (std::size_t const col_index : this->upper[static_cast<uint64_t>(row_index)]) {
                row_sum ^= x[col_index];
            }
            x[this->pivotCols[static_cast<uint64_t>(row_index)]] = row_sum ^ b[static_cast<uint64_t>(row_index)];
        }
        return x;
    }
};
// namespace ldpc::gf2dense
} // namespace ldpc::gf2dense

#endif
