#ifndef GF2DENSE_H
#define GF2DENSE_H


#include <chrono>
#include <climits>
#include <iterator>
#include <random>
#include <vector>

namespace ldpc::gf2dense {

    template<class T>
    int vector_find(std::vector<T> vec, T value) {
        int index = 0;
        for (auto val: vec) {
            if (val == value) {
                return index;
            }
            index++;
        }
        return -1;
    }

    enum PluMethod {
        SPARSE_ELIMINATION = 0,
        DENSE_ELIMINATION = 1
    };

    typedef std::vector<std::vector<int>> CscMatrix;
    using CsrMatrix = std::vector<std::vector<int>>;


/**
 * A class to represent the PLU decomposition.
 * That is for a matrix A, we have P*A = L*U, where P is a permutation matrix.
 */
    class PluDecomposition {
    private:
        CscMatrix csc_mat; // csc_mat[column][row]
        std::vector<uint8_t> y_image_check_vector;

    public:
        CsrMatrix L;
        CsrMatrix U;
        CscMatrix P;
        int matrix_rank{};
        int cols_eliminated{};
        int row_count{};
        int col_count{};
        std::vector<int> rows;
        std::vector<int> swap_rows;
        std::vector<std::vector<int>> elimination_rows;
        std::vector<int> pivot_cols;
        std::vector<int> not_pivot_cols;
        bool LU_constructed = false;

        PluDecomposition(int row_count, int col_count, std::vector<std::vector<int>> &csc_mat)
                : row_count(row_count),
                  col_count(col_count),
                  csc_mat(csc_mat) {
        }

        PluDecomposition() = default;

        ~PluDecomposition() =
        default;

        /**
         * Reset all internal information.
         */
        void reset() {
            this->matrix_rank = 0;
            this->cols_eliminated = 0;
            this->rows.clear();
            this->swap_rows.clear();
            this->pivot_cols.clear();
            this->not_pivot_cols.clear();
            this->y_image_check_vector.clear();

            for (auto &col: this->L) {
                col.clear();
            }
            this->L.clear();

            for (auto &col: this->elimination_rows) {
                col.clear();
            }
            this->elimination_rows.clear();

            for (auto &col: this->U) {
                col.clear();
            }
            this->U.clear();

            for (auto &col: this->P) {
                col.clear();
            }
            this->P.clear();

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
            for (auto i = 0; i < this->row_count; i++) {
                this->rows.push_back(i);
            }

            auto max_rank = std::min(this->row_count, this->col_count);

            if (construct_L) {
                this->L.resize(this->row_count, std::vector<int>{});
            }

            // if (construct_U){
            //     this->U.resize(this->row_count, std::vector<int>{});
            // }

            for (auto col_idx = 0; col_idx < this->col_count; col_idx++) {
                this->eliminate_column(col_idx, construct_L, construct_U);
                if (this->matrix_rank == max_rank) {
                    break;
                }
            }

            if (construct_L && construct_U) {
                this->LU_constructed = true;
            }
        }

        bool eliminate_column(int col_idx, const bool construct_L = true, const bool construct_U = true) {
            auto rr_col = std::vector<uint8_t>(this->row_count, 0);
            this->cols_eliminated = col_idx + 1;

            for (auto row_index: this->csc_mat[col_idx]) {
                rr_col[row_index] = 1;
            }
            // apply previous operations to current column
            for (auto i = 0; i < this->matrix_rank; i++) {
                std::swap(rr_col[i], rr_col[this->swap_rows[i]]);
                if (rr_col[i] == 1) {
                    // if row elem is one, do elimination for current column below the pivot
                    // elimination operations to apply are stored in the `elimination_rows` attribute,
                    for (auto row_idx: this->elimination_rows[i]) {
                        rr_col[row_idx] ^= 1;
                    }
                }
            }
            bool PIVOT_FOUND = false;
            for (auto i = this->matrix_rank; i < this->row_count; i++) {
                if (rr_col[i] == 1) {
                    PIVOT_FOUND = true;
                    this->swap_rows.push_back(i);
                    this->pivot_cols.push_back(col_idx);
                    break;
                }
            }
            // if no pivot was found, we go to next column
            if (!PIVOT_FOUND) {
                this->not_pivot_cols.push_back(col_idx);
                return false;
            }

            std::swap(rr_col[this->matrix_rank], rr_col[this->swap_rows[this->matrix_rank]]);
            std::swap(this->rows[this->matrix_rank], this->rows[this->swap_rows[this->matrix_rank]]);
            this->elimination_rows.emplace_back();

            if (construct_L) {
                std::swap(this->L[this->matrix_rank], this->L[this->swap_rows[this->matrix_rank]]);
                this->L[this->matrix_rank].push_back(this->matrix_rank);
            }

            for (auto i = this->matrix_rank + 1; i < this->row_count; i++) {
                if (rr_col[i] == 1) {
                    this->elimination_rows[this->matrix_rank].push_back(i);
                    if (construct_L) {
                        this->L[i].push_back(this->matrix_rank);
                    }
                }
            }

            if (construct_U) {
                this->U.emplace_back();
                for (auto i = 0; i <= this->matrix_rank; i++) {
                    if (rr_col[i] == 1) {
                        this->U[i].push_back(col_idx);
                    }
                }
            }

            this->matrix_rank++;
            return true;
        }

        std::vector<uint8_t> lu_solve(std::vector<uint8_t> &y) {
            /*
            Equation: Ax = y

            We use LU decomposition to arrange the above into the form:
            LU(Qx) = PAQ^T(Qx)=Py

            We can then solve for x using forward-backward substitution:
            1. Forward substitution: Solve Lb = Py for b
            2. Backward substitution: Solve UQx = b for x
            */
            if (y.size() != this->row_count) {
                throw std::invalid_argument("Input parameter `y` is of the incorrect length for lu_solve.");
            }
            if (!this->LU_constructed) {
                throw std::invalid_argument("LU decomposition has not been constructed. Please call rref() first.");
            }
            auto x = std::vector<uint8_t>(this->col_count, 0);
            auto b = std::vector<uint8_t>(this->matrix_rank, 0);
            // First we solve Lb = y, where b = Ux
            // Solve Lb=y with forwared substitution
            for (auto row_index = 0; row_index < this->matrix_rank; row_index++) {
                int row_sum = 0;
                for (auto col_index: this->L[row_index]) {
                    row_sum ^= b[col_index];
                }
                b[row_index] = row_sum ^ y[this->rows[row_index]];
            }
            // Solve Ux = b with backwards substitution
            for (auto row_index = this->matrix_rank - 1; row_index >= 0; row_index--) {
                int row_sum = 0;
                for (auto col_index: this->U[row_index]) {
                    row_sum ^= x[col_index];
                }
                x[this->pivot_cols[row_index]] = row_sum ^ b[row_index];
            }
            return x;
        }


    };
// namespace ldpc::gf2dense
}// end namespace ldpc

#endif
