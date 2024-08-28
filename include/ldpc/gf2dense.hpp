#ifndef GF2DENSE_H
#define GF2DENSE_H

#include <vector>
#include <iterator>
#include <chrono>
#include <climits>
#include <random>

#include "ldpc/util.hpp"
#include "ldpc/gf2sparse.hpp"
#include "ldpc/sparse_matrix_util.hpp"
#include "ldpc/rng.hpp"


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
    using CsrMatrix = std::vector<std::vector<int> >;


    inline CsrMatrix csc_to_csr(CscMatrix csc_mat) {
        int row_count = -1;
        for (auto &col: csc_mat) {
            for (int entry: col) {
                if (entry > row_count) {
                    row_count = entry;
                }
            }
        }

        auto csr_mat = CsrMatrix(row_count + 1, std::vector<int>{});

        for (int col_index = 0; col_index < csc_mat.size(); col_index++) {
            for (int row_index: csc_mat[col_index]) {
                csr_mat[row_index].push_back(col_index);
            }
        }
        return csr_mat;
    }


    inline void print_csr(const std::vector<std::vector<int>> &csr_mat) {

        int col_count = -1;
        for (const auto &row: csr_mat) {
            for (int entry: row) {
                if (entry > col_count) {
                    col_count = entry;
                }
            }
        }

        col_count++;

        for (const auto &row: csr_mat) {
            auto row_dense = std::vector<int>(col_count, 0);
            for (auto entry: row) {
                row_dense[entry] = 1;
            }
            for (auto entry: row_dense) {
                std::cout << unsigned(entry) << " ";
            }
            std::cout << std::endl;
        }

    }


    inline void print_csc(const std::vector<std::vector<int>> &csc_mat) {
        CsrMatrix csr_matrix;
        int col_index = 0;
        for (const auto &col: csc_mat) {
            for (auto entry: col) {
                if (entry >= csr_matrix.size()) {
                    csr_matrix.resize(entry + 1, std::vector<int>{});
                }
                csr_matrix[entry].push_back(col_index);
            }
            col_index++;
        }

        print_csr(csr_matrix);

    }

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
            2. Backward subsitution: Solve UQx = b for x
            */
            if (y.size() != this->row_count) {
                throw std::invalid_argument("Input parameter `y` is of the incorrect length for lu_solve.");
            }
            if (!this->LU_constructed) {
                throw std::invalid_argument("LU decomposition has not been constructed. Please call rref() first.");
            }
            auto x = std::vector<uint8_t>(this->col_count, 0);
            auto b = std::vector<uint8_t>(this->matrix_rank, 0);
            //First we solve Lb = y, where b = Ux
            //Solve Lb=y with forwared substitution
            for (auto row_index = 0; row_index < this->matrix_rank; row_index++) {
                int row_sum = 0;
                for (auto col_index: this->L[row_index]) {
                    row_sum ^= b[col_index];
                }
                b[row_index] = row_sum ^ y[this->rows[row_index]];
            }
            //Solve Ux = b with backwards substitution
            for (auto row_index = this->matrix_rank - 1; row_index >= 0; row_index--) {
                int row_sum = 0;
                for (auto col_index: this->U[row_index]) {
                    row_sum ^= x[col_index];
                }
                x[this->pivot_cols[row_index]] = row_sum ^ b[row_index];
            }
            return x;
        }

        /**
         * Computes PLU factorization that breaks off as soon as y is in the image of the matrix.
         * Elimination is column wise and starts from column `start_col_idx`.
         * @param y
         * @param start_col_idx
         * @return true if y is in the image of the matrix, false otherwise.
         */
        bool rref_with_y_image_check(std::vector<uint8_t> &y, int start_col_idx = 0) {
            if (start_col_idx < 0) {
                throw std::invalid_argument("start_col_idx must be non-negative.");
            }
            if (start_col_idx == this->col_count) {
                auto in_image = true;
                for (int i = matrix_rank; i < this->row_count; i++) {
                    if (this->y_image_check_vector[i] == 1) {
                        in_image = false;
                        break;
                    }
                }
                return in_image;
            }
            if (start_col_idx == 0) {
                this->reset();
                this->y_image_check_vector = y;
                for (auto i = 0; i < this->row_count; i++) {
                    this->rows.push_back(i);
                }
            }

            auto previous_syndrome_size = this->y_image_check_vector.size();
            if (previous_syndrome_size < this->row_count) {
                this->y_image_check_vector.resize(this->row_count, 0);
                for (auto i = previous_syndrome_size; i < this->row_count; i++) {
                    this->y_image_check_vector[i] = y[i];
                }
            }

            auto y_sum = 0;
            for (auto i = 0; i < this->row_count; i++) {
                y_sum += y[i];
            }
            /*Trivial syndrome is always in image? What is the convention here?*/
            if (y_sum == 0) {
                this->LU_constructed = true;
                return true;
            }
            /*The vector we use to check whether y is in the image of the LU
            decomposition up to the current point of elimination.*/
            auto max_rank = std::min(this->row_count, this->col_count);
            /*Check whether the L matrix has the correct number of rows*/
            if (this->L.size() != this->row_count) {
                this->L.resize(this->row_count, std::vector<int>{});
            }

            bool in_image = false;
            //iterate over the columnsm starting from column `start_col_idx`
            for (auto col_idx = start_col_idx; col_idx < this->col_count; col_idx++) {
                //eliminate the column

                //exit if the maximum rank has been reached
                if (this->matrix_rank == max_rank) {
                    in_image = true;
                    break;
                }

                bool pivot = this->eliminate_column(col_idx, true, true);

                //check if y is in the image of the matrix
                if (pivot) {
                    std::swap(this->y_image_check_vector[this->matrix_rank - 1],
                              this->y_image_check_vector[this->swap_rows[this->matrix_rank - 1]]);
                    for (auto row_index: this->elimination_rows[this->matrix_rank - 1]) {
                        this->y_image_check_vector[row_index] ^= this->y_image_check_vector[this->matrix_rank - 1];
                    }
                    in_image = true;
                    for (int i = matrix_rank; i < this->row_count; i++) {
                        if (this->y_image_check_vector[i] == 1) {
                            in_image = false;
                            break;
                        }
                    }
                }
                //if y is in the image, we can stop eliminating
                if (in_image) {
                    break;
                }
            }
            this->LU_constructed = true;
            return in_image;
        }

        std::vector<uint8_t> fast_lu_solve(std::vector<uint8_t> &y) {
            bool y_in_image = this->rref_with_y_image_check(y);
            // this->LU_constructed = true;
            if (y_in_image) {
                return this->lu_solve(y);
            }
            throw std::invalid_argument("y is not in the image of the matrix.");

        }

        /**
         * Adds new, not eliminated column and corresponding rows to the matrix.
         * Updates row_count and col_count
         * @param new_matrix
         */
        void add_column_to_matrix(const std::vector<int> &new_col) {
            if (new_col.empty()) {
                return;
            }
            this->csc_mat.push_back(new_col);
            auto max_elem = std::max_element(new_col.begin(), new_col.end());
            // max_elem is highest row index in new_col so we need to add rows up to max_elem
            if (*max_elem >= this->row_count) {
                for (auto i = this->row_count; i <= *max_elem; i++) {
                    this->rows.push_back(i);
                }
                this->row_count = (*max_elem + 1);
            }
            this->col_count++;
        }
    };


    inline int rank(int row_count, int col_count, CscMatrix &csc_mat) {
        auto plu = PluDecomposition(row_count, col_count, csc_mat);
        plu.rref(false, false);
        return plu.matrix_rank;
    }

    inline CscMatrix kernel(int row_count, int col_count, CsrMatrix &csr_mat) {

        // To compute the kernel, we need to do PLU decomposition on the transpose of the matrix.
        // The CSR representation of mat is the CSC representation of mat.transpose().

        auto plu = PluDecomposition(col_count, row_count, csr_mat);
        plu.rref(false, false);

        std::vector<size_t> rr_col(col_count, 0);
        std::vector<std::vector<int>> ker;

        for (int j = 0; j < col_count; j++) {

            ker.emplace_back();

            std::fill(rr_col.begin(), rr_col.end(), 0);

            rr_col[j] = 1;

            for (int i = 0; i < plu.matrix_rank; i++) {

                std::swap(rr_col[i], rr_col[plu.swap_rows[i]]);
                if (rr_col[i] == 1) {
                    for (int add_row: plu.elimination_rows[i]) {
                        rr_col[add_row] ^= 1;
                    }
                }
            }

            for (int i = plu.matrix_rank; i < col_count; i++) {
                if (rr_col[i] == 1) {
                    ker[j].push_back(i - plu.matrix_rank);
                }
            }
        }
        return ker; // returns the kernel as a csc matrix
    }

    inline std::vector<int> pivot_rows(int row_count, int col_count, CsrMatrix &csr_mat) {
        auto plu = PluDecomposition(col_count, row_count, csr_mat);
        plu.rref(false, false);
        return plu.pivot_cols;
    }

    inline CsrMatrix row_span(int row_count, int col_count, CsrMatrix &csr_mat) {
        int row_permutations = std::pow(2, row_count);
        CsrMatrix row_span;

        for (int i = 0; i < row_permutations; i++) {

            std::vector<size_t> current_row(col_count, 0);
            auto row_add_indices = ldpc::util::decimal_to_binary_sparse(i, row_count);
            for (auto row_index: row_add_indices) {
                for (auto col_index: csr_mat[row_index]) {
                    current_row[col_index] ^= 1;
                }
            }
            std::vector<int> current_row_sparse;
            for (int j = 0; j < col_count; j++) {
                if (current_row[j] == 1) {
                    current_row_sparse.push_back(j);
                }
            }
            row_span.push_back(current_row_sparse);
        }
        return row_span;
    }

    struct DistanceStruct {
        int min_distance = INT_MAX;
        int samples_searched = 0;
        std::vector<std::vector<int>> min_weight_words;
    };

    inline DistanceStruct estimate_minimum_linear_row_combination(int row_count, int col_count, CsrMatrix &csr_mat,
                                                                  double timeout_seconds = 0,
                                                                  int number_of_words_to_save = 100) {

        DistanceStruct distance_struct;
        distance_struct.min_weight_words.resize(number_of_words_to_save, std::vector<int>{});
        int max_weight_saved_word = INT_MAX;
        int cc = 0;

        for (const auto &word: csr_mat) {
            int word_size = word.size();
            if (word_size < max_weight_saved_word) {
                int max1 = -10;
                int max2 = -10;
                int replace_word_index = 0;
                int count_index = 0;

                for (const auto &saved_word: distance_struct.min_weight_words) {
                    int saved_word_size = saved_word.size();
                    if (saved_word_size == 0) {
                        replace_word_index = count_index;
                        break;
                    }
                    if (saved_word_size > max1) {
                        max1 = saved_word_size;
                        replace_word_index = count_index;
                    } else if (saved_word_size > max2) {
                        max2 = saved_word_size;
                    }
                    count_index++;
                }

                distance_struct.min_weight_words[replace_word_index] = word;
                if (word_size > max2) {
                    max_weight_saved_word = word_size;
                } else {
                    max_weight_saved_word = max2;
                }

            }

        }


        double sample_prob = 2.0 / static_cast<double>(row_count);


        int row_permutations = std::pow(2, row_count);
        auto rand_gen = ldpc::rng::RandomNumberGenerator();

        auto start = std::chrono::high_resolution_clock::now();

        int count = 0;
        while (true) {
            count++;

            // auto rand = rand_gen.random_int(row_permutations-1);
            // auto row_add_indices = ldpc::util::decimal_to_binary_sparse(rand,row_count);

            std::vector<int> row_add_indices;
            for (int i = 0; i < row_count; i++) {
                if (rand_gen.random_double() < sample_prob) {
                    row_add_indices.push_back(i);
                }
            }

            std::vector<size_t> current_row(col_count, 0);
            for (auto row_index: row_add_indices) {
                for (auto col_index: csr_mat[row_index]) {
                    current_row[col_index] ^= 1;
                }
            }


            std::vector<int> current_row_sparse;
            for (int j = 0; j < col_count; j++) {
                if (current_row[j] == 1) {
                    current_row_sparse.push_back(j);
                }
            }

            int current_row_size = current_row_sparse.size();

            if (current_row_size == 0) {
                continue;
            }

            if (current_row_size < distance_struct.min_distance) {
                distance_struct.min_distance = current_row_size;
            }

            if (current_row_size <= max_weight_saved_word) {

                int max1 = -10;
                int max2 = -10;
                int replace_word_index = 0;
                int count_index = 0;
                for (const auto &word: distance_struct.min_weight_words) {

                    int word_size = word.size();

                    if (word_size == 0) {
                        replace_word_index = count_index;
                        break;
                    }
                    if (word_size > max1) {
                        max1 = word_size;
                        replace_word_index = count_index;
                    } else if (word_size > max2) {
                        max2 = word_size;
                    }
                    count_index++;
                }

                distance_struct.min_weight_words[replace_word_index] = current_row_sparse;
                if (current_row_size > max2) {
                    max_weight_saved_word = current_row_size;
                } else {
                    max_weight_saved_word = max2;
                }
            }

            auto now = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count() / 1000.0;
            if (elapsed >= timeout_seconds) {
                break;
            }

        }

        distance_struct.samples_searched = count;

        return distance_struct;
    }

    inline DistanceStruct
    estimate_code_distance(int row_count, int col_count, CsrMatrix &csr_mat, double timeout_seconds = 0,
                           int number_of_words_to_save = 100) {

        auto ker = ldpc::gf2dense::kernel(row_count, col_count, csr_mat);
        // convert to csr_matrix
        int max_row = -1;
        for (auto &col_index: ker) {
            for (int row_index: col_index) {
                if (row_index > max_row) {
                    max_row = row_index;
                }
            }
        }

        CsrMatrix ker_csr = CscMatrix(max_row + 1, std::vector<int>{});

        for (int col_index = 0; col_index < ker.size(); col_index++) {
            for (int row_index: ker[col_index]) {
                ker_csr[row_index].push_back(col_index);
            }
        }

        return ldpc::gf2dense::estimate_minimum_linear_row_combination(max_row + 1, col_count, ker_csr,
                                                                       timeout_seconds,
                                                                       number_of_words_to_save);

    }

    inline int compute_exact_code_distance(int row_count, int col_count, CsrMatrix &csr_mat) {

        CscMatrix ker = ldpc::gf2dense::kernel(row_count, col_count, csr_mat);


        CsrMatrix ker_csr = csc_to_csr(ker);

        int row_permutations = std::pow(2, ker_csr.size());

        int distance = col_count;

        for (int i = 1; i < row_permutations; i++) {

            std::vector<size_t> current_row(col_count, 0);

            auto row_add_indices = ldpc::util::decimal_to_binary_sparse(i, ker_csr.size());

            int row_count = 0;
            for (auto row_index: row_add_indices) {
                for (auto col_index: ker_csr[row_index]) {
                    if (current_row[col_index] == 0) {
                        row_count++;
                    } else {
                        row_count--;
                    }
                }
            }

            if (row_count < distance) {
                distance = row_count;
            }

        }

        return distance;

    }

    inline int count_non_zero_matrix_entries(const CscMatrix &csc_mat) {
        auto count = 0;
        for (const auto &col: csc_mat) {
            count += col.size();
        }
        return count;
    }
}   // namespace ldpc::gf2dense
//end namespace ldpc

#endif