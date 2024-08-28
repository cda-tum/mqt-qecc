/**
 * @file gf2sparse.h
 * @brief A C++ library for sparse matrices in GF(2)
 */

#ifndef GF2SPARSE_H
#define GF2SPARSE_H

#include <vector>
#include <iterator>
#include "ldpc/sparse_matrix_base.hpp"


    namespace ldpc::gf2sparse {

/**
 * @brief An entry in a sparse matrix over GF(2)
 */
        class GF2Entry : public ldpc::sparse_matrix_base::EntryBase<GF2Entry> {
        public:
            ~GF2Entry() = default;
        };

/**
 * @brief A sparse matrix over GF(2)
 */
        template<class ENTRY_OBJ = GF2Entry>
        class GF2Sparse : public ldpc::sparse_matrix_base::SparseMatrixBase<ENTRY_OBJ> {
        public:
            typedef ldpc::sparse_matrix_base::SparseMatrixBase<ENTRY_OBJ> BASE;

            /**
             * @brief Constructor for creating a new GF2Sparse object with the given dimensions
             * @param rows The number of rows in the matrix
             * @param cols The number of columns in the matrix
             */
            GF2Sparse(int rows, int cols, int entry_count = 0) : ldpc::sparse_matrix_base::SparseMatrixBase<ENTRY_OBJ>(
                    rows, cols, entry_count) {}


            /**
             * @brief Default constructor for creating a new GF2Sparse object.
             *        Creates an empty GF2Sparse matrix with zero rows and columns.
             */
            GF2Sparse() = default;


            /**
             * @brief Destructor for GF2Sparse object
             */
            ~GF2Sparse() = default;

            /**
             * @brief Inserts a row of entries in compressed sparse row (CSR) format
             * @param row_index The index of the row to insert
             * @param column_indices The indices of the non-zero entries in the row
             */
            void csr_row_insert(int row_index, std::vector<int> &column_indices);

            /**
             * @brief Inserts a matrix in CSR format
             * @param csr_matrix The matrix to insert
             */
            void csr_insert(std::vector<std::vector<int>> &csr_matrix);

            /**
             * @brief Inserts a matrix in CSR format
             * @param csr_matrix The matrix to insert
             */
            void csc_insert(std::vector<std::vector<int>> &csc_matrix);

            /**
             * @brief Multiplies the matrix by a vector and stores the result in another vector
             * @param input_vector The vector to multiply the matrix with
             * @param output_vector The vector to store the result in
             * @return A reference to the output vector
             */
            std::vector<uint8_t> &mulvec(std::vector<uint8_t> &input_vector, std::vector<uint8_t> &output_vector);

            /**
             * @brief Multiplies the matrix by a vector and returns the result as a new vector
             * @param input_vector The vector to multiply the matrix with
             * @return The resulting vector
             */
            std::vector<uint8_t> mulvec(std::vector<uint8_t> &input_vector);

            // /**
            //  * @brief Multiplies the matrix by a vector with parallel execution using OpenMP and stores the result in another vector
            //  * @param input_vector The vector to multiply the matrix with
            //  * @param output_vector The vector to store the result in
            //  * @return A reference to the output vector
            //  */
            // std::vector<uint8_t>& mulvec_parallel(std::vector<uint8_t>& input_vector, std::vector<uint8_t>& output_vector);

            /**
             * @brief Multiplies the matrix by another matrix and returns the result as a new matrix
             * @tparam ENTRY_OBJ2 The type of entries in the matrix to be multiplied with
             * @param mat_right The matrix to multiply with
             * @return The resulting matrix
             */
            template<typename ENTRY_OBJ2>
            GF2Sparse<ENTRY_OBJ> matmul(GF2Sparse<ENTRY_OBJ2> &mat_right);

            /**
             * @brief Adds two rows together
             * @param i The row to overwrite
             * @param j The row to add to row i
             */
            void add_rows(int i, int j);

            /**
             * @brief Transposes the matrix
             * @return The transposed matrix
             */
            GF2Sparse<ENTRY_OBJ> transpose();

            /**
             * @brief Compares two matrices for equality
             * @tparam ENTRY_OBJ2 The type of entries in the matrix to compare with
             * @param matrix2 The matrix to compare with
             * @return True if the matrices are equal, false otherwise
             */
            template<typename ENTRY_OBJ2>
            bool gf2_equal(GF2Sparse<ENTRY_OBJ2> &matrix2);

            /**
             * @brief Allocates memory for the sparse matrix
             *
             * @param m The row count
             * @param n The column count
             * @param entry_count The number of entries in the matrix (optional)
             * @return ** void
             */
            void allocate(int m, int n, int entry_count = 0);

        };

        template<class ENTRY_OBJ>
        void GF2Sparse<ENTRY_OBJ>::allocate(int m, int n, int entry_count) {
            this->initialise_sparse_matrix(m, n, entry_count);
        }

        template<class ENTRY_OBJ>
        void GF2Sparse<ENTRY_OBJ>::csr_row_insert(int row_index, std::vector<int> &column_indices) {
            // Iterate through each column index in the vector
            for (auto col_index: column_indices) {
                // Insert a new entry with the given row and column indices
                this->insert_entry(row_index, col_index);
            }
        }

        template<class ENTRY_OBJ>
        void GF2Sparse<ENTRY_OBJ>::csr_insert(std::vector<std::vector<int>> &csr_matrix) {
            int i = 0;
            // Iterate through each row of the matrix
            for (auto row: csr_matrix) {
                // Insert the row of entries in compressed sparse row (CSR) format
                this->csr_row_insert(i, row);
                i++;
            }
        }

        template<class ENTRY_OBJ>
        void GF2Sparse<ENTRY_OBJ>::csc_insert(std::vector<std::vector<int>> &csc_matrix) {
            int i = 0;
            // Iterate through each column of the matrix
            for (const auto &col: csc_matrix) {
                // Insert the column of entries in compressed sparse column (CSC) format
                for (auto row_index: col) {
                    this->insert_entry(row_index, i);
                }
                i++;
            }
        }


        template<class ENTRY_OBJ>
        std::vector<uint8_t> &
        GF2Sparse<ENTRY_OBJ>::mulvec(std::vector<uint8_t> &input_vector, std::vector<uint8_t> &output_vector) {
            // Initialize the output vector to all zeros
            for (int i = 0; i < this->m; i++) {
                output_vector[i] = 0;
            }

            // Iterate through each row of the matrix
            for (int i = 0; i < this->m; i++) {
                // Iterate through each non-zero entry in the row
                for (auto &e: this->iterate_row(i)) {
                    // Compute the XOR of the current output value with the value in the input vector at the entry's column index
                    output_vector[i] ^= input_vector[e.col_index];
                }
            }

            // Return the output vector
            return output_vector;
        }

        template<class ENTRY_OBJ>
        std::vector<uint8_t> GF2Sparse<ENTRY_OBJ>::mulvec(std::vector<uint8_t> &input_vector) {
            // Initialize the output vector to all zeros
            std::vector<uint8_t> output_vector(this->m, 0);

            // Iterate through each row of the matrix
            for (int i = 0; i < this->m; i++) {
                // Iterate through each non-zero entry in the row
                for (auto &e: this->iterate_row(i)) {
                    // Compute the XOR of the current output value with the value in the input vector at the entry's column index
                    output_vector[i] ^= input_vector[e.col_index];
                }
            }

            // Return the output vector
            return output_vector;
        }


// template<class ENTRY_OBJ>
// std::vector<uint8_t>& GF2Sparse<ENTRY_OBJ>::mulvec_parallel(std::vector<uint8_t>& input_vector, std::vector<uint8_t>& output_vector){
//     // Initialize the output vector to all zeros
//     #pragma omp for
//     for(int i = 0; i<this->m; i++) output_vector[i] = 0;

//     // Iterate through each row of the matrix
//     #pragma omp for
//     for(int i = 0; i < this->m; i++){
//         // Iterate through each non-zero entry in the row
//         for(auto& e: this->iterate_row(i)){
//             // Compute the XOR of the current output value with the value in the input vector at the entry's column index
//             output_vector[i] ^= input_vector[e.col_index];
//         }
//     }

//     // Return the output vector
//     return output_vector;
// }


        template<typename ENTRY_OBJ>
        template<typename ENTRY_OBJ2>
        GF2Sparse<ENTRY_OBJ> GF2Sparse<ENTRY_OBJ>::matmul(GF2Sparse<ENTRY_OBJ2> &mat_right) {

            // Check if the dimensions of the input matrices are valid for multiplication
            if (this->n != mat_right.m) {
                throw std::invalid_argument("Input matrices have invalid dimensions!");
            }

            // Create a new GF2Sparse matrix to store the output
            auto output_mat = GF2Sparse<ENTRY_OBJ>(this->m, mat_right.n);

            // Iterate over each row and column of the output matrix
            for (int i = 0; i < output_mat.m; i++) {
                for (int j = 0; j < output_mat.n; j++) {
                    int sum = 0;
                    // Iterate over the non-zero entries in the column of the right-hand matrix
                    for (auto &e: mat_right.iterate_column(j)) {
                        // Iterate over the non-zero entries in the row of this matrix
                        for (auto &g: this->iterate_row(i)) {
                            // Check if the column index of this matrix matches the row index of the right-hand matrix
                            if (g.col_index == e.row_index) {
                                sum ^= 1;
                            }
                        }
                    }
                    // Insert an entry in the output matrix if the sum is non-zero
                    if (sum) {
                        output_mat.insert_entry(i, j);
                    }
                }
            }

            // Return a shared pointer to the output matrix
            return output_mat;
        }


        template<class ENTRY_OBJ>
        void GF2Sparse<ENTRY_OBJ>::add_rows(int i, int j) {

            // Set row i as the row that will be overwritten
            // Initialize variables
            bool intersection = 0;
            std::vector<ENTRY_OBJ *> entries_to_remove;

            for (auto &g: this->iterate_row(j)) {
                intersection = false;
                int col_index = g.col_index;
                for (auto &e: this->iterate_column(col_index)) {
                    if (e.row_index == i) {
                        entries_to_remove.push_back(&e);
                        intersection = true;
                    }
                }
                // If there was no intersection between the entries, insert the entry from row j into row i
                if (!intersection) {
                    this->insert_entry(i, col_index);
                }

            }

            // Remove all the entries from row i that were marked for removal
            for (auto e: entries_to_remove) {
                this->remove(*e);
            }

        }

        template<class ENTRY_OBJ>
        GF2Sparse<ENTRY_OBJ> GF2Sparse<ENTRY_OBJ>::transpose() {

            // Create a new GF2Sparse matrix with the dimensions of the transposed matrix
            auto pcmT = GF2Sparse<ENTRY_OBJ>(this->n, this->m);

            // Iterate over each row of this matrix
            for (int i = 0; i < this->m; i++) {
                // Iterate over each non-zero entry in the row
                for (auto &e: this->iterate_row(i)) {
                    // Insert the entry into the transposed matrix with the row and column indices swapped
                    pcmT.insert_entry(e.col_index, e.row_index);
                }
            }

            // Return a shared pointer to the transposed matrix
            return pcmT;
        }

        template<class ENTRY_OBJ>
        template<typename ENTRY_OBJ2>
        bool GF2Sparse<ENTRY_OBJ>::gf2_equal(GF2Sparse<ENTRY_OBJ2> &matrix2) {

            // Check if the dimensions of the matrices are equal
            if (this->n != matrix2.n || this->m != matrix2.m) {
                return false;
            }

            // Iterate over each row of this matrix
            for (int i = 0; i < this->m; i++) {

                // Iterate over each non-zero entry in the row of this matrix
                for (auto &e: this->iterate_row(i)) {

                    // Get the corresponding entry in the same position of the other matrix
                    auto &g = matrix2.get_entry(e.row_index, e.col_index);

                    // If there is no corresponding entry in the same position of the other matrix, the matrices are not equal
                    if (g.at_end()) {
                        return false;
                    }
                }

                // Iterate over each non-zero entry in the row of the other matrix
                for (auto &e: matrix2.iterate_row(i)) {

                    // Get the corresponding entry in the same position of this matrix
                    auto &g = this->get_entry(e.row_index, e.col_index);

                    // If there is no corresponding entry in the same position of this matrix, the matrices are not equal
                    if (g.at_end()) {
                        return false;
                    }
                }
            }

            // If all non-zero entries in both matrices are the same, the matrices are equal
            return true;
        }

/**
 * @brief Compares two GF2Sparse matrices for equality
 * @tparam ENTRY_OBJ1 The type of entries in the first matrix
 * @tparam ENTRY_OBJ2 The type of entries in the second matrix
 * @param matrix1 A shared pointer to the first matrix
 * @param matrix2 A shared pointer to the second matrix
 * @return True if the matrices are equal, false otherwise
 */
        template<typename ENTRY_OBJ1, typename ENTRY_OBJ2>
        bool operator==(GF2Sparse<ENTRY_OBJ1> &matrix1, GF2Sparse<ENTRY_OBJ2> &matrix2) {
            return matrix1.gf2_equal(matrix2);
        }

/**
 * @brief Creates a new GF2Sparse identity matrix with the given dimensions
 * @tparam ENTRY_OBJ The type of entries in the matrix
 * @param n The number of rows/columns in the matrix
 * @return The GF2Sparse identity matrix
 */
        template<class ENTRY_OBJ = GF2Entry>
        GF2Sparse<ENTRY_OBJ> identity(int n) {
            // Create a new GF2Sparse matrix with the given dimensions
            auto matrix = GF2Sparse<ENTRY_OBJ>(n, n);

            // Insert non-zero entries along the diagonal
            for (int i = 0; i < n; i++) {
                matrix.insert_entry(i, i);
            }

            // Return a shared pointer to the new GF2Sparse identity matrix
            return matrix;
        }


        template<class GF2SPARSE_MATRIX_CLASS>
        GF2SPARSE_MATRIX_CLASS copy_cols(GF2SPARSE_MATRIX_CLASS &mat, const std::vector<int> &cols) {
            int m = 0;
            int n = 0;
            int i = 0;
            int j = 0;
            m = mat.m;
            n = cols.size();
            auto copy_mat = GF2SPARSE_MATRIX_CLASS(m, n);
            int new_col_index = -1;
            for (auto col_index: cols) {
                new_col_index += 1;
                for (auto &e: mat.iterate_column(col_index)) {
                    copy_mat.insert_entry(e.row_index, new_col_index);
                }
            }
            return copy_mat;
        }


        template<class GF2MATRIX>
        GF2MATRIX vstack(std::vector<GF2MATRIX> &mats) {

            int mat_count = mats.size();
            int n = mats[0].n;
            int m0 = mats[0].m;

            int m = 0;
            for (auto mat: mats) {
                if (mat.n != n) {
                    throw std::invalid_argument("All matrices must have the same number of columns");
                }
                m += mat.m;
            }

            auto stacked_mat = GF2MATRIX(m, n);

            int row_offset = 0;
            for (auto mat: mats) {
                for (auto i = 0; i < mat.m; i++) {
                    for (auto &e: mat.iterate_row(i)) {
                        stacked_mat.insert_entry(row_offset + e.row_index, e.col_index);
                    }
                }
                row_offset += mat.m;
            }

            return stacked_mat;

        }

        template<class ENTRY_OBJ>
        GF2Sparse<ENTRY_OBJ> hstack(std::vector<GF2Sparse<ENTRY_OBJ>> &mats) {

            int mat_count = mats.size();
            int m = mats[0].m;

            int n = 0;

            for (auto mat: mats) {
                if (mat.m != m) {
                    throw std::invalid_argument("All matrices must have the same number of rows");
                }
                n += mat.n;
            }

            auto stacked_mat = GF2Sparse<ENTRY_OBJ>(m, n);

            int col_offset = 0;
            for (auto mat: mats) {
                for (auto i = 0; i < mat.m; i++) {
                    for (auto &e: mat.iterate_row(i)) {
                        stacked_mat.insert_entry(e.row_index, col_offset + e.col_index);
                    }
                }
                col_offset += mat.n;
            }

            return stacked_mat;


        }

        template<class ENTRY_OBJ>
        GF2Sparse<ENTRY_OBJ> hstack(GF2Sparse<ENTRY_OBJ> &mat1, GF2Sparse<ENTRY_OBJ> &mat2) {

            int m = 0;
            int n = 0;

            n = mat1.n + mat2.n;
            m = mat1.m;
            if (mat1.m != mat2.m) {
                throw std::invalid_argument("All matrices must have the same number of rows");
            }

            auto stacked_mat = GF2Sparse<ENTRY_OBJ>(m, n);

            for (int i = 0; i < m; i++) {
                for (auto &e: mat1.iterate_row(i)) {
                    stacked_mat.insert_entry(e.row_index, e.col_index);
                }

                for (auto &e: mat2.iterate_row(i)) {
                    stacked_mat.insert_entry(e.row_index, e.col_index + mat1.n);
                }
            }

            return stacked_mat;


        }

        template<class GF2MATRIX>
        GF2MATRIX kron(GF2MATRIX &mat1, GF2MATRIX &mat2) {

            int m1 = 0;
            int n1 = 0;
            int m2 = 0;
            int n2 = 0;
            m1 = mat1.m;
            n1 = mat1.n;
            m2 = mat2.m;
            n2 = mat2.n;

            auto kron_mat = GF2MATRIX(m1 * m2, n1 * n2);

            for (auto i = 0; i < m1; i++) {
                for (auto &e: mat1.iterate_row(i)) {
                    int row_offset = e.row_index * m2;
                    int col_offset = e.col_index * n2;

                    for (auto j = 0; j < m2; j++) {
                        for (auto &f: mat2.iterate_row(j)) {
                            kron_mat.insert_entry(row_offset + f.row_index, col_offset + f.col_index);
                        }
                    }

                }

            }

            return kron_mat;

        }

        template<class GF2MATRIX>
        GF2MATRIX operator+(GF2MATRIX &mat1, GF2MATRIX &mat2) {

            if (mat1.m != mat2.m || mat1.n != mat2.n) {
                throw std::invalid_argument("Matrices must have the same dimensions");
            }

            //Deepcopy mat2
            auto sum_mat = GF2Sparse<GF2MATRIX>(mat2.m, mat2.n);
            for (int i = 0; i < mat2.m; i++) {
                for (auto &e: mat2.iterate_row(i)) {
                    sum_mat.insert_entry(e.row_index, e.col_index);
                }
            }

            // Add mat1 to mat2
            for (int i = 0; i < mat1.m; i++) {
                for (auto &e: mat1.iterate_row(i)) {
                    auto &g = sum_mat.get_entry(e.row_index, e.col_index);
                    if (!g.at_end()) {
                        sum_mat.remove(g);
                    } else {
                        sum_mat.insert_entry(e.row_index, e.col_index);
                    }
                }
            }

            return sum_mat;
        }


        template<class ENTRY_OBJ = GF2Entry>
        GF2Sparse<ENTRY_OBJ> csr_to_gf2sparse(std::vector<std::vector<int>> &csr_matrix) {
            int row_count = csr_matrix.size();
            int col_count = 0;
            for (const auto &row: csr_matrix) {
                for (auto col: row) {
                    if (col > col_count) {
                        col_count = col;
                    }
                }
            }

            auto gf2sparse_mat = GF2Sparse<ENTRY_OBJ>(row_count, col_count + 1);

            for (auto row_index = 0; row_index < row_count; row_index++) {
                for (auto col_index: csr_matrix[row_index]) {
                    gf2sparse_mat.insert_entry(row_index, col_index);
                }
            }

            return gf2sparse_mat;

        }

        template<class ENTRY_OBJ = GF2Entry>
        GF2Sparse<ENTRY_OBJ> csc_to_gf2sparse(std::vector<std::vector<int>> &csc_matrix) {
            int col_count = csc_matrix.size();
            int row_count = 0;
            for (const auto &col: csc_matrix) {
                for (auto row: col) {
                    if (row > row_count) {
                        row_count = row;
                    }
                }
            }

            auto gf2sparse_mat = GF2Sparse<ENTRY_OBJ>(row_count + 1, col_count);

            for (auto col_index = 0; col_index < col_count; col_index++) {
                for (auto row_index: csc_matrix[col_index]) {
                    gf2sparse_mat.insert_entry(row_index, col_index);
                }
            }

            return gf2sparse_mat;

        }


    } // end namespace gf2sparse
// end namespace ldpc


using cygf2_entry = ldpc::gf2sparse::GF2Entry;
using cygf2_sparse = ldpc::gf2sparse::GF2Sparse<ldpc::gf2sparse::GF2Entry>;

#endif