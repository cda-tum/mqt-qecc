#ifndef SPARSE_MATRIX_BASE_H
#define SPARSE_MATRIX_BASE_H

#include <vector>
#include <iterator>

namespace ldpc {

    namespace sparse_matrix_base {

/**
 * @brief Base class for defining the node types for Sparse Matrices
 * 
 * This class defines the basic properties of a node in a sparse matrix such as its row index, column index,
 * and pointers to its neighboring nodes in the same row and column. Each node class that derives from this
 * base class will inherit these properties and add any additional properties as required by the specific
 * sparse matrix implementation.
 * 
 * @tparam DERIVED The derived class from EntryBase.
 */
        template<class DERIVED>
        class EntryBase {
        public:
            int row_index = -100; /**< The row index of the matrix element */
            int col_index = -100; /**< The column index of the matrix element */
            DERIVED *left = static_cast<DERIVED *>(this); /**< Pointer to the previous element in the row */
            DERIVED *right = static_cast<DERIVED *>(this); /**< Pointer to the next element in the row */
            DERIVED *up = static_cast<DERIVED *>(this); /**< Pointer to the previous element in the column */
            DERIVED *down = static_cast<DERIVED *>(this); /**< Pointer to the next element in the column */

            /**
             * @brief Resets the values of the entry to their default values.
             *
             * This function resets the values of the entry to their default values. This is useful for when an
             * entry is removed from the matrix and needs to be returned to its default state for re-use.
             */
            void reset() {
                row_index = -100;
                col_index = -100;
                left = static_cast<DERIVED *>(this);
                right = static_cast<DERIVED *>(this);
                up = static_cast<DERIVED *>(this);
                down = static_cast<DERIVED *>(this);
            }

            /**
             * @brief Checks if the entry is at the end of the list
             *
             * This function checks if the entry is at the end of the list by checking if its row index is equal
             * to -100. If it is, then the function returns true to indicate that the entry is at the end of the
             * list.
             *
             * @return True if the entry is at the end, false otherwise
             */
            bool at_end() {
                if (row_index == -100) {
                    {
                        return true;
                    }
                } else {
                    {
                        return false;
                    }
                }
            }

            /**
             * @brief Returns a string representation of the entry
             *
             * This function returns a string representation of the entry. In this implementation, the function
             * always returns "1", but in other implementations, this function could be used to return the value
             * stored in the entry or some other useful information.
             *
             * @return The string representation of the entry
             */
            std::string str() {
                return "1";
            }

        };


        struct CsrMatrix {
            int m;
            int n;
            int entry_count;
            std::vector<std::vector<int>> row_adjacency_list;
        };


/**
 * @brief Template base class for implementing sparse matrices in a doubly linked list format
 * 
 * @tparam ENTRY_OBJ The entry object class that the sparse matrix will use
 * for its entries. This class should contain fields for row index, column index,
 * and value.
 * 
 * This class allows for the construction of sparse matrices with custom data types by
 * passing node objects via the ENTRY_OBJ template parameter. The matrix is stored as a
 * doubly linked list, where each row and column is represented by a linked list of entries.
 * Each entry contains a reference to the next and previous entries in its row and column,
 * respectively. This format allows for efficient insertion and deletion of entries in the matrix,
 * especially when the matrix is large and sparse.
 */
        template<class ENTRY_OBJ>
        class SparseMatrixBase {
        public:
            int m, n; //m: check-count; n: bit-count // todo refactor to: check_cnt, bit_cnt
            int node_count;
            int entry_block_size;
            int allocated_entry_count;
            int released_entry_count;
            int block_position;
            int block_idx;
            std::vector<std::vector<ENTRY_OBJ>> entries;
            std::vector<ENTRY_OBJ *> removed_entries;
            std::vector<ENTRY_OBJ *> row_heads; //starting point for each row
            std::vector<ENTRY_OBJ *> column_heads; //starting point for each column
            bool MEMORY_ALLOCATED = false;

            // std::vector<ENTRY_OBJ*> matrix_entries;

            SparseMatrixBase() = default;


            SparseMatrixBase(int check_count, int bit_count, int entry_count = 0) {
                this->initialise_sparse_matrix(check_count, bit_count, entry_count);
            }

            /**
             * @brief Constructs a sparse matrix with the given dimensions
             *
             * @param check_count The number of rows in the matrix
             * @param bit_count The number of columns in the matrix
             */
            void initialise_sparse_matrix(int check_count, int bit_count, int entry_count = 0) {

                this->reset_matrix();
                this->m = check_count;
                this->n = bit_count;
                this->block_idx = -1;
                this->released_entry_count = 0;
                this->allocated_entry_count = 0;
                this->entry_block_size = this->m + this->n + entry_count;
                this->allocate_memory();
                this->entry_block_size = this->m + this->n;
            }

            void reset_matrix() {

                if (this->MEMORY_ALLOCATED) {
                    this->column_heads.clear();
                    this->row_heads.clear();
                    this->removed_entries.clear();
                    for (auto entry_block: this->entries) {
                        entry_block.clear();
                    }
                    this->entries.clear();
                }
                this->m = 0;
                this->n = 0;
                this->block_idx = -1;
                this->released_entry_count = 0;
                this->allocated_entry_count = 0;
                this->entry_block_size = 0;
                this->entry_block_size = 0;
                this->MEMORY_ALLOCATED = false;
            }

            /**
             * @brief Destructor for SparseMatrixBase. Frees memory occupied by entries.
             */
            ~SparseMatrixBase() {
                this->reset_matrix();
            }

            /**
             * @brief Allocates a new entry object and returns a pointer to it.
             *
             * If there are any entries that have been removed from the matrix, the function returns the last
             * removed entry. Otherwise, if there is space in the entries vector, a new entry is allocated at
             * the end of the vector. If the entries vector is full, the function allocates a new block of
             * entries and returns the first entry in the new block.
             *
             * @return A pointer to a new entry object.
             */
            ENTRY_OBJ *allocate_new_entry() {
                // if there are any previously removed entries, use them instead of creating new ones
                if (this->removed_entries.size() != 0) {
                    auto e_ptr = this->removed_entries.back();
                    this->removed_entries.pop_back();
                    return e_ptr;
                }
                // if there are no previously removed entries, create a new one
                // if there is no space for the new entry, add a new block of entries
                if (this->released_entry_count == this->allocated_entry_count) {
                    this->entries.push_back(std::vector<ENTRY_OBJ>(this->entry_block_size));
                    this->allocated_entry_count += this->entry_block_size;
                    this->block_idx++;
                    this->block_position = 0;
                }



                // use the next available entry in the pool
                auto e_ptr = &this->entries[block_idx][this->block_position];
                this->block_position++;
                this->released_entry_count++;
                return e_ptr;
            }


            /**
             * @brief Returns the number of non-zero entries in the matrix.
             *
             * @return The number of non-zero entries in the matrix.
             */
            int entry_count() {
                return this->released_entry_count - this->n - this->m - this->removed_entries.size();
            }

            /**
             * @brief Computes the sparsity of the matrix
             *
             * The sparsity of a matrix is defined as the ratio of the number of zero elements in the
             * matrix to the total number of elements in the matrix. This function computes the
             * sparsity of the matrix represented by this object, and returns the result as a double value.
             *
             * @return The sparsity of the matrix as a double value.
             */
            double sparsity() {
                return double(this->entry_count()) / double(this->m * this->n);
            }

            /**
            * @brief Allocates memory for the row and column header nodes.
            *
            * This function resizes the row_heads and column_heads vectors to m and n respectively.
            * For each row and column header node, it allocates memory for a new entry object, sets
            * the row and column indices to -100 to indicate it is not a value element, and sets the right,
            * left, up, and down pointers to point to itself since there are no other nodes in the same row or column yet.
            */
            void allocate_memory() {

                this->MEMORY_ALLOCATED = true;

                this->row_heads.resize(this->m); // resize row_heads vector to m
                this->column_heads.resize(this->n); // resize column_heads vector to n

                // create and initialize each row header node
                for (int i = 0; i < this->m; i++) {
                    ENTRY_OBJ *row_entry_ptr = this->allocate_new_entry(); // allocate memory for a new entry object
                    row_entry_ptr->row_index = -100; // set row index to -100 to indicate it is not a value element
                    row_entry_ptr->col_index = -100; // set col index to -100 to indicate it is not a value element
                    row_entry_ptr->right = row_entry_ptr; // point to itself since there are no other nodes in the same row yet
                    row_entry_ptr->left = row_entry_ptr; // point to itself since there are no other nodes in the same row yet
                    row_entry_ptr->up = row_entry_ptr; // point to itself since there are no other nodes in the same column yet
                    row_entry_ptr->down = row_entry_ptr; // point to itself since there are no other nodes in the same column yet
                    this->row_heads[i] = row_entry_ptr; // add the new row header node to the row_heads vector
                }

                // create and initialize each column header node
                for (int i = 0; i < this->n; i++) {
                    ENTRY_OBJ *col_entry_ptr = this->allocate_new_entry(); // allocate memory for a new entry object
                    col_entry_ptr->row_index = -100; // set row index to -100 to indicate it is not a value element
                    col_entry_ptr->col_index = -100; // set col index to -100 to indicate it is not a value element
                    col_entry_ptr->right = col_entry_ptr; // point to itself since there are no other nodes in the same column yet
                    col_entry_ptr->left = col_entry_ptr; // point to itself since there are no other nodes in the same column yet
                    col_entry_ptr->up = col_entry_ptr; // point to itself since there are no other nodes in the same row yet
                    col_entry_ptr->down = col_entry_ptr; // point to itself since there are no other nodes in the same row yet
                    this->column_heads[i] = col_entry_ptr; // add the new column header node to the column_heads vector
                }
            }


            /**
             * @brief Swaps two rows of the matrix
             *
             * This function swaps rows i and j of the matrix. All row entries must be updated accordingly.
             *
             * @param i The index of the first row to swap
             * @param j The index of the second row to swap
             */
            void swap_rows(int i, int j) {
                auto tmp1_ptr = this->row_heads[i]; // store the head element of row i in a temporary variable
                auto tmp2_ptr = this->row_heads[j]; // store the head element of row j in a temporary variable
                this->row_heads[j] = tmp1_ptr; // set the head element of row j to the head element of row i
                this->row_heads[i] = tmp2_ptr; // set the head element of row i to the head element of row j
                for (auto &e: this->iterate_row(i)) {
                    {
                        e.row_index = i; // update the row index of all elements in row i to j
                    }
                }
                for (auto &e: this->iterate_row(j)) {
                    {
                        e.row_index = j; // update the row index of all elements in row j to i
                    }
                }
            }


            void reorder_rows(std::vector<int> rows) {

                std::vector<ENTRY_OBJ *> temp_row_heads = this->row_heads;
                for (int i = 0; i < m; i++) {
                    this->row_heads[i] = temp_row_heads[rows[i]];
                    for (auto &e: this->iterate_row(i)) {
                        e.row_index = i;
                    }
                }

            }

            /**
             * Swaps two columns in the sparse matrix.
             *
             * @param i The index of the first column to swap.
             * @param j The index of the second column to swap.
             */
            void swap_columns(int i, int j) {
                auto tmp1 = this->column_heads[i];
                auto tmp2 = this->column_heads[j];
                this->column_heads[j] = tmp1;
                this->column_heads[i] = tmp2;
                // update the column indices for all entries in columns i and j
                for (auto &e: this->iterate_column(i)) {
                    e.col_index = i;
                }
                for (auto &e: this->iterate_column(j)) {
                    e.col_index = j;
                }
            }

            /**
             * @brief Gets the number of non-zero entries in a row of the matrix.
             *
             * This function returns the degree of a given row in the matrix, i.e., the number of non-zero entries in the row.
             *
             * @param row The index of the row to get the degree of.
             * @return The number of non-zero entries in the row.
             */
            int get_row_degree(int row) {
                return abs(this->row_heads[row]->col_index) - 100;
            }

            /**
             * @brief Gets the number of non-zero entries in a column of the matrix.
             *
             * This function returns the degree of a given column in the matrix, i.e., the number of non-zero entries in the column.
             *
             * @param col The index of the column to get the degree of.
             * @return The number of non-zero entries in the column.
             */
            int get_col_degree(int col) {
                return abs(this->column_heads[col]->col_index) - 100;
            }

            /**
             * @brief Removes an entry from the matrix.
             *
             * This function removes a given entry from the matrix. The entry is identified by its row and column indices.
             *
             * @param i The row index of the entry to remove.
             * @param j The column index of the entry to remove.
             */
            void remove_entry(int i, int j) {
                auto &e = this->get_entry(i, j);
                this->remove(e);
            }

            /**
             * @brief Removes an entry from the matrix and updates the row/column weights
             *
             * @param e Pointer to the entry object to be removed
             */
            void remove(ENTRY_OBJ &e_ref) {
                ENTRY_OBJ *e_ptr = &e_ref;
                // Check if the entry is not already at the end of the row or column.
                if (!e_ptr->at_end()) {
                    // Store pointers to the adjacent entries.
                    auto e_left_ptr = e_ptr->left;
                    auto e_right_ptr = e_ptr->right;
                    auto e_up_ptr = e_ptr->up;
                    auto e_down_ptr = e_ptr->down;

                    // Update pointers of the adjacent entries to remove the entry from the linked list.
                    e_left_ptr->right = e_right_ptr;
                    e_right_ptr->left = e_left_ptr;
                    e_up_ptr->down = e_down_ptr;
                    e_down_ptr->up = e_up_ptr;

                    /* Update the row/column weights. Note that this information is stored in the
                    ENTRY_OBJ.col_index field as a negative number (to avoid confusion with
                    an actual column index). To get the row/column weights, use the get_row_weight()
                    and get_col_weight() functions. */
                    this->row_heads[e_ptr->row_index]->col_index++;
                    this->column_heads[e_ptr->col_index]->col_index++;

                    // Reset the entry to default values before pushing it to the buffer.
                    e_ptr->reset();
                    // Store the removed entry in the buffer for later reuse.
                    this->removed_entries.push_back(e_ptr);
                }
            }

            /**
            * @brief Inserts a new entry in the matrix at position (j, i).
            *
            * @param j The row index of the new entry.
            * @param i The column index of the new entry.
            * @return A reference to the newly created ENTRY_OBJ object.
            * @throws std::std::invalid_argument if either i or j is out of bounds.
            *
            * This function inserts a new entry in the matrix at position (j, i). If an entry
            * already exists at that position, this function simply returns a reference to it.
            * Otherwise, it creates a new entry and inserts it into the matrix, linking it to
            * the neighboring entries to maintain the doubly linked structure. This function
            * also updates the row and column weights of the matrix. The row and column weights
            * are stored as negative integers in the col_index field of the ENTRY_OBJ object.
            * To retrieve the actual weights, use the get_row_weight() and get_col_weight()
            * functions.
            */
            ENTRY_OBJ &insert_entry(int j, int i) {
                // std::cout<<i<<" "<<j<<std::endl;
                // Check if indices are within bounds
                if (j >= this->m || i >= this->n || j < 0 || i < 0) {
                    {
                        throw std::invalid_argument("Index i or j is out of bounds");
                    }
                }

                // Find the left and right entries in the jth row of the matrix
                auto left_entry_ptr = this->row_heads[j];
                auto right_entry_ptr = this->row_heads[j];
                for (auto &entry: reverse_iterate_row(j)) {
                    if (entry.col_index == i) {
                        // Entry already exists at this position
                        return entry;
                    }
                    if (entry.col_index > i) {
                        right_entry_ptr = &entry;
                    }
                    if (entry.col_index < i) {
                        left_entry_ptr = &entry;
                        break;
                    }
                }

                // Find the up and down entries in the ith column of the matrix
                auto up_entry_ptr = this->column_heads[i];
                auto down_entry_ptr = this->column_heads[i];
                for (auto &entry: this->reverse_iterate_column(i)) {
                    if (entry.row_index > j) {
                        down_entry_ptr = &entry;
                    }
                    if (entry.row_index < j) {
                        up_entry_ptr = &entry;
                        break;
                    }
                }

                // Create and link the new entry
                auto e_ptr = this->allocate_new_entry();
                node_count++;
                e_ptr->row_index = j;
                e_ptr->col_index = i;
                e_ptr->right = right_entry_ptr;
                e_ptr->left = left_entry_ptr;
                e_ptr->up = up_entry_ptr;
                e_ptr->down = down_entry_ptr;
                left_entry_ptr->right = e_ptr;
                right_entry_ptr->left = e_ptr;
                up_entry_ptr->down = e_ptr;
                down_entry_ptr->up = e_ptr;

                // Update row and column weights
                this->row_heads[e_ptr->row_index]->col_index--;
                this->column_heads[e_ptr->col_index]->col_index--;

                // Return a reference to the newly created entry
                return *e_ptr;
            }


            /**
             * Get an entry at row j and column i.
             *
             * @param j The row index
             * @param i The column index
             * @return a pointer to the entry, or a pointer to the corresponding column head.
             * @throws std::std::invalid_argument if j or i are out of bounds.
             */
            ENTRY_OBJ &get_entry(int j, int i) {
                if (j >= this->m || i >= this->n || j < 0 || i < 0) {
                    {
                        throw std::invalid_argument("Index i or j is out of bounds");
                    }
                }

                // Iterate over the column at index i and check each entry's row index.
                // If the row index matches j, return that entry.
                for (auto &e: this->reverse_iterate_column(i)) {
                    if (e.row_index == j) {
                        return e;
                    }
                }

                // If no entry is found, return the column head at index i.
                return *this->column_heads[i];
            }


            /**
             * Insert a new row at row_index with entries at column indices col_indices.
             *
             * @param row_index The index of the row to insert.
             * @param col_indices A vector of indices indicating which columns to insert entries into.
             * @return a pointer to the row head entry for the newly inserted row.
             */
            ENTRY_OBJ *insert_row(int row_index, std::vector<int> &col_indices) {
                // Insert an entry at each specified column index.
                for (auto j: col_indices) {
                    this->insert_entry(row_index, j);
                }

                // Return the row head at row_index.
                return this->row_heads[row_index];
            }


            /**
             * @brief Returns the coordinates of all non-zero entries in the matrix.
             *
             * @return Vector of vectors, where each inner vector represents the row and column indices of a non-zero entry.
             */
            std::vector<std::vector<int>> nonzero_coordinates() {

                std::vector<std::vector<int>> nonzero;

                // Initialize node count to 0
                this->node_count = 0;

                // Iterate through all rows and columns to find non-zero entries
                for (int i = 0; i < this->m; i++) {
                    for (auto &e: this->iterate_row(i)) {
                        // Increment node count and add non-zero entry coordinates to vector
                        this->node_count += 1;
                        std::vector<int> coord;
                        coord.push_back(e.row_index);
                        coord.push_back(e.col_index);
                        nonzero.push_back(coord);
                    }
                }

                // Return vector of non-zero entry coordinates
                return nonzero;

            }

            /**
             * Returns row adjacency list as vector of vectors.
             * @return
             */
            std::vector<std::vector<int>> row_adjacency_list() {
                std::vector<std::vector<int>> adj_list;
                for (int i = 0; i < this->m; i++) {
                    std::vector<int> row;
                    for (auto &e: this->iterate_row(i)) {
                        row.push_back(e.col_index);
                    }
                    adj_list.push_back(row);
                }
                return adj_list;
            }

            std::vector<std::vector<int>> col_adjacency_list() {
                std::vector<std::vector<int>> adj_list;
                for (int i = 0; i < this->n; i++) {
                    std::vector<int> col;
                    for (auto &e: this->iterate_column(i)) {
                        col.push_back(e.row_index);
                    }
                    adj_list.push_back(col);
                }
                return adj_list;
            }

            /**
             * Return a single column as 1D csc_matrix
             * @param col_index
             * @return
             */
            std::vector<int> get_column_csc(const std::size_t col_index) {
                std::vector<int> col;
                for (auto entry: this->iterate_column(col_index)) {
                    col.push_back(entry.row_index);
                }
                return col;
            }

            CsrMatrix to_csr_matrix() {
                CsrMatrix csr;
                csr.m = this->m;
                csr.n = this->n;
                csr.entry_count = this->entry_count();
                csr.row_adjacency_list = this->row_adjacency_list();
                return csr;
            }

            /**
            * Returns a vector of vectors, where each vector contains the column indices of the non-zero entries in a row.
            * @return A vector of vectors, where each vector contains the column indices of the non-zero entries in a row.
            */
            std::vector<std::vector<int>> nonzero_rows() {

                std::vector<std::vector<int>> nonzero;

                this->node_count = 0; //reset node_count to 0

                //iterate over the rows of the matrix
                for (int i = 0; i < this->m; i++) {
                    std::vector<int> row;

                    //iterate over the entries in the current row
                    for (auto &e: this->iterate_row(i)) {
                        this->node_count += 1; //increment node_count
                        row.push_back(
                                e.col_index); //add the column index of the current entry to the current row vector
                    }
                    nonzero.push_back(row); //add the current row vector to the vector of non-zero rows
                }

                return nonzero;

            }


            /**
             * @brief An iterator class that iterates over rows of a sparse matrix in a doubly linked list format.
             *
             * This class inherits from the Iterator class and is designed to work with SparseMatrixBase and its
             * subclasses. It is used to iterate over the rows of a sparse matrix in a doubly linked list format.
             *
             * @tparam ENTRY_OBJ The entry object class that the sparse matrix will use for its entries. This class
             * should contain fields for row index, column index, and value.
             */
            class RowIterator {
            public:
                using iterator_category = std::forward_iterator_tag;
                using difference_type = std::ptrdiff_t;
                SparseMatrixBase<ENTRY_OBJ> *matrix;
                int it_count;
                int entry_count;
                ENTRY_OBJ *e;

                RowIterator(SparseMatrixBase<ENTRY_OBJ> *mat, int i) {
                    matrix = mat;
                    entry_count = matrix->get_row_degree(i);
                    it_count = 0;
                    e = matrix->row_heads[i];
                }

                ~RowIterator() = default;;

                RowIterator &end() {
                    return *this;
                }

                RowIterator &begin() {
                    e = e->right;
                    ++it_count;
                    return *this;
                }

                ENTRY_OBJ &operator*() {
                    return *e;
                }

                // ENTRY_OBJ& operator*() {
                //     return *e;
                // }
                friend bool operator==(const RowIterator &a, const RowIterator &b) {
                    return a.it_count > b.entry_count;
                };

                friend bool operator!=(const RowIterator &a, const RowIterator &b) {
                    return a.it_count <= b.entry_count;
                };


                RowIterator &operator++() {
                    e = e->right;
                    ++it_count;
                    return *this;
                }

            };

            /**
             * @brief A reverse iterator for iterating over the rows of a SparseMatrixBase
             *
             * This iterator inherits from the `Iterator` base class using the CRTP pattern.
             * It is designed to be used with a SparseMatrixBase object to iterate over the rows
             * of the matrix in reverse order. It iterates over the rows in reverse order by
             * starting at the `head->left` entry and moving to the left using the `operator++()`
             * method. It can be indexed to start at any row of the matrix using the `operator[]`
             * method.
             *
             * @tparam ENTRY_OBJ The entry object class that the iterator will use for its entries.
             */
            class ReverseRowIterator {
            public:
                using iterator_category = std::forward_iterator_tag;
                using difference_type = std::ptrdiff_t;
                SparseMatrixBase<ENTRY_OBJ> *matrix;
                int it_count;
                int entry_count;
                ENTRY_OBJ *e;

                ReverseRowIterator(SparseMatrixBase<ENTRY_OBJ> *mat, int i) {
                    matrix = mat;
                    entry_count = matrix->get_row_degree(i);
                    it_count = 0;
                    e = matrix->row_heads[i];
                }

                ~ReverseRowIterator() = default;;

                ReverseRowIterator &end() {
                    return *this;
                }

                ReverseRowIterator &begin() {
                    e = e->left;
                    ++it_count;
                    return *this;
                }

                ENTRY_OBJ &operator*() {
                    return *e;
                }

                friend bool operator==(const ReverseRowIterator &a, const ReverseRowIterator &b) {
                    return a.it_count > b.entry_count;
                };

                friend bool operator!=(const ReverseRowIterator &a, const ReverseRowIterator &b) {
                    return a.it_count <= b.entry_count;
                };


                ReverseRowIterator &operator++() {
                    e = e->left;
                    ++it_count;
                    return *this;
                }

            };


            /**
             * @brief A forward iterator class that iterates over columns of a sparse matrix in a doubly linked list format.
             *
             * This class inherits from the Iterator class and is designed to work with SparseMatrixBase and its
             * subclasses. It is used to iterate over the columns of a sparse matrix in a doubly linked list format.
             *
             * @tparam ENTRY_OBJ The entry object class that the sparse matrix will use for its entries. This class
             * should contain fields for row index, column index, and value.
             */
            class ColumnIterator {
            public:
                using iterator_category = std::forward_iterator_tag;
                using difference_type = std::ptrdiff_t;
                SparseMatrixBase<ENTRY_OBJ> *matrix;
                int it_count;
                int entry_count;
                ENTRY_OBJ *e;

                ColumnIterator(SparseMatrixBase<ENTRY_OBJ> *mat, int i) {
                    matrix = mat;
                    entry_count = matrix->get_col_degree(i);
                    it_count = 0;
                    e = matrix->column_heads[i];
                }

                ~ColumnIterator() = default;;

                ColumnIterator &end() {
                    return *this;
                }

                ColumnIterator &begin() {
                    e = e->down;
                    ++it_count;
                    return *this;
                }

                ENTRY_OBJ &operator*() {
                    return *e;
                }

                friend bool operator==(const ColumnIterator &a, const ColumnIterator &b) {
                    return a.it_count > b.entry_count;
                };

                friend bool operator!=(const ColumnIterator &a, const ColumnIterator &b) {
                    return a.it_count <= b.entry_count;
                };


                ColumnIterator &operator++() {
                    e = e->down;
                    ++it_count;
                    return *this;
                }

            };

            /**
             * @brief A reverse iterator class that iterates over rows of a sparse matrix in a doubly linked list format.
             *
             * This class inherits from the Iterator class and is designed to work with SparseMatrixBase and its
             * subclasses. It is used to iterate over the rows of a sparse matrix in a doubly linked list format
             * starting from the rightmost element in the row.
             *
             * @tparam ENTRY_OBJ The entry object class that the sparse matrix will use for its entries. This class
             * should contain fields for row index, column index, and value.
             */
            class ReverseColumnIterator {
            public:
                using iterator_category = std::forward_iterator_tag;
                using difference_type = std::ptrdiff_t;
                SparseMatrixBase<ENTRY_OBJ> *matrix;
                int it_count;
                int entry_count;
                ENTRY_OBJ *e;

                ReverseColumnIterator(SparseMatrixBase<ENTRY_OBJ> *mat, int i) {
                    matrix = mat;
                    entry_count = matrix->get_col_degree(i);
                    it_count = 0;
                    e = matrix->column_heads[i];
                }

                ~ReverseColumnIterator() = default;;

                ReverseColumnIterator &end() {
                    return *this;
                }

                ReverseColumnIterator &begin() {
                    e = e->up;
                    ++it_count;
                    return *this;
                }

                ENTRY_OBJ &operator*() {
                    return *e;
                }

                friend bool operator==(const ReverseColumnIterator &a, const ReverseColumnIterator &b) {
                    return a.it_count > b.entry_count;
                };

                friend bool operator!=(const ReverseColumnIterator &a, const ReverseColumnIterator &b) {
                    return a.it_count <= b.entry_count;
                };


                ReverseColumnIterator &operator++() {
                    e = e->up;
                    ++it_count;
                    return *this;
                }

            };


            /**
             * @brief Returns an iterator that iterates over the given row of the sparse matrix in a forward direction
             *
             * @param i The row index of the matrix to iterate over
             * @throws std::invalid_argument If the given index is out of bounds
             * @return RowIterator An iterator object that iterates over the given row
             */
            RowIterator iterate_row(int i) {
                if (i < 0 || i >= m) {
                    throw std::invalid_argument("Iterator index out of bounds");
                }
                return RowIterator(this, i);
            }

            /**
             * @brief Returns an iterator that iterates over the given row of the sparse matrix in a reverse direction
             *
             * @param i The row index of the matrix to iterate over
             * @throws std::invalid_argument If the given index is out of bounds
             * @return ReverseRowIterator An iterator object that iterates over the given row in reverse
             */
            ReverseRowIterator reverse_iterate_row(int i) {
                if (i < 0 || i >= m) {
                    throw std::invalid_argument("Iterator index out of bounds");
                }
                return ReverseRowIterator(this, i);
            }

            /**
             * @brief Returns an iterator that iterates over the given column of the sparse matrix in a forward direction
             *
             * @param i The column index of the matrix to iterate over
             * @throws std::invalid_argument If the given index is out of bounds
             * @return ColumnIterator An iterator object that iterates over the given column
             */
            ColumnIterator iterate_column(int i) {
                if (i < 0 || i >= n) {
                    throw std::invalid_argument("Iterator index out of bounds");
                }
                return ColumnIterator(this, i);
            }

            /**
             * @brief Returns an iterator that iterates over the given column of the sparse matrix in a reverse direction
             *
             * @param i The column index of the matrix to iterate over
             * @throws std::invalid_argument If the given index is out of bounds
             * @return ReverseColumnIterator An iterator object that iterates over the given column in reverse
             */
            ReverseColumnIterator reverse_iterate_column(int i) {
                if (i < 0 || i >= n) {
                    throw std::invalid_argument("Iterator index out of bounds");
                }
                return ReverseColumnIterator(this, i);
            }


        };


    }  // namespace sparse_matrix_base
}  // namespace ldpc


#endif