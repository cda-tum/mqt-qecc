#ifndef GF2LINALG_H
#define GF2LINALG_H

#include <vector>
#include <iterator>
#include "gf2sparse.hpp"
#include "sparse_matrix_util.hpp"

namespace ldpc {
    namespace gf2sparse_linalg {

        std::vector<int> NULL_INT_VECTOR = {};

        template<class ENTRY_OBJ = ldpc::gf2sparse::GF2Entry>
        class RowReduce {

        public:
            ldpc::gf2sparse::GF2Sparse<ENTRY_OBJ> &A;
            ldpc::gf2sparse::GF2Sparse<> L;
            ldpc::gf2sparse::GF2Sparse<> U;
            ldpc::gf2sparse::GF2Sparse<> P;
            std::vector<int> rows;
            std::vector<int> cols;
            std::vector<int> inv_rows;
            std::vector<int> inv_cols;
            std::vector<bool> pivots;
            std::vector<uint8_t> x;
            std::vector<uint8_t> b;

            bool FULL_REDUCE{};

            bool LU_ALLOCATED{};

            bool LOWER_TRIANGULAR{};

            int rank{};

            RowReduce() = default;

            explicit RowReduce(ldpc::gf2sparse::GF2Sparse<ENTRY_OBJ> &A) : A(A) {

                this->pivots.resize(this->A.n, false);
                this->cols.resize(this->A.n);
                this->rows.resize(this->A.m);
                this->inv_cols.resize(this->A.n);
                this->inv_rows.resize(this->A.m);
                this->x.resize(this->A.n);
                this->b.resize(this->A.m);
                this->LU_ALLOCATED = false;
                this->LOWER_TRIANGULAR = false;

            }

            void initialise() {

                this->pivots.resize(this->A.n, false);
                this->cols.resize(this->A.n);
                this->rows.resize(this->A.m);
                this->inv_cols.resize(this->A.n);
                this->inv_rows.resize(this->A.m);
                this->x.resize(this->A.n);
                this->b.resize(this->A.m);
                this->LU_ALLOCATED = false;
                this->LOWER_TRIANGULAR = false;

            }

            ~RowReduce() = default;

            void initialise_LU() {

                this->U.allocate(this->A.m, this->A.n);
                this->L.allocate(this->A.m, this->A.m);

                for (int i = 0; i < this->A.m; i++) {
                    for (auto &e: this->A.iterate_row(i)) {
                        this->U.insert_entry(e.row_index, e.col_index);
                    }
                    if (!this->LOWER_TRIANGULAR) {
                        this->L.insert_entry(i, i);
                    }
                }

                this->LU_ALLOCATED = true;

            }

            void set_column_row_orderings(std::vector<int> &cols = NULL_INT_VECTOR,
                                          std::vector<int> &rows = NULL_INT_VECTOR) {

                if (cols == NULL_INT_VECTOR) {
                    for (int i = 0; i < this->A.n; i++) {
                        this->cols[i] = i;
                        this->inv_cols[this->cols[i]] = i;
                    }
                } else {
                    if (cols.size() != this->A.n) {
                        {
                            throw std::invalid_argument("Input parameter `cols`\
                describing the row ordering is of the incorrect length");
                        }
                    }
                    // this->cols=cols;
                    for (int i = 0; i < this->A.n; i++) {
                        this->cols[i] = cols[i];
                        inv_cols[cols[i]] = i;
                    }
                }

                if (rows == NULL_INT_VECTOR) {
                    for (int i = 0; i < this->A.m; i++) {
                        this->rows[i] = i;
                        this->inv_rows[this->rows[i]] = i;
                    }
                } else {
                    if (rows.size() != this->A.m) {
                        {
                            throw std::invalid_argument("Input parameter `rows`\
                describing the row ordering is of the incorrect length");
                        }
                    }
                    // this->rows=rows;
                    for (int i = 0; i < this->A.m; i++) {
                        this->rows[i] = rows[i];
                        this->inv_rows[rows[i]] = i;
                    }
                }

            }


            int rref(bool full_reduce = false, bool lower_triangular = false, std::vector<int> &cols = NULL_INT_VECTOR,
                     std::vector<int> &rows = NULL_INT_VECTOR) {

                this->LOWER_TRIANGULAR = lower_triangular;
                this->FULL_REDUCE = full_reduce;
                this->set_column_row_orderings(cols, rows);
                this->initialise_LU();
                int max_rank = std::min(this->U.m, this->U.n);
                this->rank = 0;
                std::fill(this->pivots.begin(), this->pivots.end(), false);

                for (int j = 0; j < this->U.n; j++) {
                    int pivot_index = this->cols[j];
                    if (this->rank == max_rank) {
                        break;
                    }

                    bool PIVOT_FOUND = false;
                    int max_row_weight = std::numeric_limits<int>::max();
                    int swap_index = 0;
                    for (auto &e: this->U.iterate_column(pivot_index)) {
                        int row_index = e.row_index;
                        if (row_index < this->rank) {
                            continue;
                        }

                        //finds the total row weight across both L and U for the given row
                        int row_weight = this->U.get_row_degree(row_index) + this->L.get_row_degree(row_index);
                        if (row_weight < max_row_weight) {
                            swap_index = e.row_index;
                            max_row_weight = row_weight;
                        }
                        PIVOT_FOUND = true;
                        this->pivots[j] = true;
                    }

                    if (!PIVOT_FOUND) {
                        continue;
                    }

                    if (swap_index != this->rank) {
                        U.swap_rows(swap_index, this->rank);
                        L.swap_rows(swap_index, this->rank);
                        auto temp1 = this->rows[swap_index];
                        auto temp2 = this->rows[this->rank];
                        this->rows[swap_index] = temp2;
                        this->rows[this->rank] = temp1;
                    }

                    if (this->LOWER_TRIANGULAR) {
                        this->L.insert_entry(this->rank, this->rank);
                    }

                    std::vector<int> add_rows;
                    for (auto &e: this->U.iterate_column(pivot_index)) {
                        int row_index = e.row_index;
                        if (row_index > this->rank || row_index < this->rank && static_cast<bool>(this->FULL_REDUCE)) {
                            add_rows.push_back(row_index);
                        }
                    }

                    for (int row: add_rows) {
                        this->U.add_rows(row, this->rank);
                        if (this->LOWER_TRIANGULAR) {
                            {
                                this->L.insert_entry(row, this->rank);
                            }
                        } else {
                            {
                                this->L.add_rows(row, this->rank);
                            }
                        }
                    }
                    this->rank++;
                }

                int pivot_count = 0;
                int non_pivot_count = 0;
                auto orig_cols = this->cols;
                for (int i = 0; i < this->U.n; i++) {
                    if (this->pivots[i]) {
                        this->cols[pivot_count] = orig_cols[i];
                        this->inv_cols[this->cols[pivot_count]] = pivot_count;
                        pivot_count++;
                    } else {
                        this->cols[this->rank + non_pivot_count] = orig_cols[i];
                        this->inv_cols[this->cols[this->rank + non_pivot_count]] = this->rank + non_pivot_count;
                        non_pivot_count++;
                    }
                }
                return this->rank;
            }

            void build_p_matrix() {

                if (!this->LU_ALLOCATED || !this->LOWER_TRIANGULAR) {
                    this->rref(false, true);
                }
                this->P.allocate(this->A.m, this->A.m);
                for (int i = 0; i < this->A.m; i++) {
                    this->P.insert_entry(this->rows[i], i);
                }

            }

            std::vector<uint8_t> &lu_solve(std::vector<uint8_t> &y) {

                if (y.size() != this->U.m) {
                    {
                        throw std::invalid_argument("Input parameter `y` is of the incorrect length for lu_solve.");
                    }
                }

                /*
                Equation: Ax = y

                We use LU decomposition to arrange the above into the form:
                LU(Qx) = PAQ^T(Qx)=Py

                We can then solve for x using forward-backward substitution:
                1. Forward substitution: Solve Lb = Py for b
                2. Backward subsitution: Solve UQx = b for x
                */


                if (!this->LU_ALLOCATED || !this->LOWER_TRIANGULAR) {
                    this->rref(false, true);
                }

                std::fill(this->x.begin(), this->x.end(), 0);
                std::fill(this->b.begin(), this->b.end(), 0);

                // //Solves LUx=y
                // int row_sum;



                //First we solve Lb = y, where b = Ux
                //Solve Lb=y with forwared substitution
                for (int row_index = 0; row_index < this->rank; row_index++) {
                    int row_sum = 0;
                    for (auto &e: this->L.iterate_row(row_index)) {
                        row_sum ^= this->b[e.col_index];
                    }
                    this->b[row_index] = row_sum ^ y[this->rows[row_index]];
                }

                //Solve Ux = b with backwards substitution
                for (int row_index = (this->rank - 1); row_index >= 0; row_index--) {
                    int row_sum = 0;
                    for (auto &e: this->U.iterate_row(row_index)) {
                        row_sum ^= this->x[e.col_index];
                    }
                    this->x[this->cols[row_index]] = row_sum ^ this->b[row_index];
                }
                return this->x;
            }

            /**
             * Solves Ax=y for x, where A is the matrix represented by this object. Stops elimination when y is in the
             * image of A.
             * @param y_vec
             * @param cols
             * @param rows
             * @return
             */
            std::vector<uint8_t> &fast_solve(std::vector<uint8_t> y_vec, std::vector<int> &cols = NULL_INT_VECTOR,
                                             std::vector<int> &rows = NULL_INT_VECTOR) {

                this->LOWER_TRIANGULAR = true;
                this->FULL_REDUCE = false;
                this->set_column_row_orderings(cols, rows);
                this->initialise_LU();
                ldpc::gf2sparse::GF2Sparse y(y_vec.size(), 1);
                int count = 0;

                for (auto &i: y_vec) {
                    if (i) {
                        y.insert_entry(count, 0);
                    }
                    count++;
                }

                int max_rank = std::min(this->U.m, this->U.n);
                this->rank = 0;
                std::fill(this->pivots.begin(), this->pivots.end(), false);

                for (int j = 0; j < this->U.n; j++) {
                    int pivot_index = this->cols[j];
                    if (this->rank == max_rank) {
                        break;
                    }
                    bool PIVOT_FOUND = false;
                    int max_row_weight = std::numeric_limits<int>::max();
                    int swap_index = 0;
                    for (auto &e: this->U.iterate_column(pivot_index)) {
                        int row_index = e.row_index;
                        if (row_index < this->rank) {
                            continue;
                        }
                        //finds the total row weight across both L and U for the given row
                        int row_weight = this->U.get_row_degree(row_index) + this->L.get_row_degree(row_index);
                        if (row_weight < max_row_weight) {
                            swap_index = e.row_index;
                            max_row_weight = row_weight;
                        }
                        PIVOT_FOUND = true;
                        this->pivots[j] = true;
                    }
                    if (!PIVOT_FOUND) {
                        continue;
                    }
                    if (swap_index != this->rank) {
                        U.swap_rows(swap_index, this->rank);
                        L.swap_rows(swap_index, this->rank);
                        y.swap_rows(swap_index, this->rank);
                        auto temp1 = this->rows[swap_index];
                        auto temp2 = this->rows[this->rank];
                        this->rows[swap_index] = temp2;
                        this->rows[this->rank] = temp1;
                    }
                    if (this->LOWER_TRIANGULAR) {
                        this->L.insert_entry(this->rank, this->rank);
                    }
                    std::vector<int> add_rows;
                    for (auto &e: this->U.iterate_column(pivot_index)) {
                        int row_index = e.row_index;
                        if (row_index > this->rank || row_index < this->rank && static_cast<bool>(this->FULL_REDUCE)) {
                            add_rows.push_back(row_index);
                        }
                    }
                    for (int row: add_rows) {
                        this->U.add_rows(row, this->rank);
                        if (this->LOWER_TRIANGULAR) {
                            this->L.insert_entry(row, this->rank);
                        } else {
                            this->L.add_rows(row, this->rank);
                        }
                        y.add_rows(row, this->rank);
                    }
                    this->rank++;
                    bool in_image = true;
                    for (auto &e: y.iterate_column(0)) {
                        int row_index = e.row_index;
                        if (row_index >= this->rank) {
                            in_image = false;
                            break;
                        }
                    }
                    if (in_image) {
                        break;
                    }
                }

                int pivot_count = 0;
                int non_pivot_count = 0;
                auto orig_cols = this->cols;
                for (int i = 0; i < this->U.n; i++) {
                    if (this->pivots[i]) {
                        this->cols[pivot_count] = orig_cols[i];
                        this->inv_cols[this->cols[pivot_count]] = pivot_count;
                        pivot_count++;
                    } else {
                        this->cols[this->rank + non_pivot_count] = orig_cols[i];
                        this->inv_cols[this->cols[this->rank + non_pivot_count]] = this->rank + non_pivot_count;
                        non_pivot_count++;
                    }
                }
                return this->lu_solve(y_vec);
            }


            auto
            rref_vrs(bool full_reduce = false, bool lower_triangular = false, std::vector<int> &cols = NULL_INT_VECTOR,
                     std::vector<int> &rows = NULL_INT_VECTOR) {

                if (lower_triangular) {
                    this->LOWER_TRIANGULAR = true;
                }
                this->set_column_row_orderings(cols, rows);
                this->initialise_LU();
                int max_rank = std::min(this->U.m, this->U.n);
                this->rank = 0;
                std::fill(this->pivots.begin(), this->pivots.end(), false);

                for (int j = 0; j < this->U.n; j++) {

                    int pivot_index = this->cols[j];

                    if (this->rank == max_rank) {
                        break;
                    }


                    bool PIVOT_FOUND = false;
                    int max_row_weight = std::numeric_limits<int>::max();
                    int swap_index = 0;
                    for (auto &e: this->U.iterate_column(pivot_index)) {
                        // int row_index = e.row_index;

                        if (this->inv_rows[e.row_index] < this->rank) {
                            continue;
                        }

                        int row_weight = this->U.get_row_degree(e.row_index) + this->L.get_row_degree(e.row_index);
                        if (this->inv_rows[e.row_index] >= this->rank && row_weight < max_row_weight) {
                            swap_index = this->inv_rows[e.row_index];
                            max_row_weight = row_weight;
                        }
                        PIVOT_FOUND = true;
                        this->pivots[j] = true;
                    }

                    if (!PIVOT_FOUND) {
                        continue;
                    }

                    if (swap_index != this->rank) {
                        // cout<<"Swapping rows "<<swap_index<<" and "<<this->rank<<endl;
                        // U.swap_rows(swap_index,this->inv_rows[this->rank]);
                        // L.swap_rows(swap_index,this->inv_rows[this->rank]);
                        auto temp1 = this->rows[swap_index];
                        auto temp2 = this->rows[this->rank];
                        this->rows[this->rank] = temp1;
                        this->rows[swap_index] = temp2;
                        this->inv_rows[temp1] = this->rank;
                        this->inv_rows[temp2] = swap_index;

                    }

                    if (this->LOWER_TRIANGULAR) {
                        this->L.insert_entry(this->rows[this->rank], this->rank);
                    }
                    // cout<<"Lower triangular: "<<endl;;
                    // print_sparse_matrix(*this->L);
                    // cout<<endl;


                    std::vector<int> add_rows;
                    for (auto &e: this->U.iterate_column(pivot_index)) {
                        // int row_index = e.row_index;
                        if (this->inv_rows[e.row_index] > this->rank ||
                            this->inv_rows[e.row_index] < this->rank && full_reduce) {
                            add_rows.push_back(e.row_index);
                        }
                    }

                    for (int row: add_rows) {
                        this->U.add_rows(row, this->rows[this->rank]);
                        if (lower_triangular) {
                            {
                                this->L.insert_entry(row, this->rank);
                            }
                        } else {
                            {
                                this->L.add_rows(row, this->rows[this->rank]);
                            }
                        }
                        // cout<<"Adding row "<<row<<" to row "<<this->rows[this->rank]<<endl;
                        // print_sparse_matrix(*this->U);
                        // cout<<endl;
                    }

                    this->rank++;

                }

                int pivot_count = 0;
                int non_pivot_count = 0;
                auto orig_cols = this->cols;
                for (int i = 0; i < this->U.n; i++) {
                    if (this->pivots[i]) {
                        this->cols[pivot_count] = orig_cols[i];
                        this->inv_cols[this->cols[pivot_count]] = pivot_count;
                        pivot_count++;
                    } else {
                        this->cols[this->rank + non_pivot_count] = orig_cols[i];
                        this->inv_cols[this->cols[this->rank + non_pivot_count]] = this->rank + non_pivot_count;
                        non_pivot_count++;
                    }
                }

                return this->U;

            }

            std::vector<uint8_t> &lu_solve_vrs(std::vector<uint8_t> &y) {

                if (y.size() != this->U.m) {
                    {
                        throw std::invalid_argument("Input parameter `y` is of the incorrect length for lu_solve.");
                    }
                }

                /*
                Equation: Ax = y

                We use LU decomposition to arrange the above into the form:
                LU(Qx) = PAQ^T(Qx)=Py
                We can then solve for x using forward-backward substitution:
                1. Forward substitution: Solve Lb = Py for b
                2. Backward subsitution: Solve UQx = b for x
                */


                if (!this->LU_ALLOCATED || !this->LOWER_TRIANGULAR) {
                    this->rref(false, true);
                }

                std::fill(this->x.begin(), this->x.end(), 0);
                std::fill(this->b.begin(), this->b.end(), 0);

                //Solves LUx=y
                int row_sum = 0;



                //First we solve Lb = y, where b = Ux
                //Solve Lb=y with forwared substitution
                for (int row_index = 0; row_index < this->rank; row_index++) {
                    row_sum = 0;
                    for (auto &e: L.iterate_row(this->rows[row_index])) {
                        row_sum ^= b[e.col_index];
                    }
                    b[row_index] = row_sum ^ y[this->rows[row_index]];
                }





                //Solve Ux = b with backwards substitution
                for (int row_index = (rank - 1); row_index >= 0; row_index--) {
                    row_sum = 0;
                    for (auto &e: U.iterate_row(this->rows[row_index])) {
                        row_sum ^= x[e.col_index];
                    }
                    x[this->cols[row_index]] = row_sum ^ b[row_index];
                }

                return x;

            }


        };

        template<class GF2MATRIX>
        std::vector<int> pivot_columns(GF2MATRIX &mat) {

            std::vector<int> swap_rows;
            std::vector<std::vector<int>> add_rows;
            std::vector<int> pivot_cols;
            std::vector<int> rows;

            for (int i = 0; i < mat.m; i++) {
                rows.push_back(i);
            }

            std::vector<uint8_t> rr_col(mat.m, 0);

            int max_rank = std::min(mat.m, mat.n);

            int rank = 0;

            for (int col = 0; col < mat.n; col++) {

                std::fill(rr_col.begin(), rr_col.end(), 0);

                for (auto &e: mat.iterate_column(col)) {
                    rr_col[e.row_index] = 1;
                }

                for (int i = 0; i < rank; i++) {

                    std::swap(rr_col[i], rr_col[swap_rows[i]]);
                    if (rr_col[i] == 1) {
                        for (int add_row: add_rows[i]) {
                            rr_col[add_row] ^= 1;
                        }
                    }

                }

                bool PIVOT_FOUND = false;

                for (int i = rank; i < mat.m; i++) {
                    if (rr_col[i] == 1) {
                        PIVOT_FOUND = true;
                        swap_rows.push_back(i);
                        pivot_cols.push_back(col);
                        break;
                    }
                }

                if (!PIVOT_FOUND) {
                    continue;
                }

                std::swap(rr_col[rank], rr_col[swap_rows[rank]]);

                add_rows.emplace_back();

                for (int i = rank + 1; i < mat.m; i++) {
                    if (rr_col[i] == 1) {
                        rr_col[i] ^= 1;
                        add_rows[rank].push_back(i);
                    }
                }

                rank++;

                if (rank == max_rank) {
                    break;
                }

            }

            return pivot_cols;

        }


        template<class GF2MATRIX>
        std::vector<int> pivot_rows(GF2MATRIX &mat) {

            std::vector<int> swap_cols;
            std::vector<std::vector<int>> add_cols;
            std::vector<int> pivot_rws;
            std::vector<int> rows;

            for (int i = 0; i < mat.n; i++) {
                rows.push_back(i);
            }

            std::vector<uint8_t> cr_row(mat.n, 0);

            int max_rank = std::min(mat.m, mat.n);

            int rank = 0;

            for (int row = 0; row < mat.m; row++) {

                std::fill(cr_row.begin(), cr_row.end(), 0);

                for (auto &e: mat.iterate_row(row)) {
                    cr_row[e.col_index] = 1;
                }

                for (int i = 0; i < rank; i++) {

                    std::swap(cr_row[i], cr_row[swap_cols[i]]);
                    if (cr_row[i] == 1) {
                        for (int add_col: add_cols[i]) {
                            cr_row[add_col] ^= 1;
                        }
                    }

                }

                bool PIVOT_FOUND = false;

                for (int i = rank; i < mat.n; i++) {
                    if (cr_row[i] == 1) {
                        PIVOT_FOUND = true;
                        swap_cols.push_back(i);
                        pivot_rws.push_back(row);
                        break;
                    }
                }

                if (!PIVOT_FOUND) {
                    continue;
                }

                std::swap(cr_row[rank], cr_row[swap_cols[rank]]);

                add_cols.emplace_back();

                for (int i = rank + 1; i < mat.n; i++) {
                    if (cr_row[i] == 1) {
                        cr_row[i] ^= 1;
                        add_cols[rank].push_back(i);
                    }
                }

                rank++;

                if (rank == max_rank) {
                    break;
                }

            }

            return pivot_rws;

        }


        template<class GF2MATRIX>
        std::vector<std::vector<int>> kernel2(GF2MATRIX &mat) {

            std::vector<int> swap_cols;
            std::vector<std::vector<int>> add_cols;
            std::vector<int> pivot_rws;
            std::vector<int> rows;

            for (int i = 0; i < mat.n; i++) {
                rows.push_back(i);
            }

            std::vector<uint8_t> cr_row(mat.n, 0);

            int max_rank = std::min(mat.m, mat.n);

            int rank = 0;

            for (int row = 0; row < mat.m; row++) {

                std::fill(cr_row.begin(), cr_row.end(), 0);

                for (auto &e: mat.iterate_row(row)) {
                    cr_row[e.col_index] = 1;
                }

                for (int i = 0; i < rank; i++) {

                    std::swap(cr_row[i], cr_row[swap_cols[i]]);
                    if (cr_row[i] == 1) {
                        for (int add_col: add_cols[i]) {
                            cr_row[add_col] ^= 1;
                        }
                    }

                }

                bool PIVOT_FOUND = false;

                for (int i = rank; i < mat.n; i++) {
                    if (cr_row[i] == 1) {
                        PIVOT_FOUND = true;
                        swap_cols.push_back(i);
                        pivot_rws.push_back(row);
                        break;
                    }
                }

                if (!PIVOT_FOUND) {
                    continue;
                }

                std::swap(cr_row[rank], cr_row[swap_cols[rank]]);

                add_cols.emplace_back();

                for (int i = rank + 1; i < mat.n; i++) {
                    if (cr_row[i] == 1) {
                        cr_row[i] ^= 1;
                        add_cols[rank].push_back(i);
                    }
                }

                rank++;

                if (rank == max_rank) {
                    break;
                }

            }

            std::vector<std::vector<int>> ker;

            for (int j = 0; j < mat.n; j++) {

                ker.emplace_back();

                std::fill(cr_row.begin(), cr_row.end(), 0);

                cr_row[j] = 1;

                for (int i = 0; i < rank; i++) {

                    std::swap(cr_row[i], cr_row[swap_cols[i]]);
                    if (cr_row[i] == 1) {
                        for (int add_col: add_cols[i]) {
                            cr_row[add_col] ^= 1;
                        }
                    }

                }

                for (int i = rank; i < mat.n; i++) {
                    if (cr_row[i] == 1) {
                        ker[j].push_back(i - rank);
                    }
                }

            }

            return ker;

        }


        template<class GF2MATRIX>
        int rank2(GF2MATRIX &mat) {
            return pivot_rows(mat).size();
        }


        template<class GF2MATRIX>
        std::vector<std::vector<int>> kernel_adjacency_list(GF2MATRIX &mat) {


            auto matT = mat.transpose();

            // cout<<"Transpose of input matrix: "<<endl;

            auto rr = RowReduce<>(matT);
            rr.rref(false, false);

            // cout<<"Rref done"<<endl;

            int rank = rr.rank;
            int n = mat.n;
            int k = n - rank;
            // auto ker = GF2MATRIX::New(k,n);

            std::vector<std::vector<int>> ker;
            ker.resize(k, std::vector<int>{});

            for (int i = rank; i < n; i++) {
                for (auto &e: rr.L.iterate_row(i)) {
                    ker[i - rank].push_back(e.col_index);
                }
            }

            return ker;

        }


        template<class GF2MATRIX>
        GF2MATRIX kernel(GF2MATRIX &mat) {

            auto adj_list = kernel_adjacency_list(mat);
            auto ker = GF2MATRIX(adj_list.size(), mat.n);
            ker.csr_insert(adj_list);

            return ker;

        }

//cython helper
        ldpc::sparse_matrix_base::CsrMatrix cy_kernel(ldpc::gf2sparse::GF2Sparse<ldpc::gf2sparse::GF2Entry> *mat) {
            return kernel(*mat).to_csr_matrix();
        }

        template<class GF2MATRIX>
        int rank(GF2MATRIX &mat) {
            auto rr = RowReduce(mat);
            rr.rref(false, false);
            return rr.rank;
        }

        template<class GF2MATRIX>
        GF2MATRIX row_complement_basis(GF2MATRIX &mat) {
            auto matT = mat.transpose();

            auto id_mat = GF2MATRIX(matT.m, matT.m, matT.m);
            for (int i = 0; i < matT.m; i++) {
                id_mat.insert_entry(i, i);
            }


            // print_sparse_matrix(id_mat);

            auto mat_aug = ldpc::gf2sparse::hstack(matT, id_mat);

            auto rr = RowReduce<>(mat_aug);
            rr.rref(false, false);


            std::vector<int> basis_rows;
            for (int i = 0; i < rr.rank; i++) {
                if (rr.cols[i] >= matT.n) {
                    basis_rows.push_back(rr.cols[i] - matT.n);
                }
            }

            GF2MATRIX basis = GF2MATRIX(basis_rows.size(), mat.n);
            for (int i = 0; i < basis_rows.size(); i++) {
                basis.insert_entry(i, basis_rows[i]);
            }

            // print_sparse_matrix(basis);
            // std::cout<<std::endl;

            return basis;


        }

//cython helper
        ldpc::sparse_matrix_base::CsrMatrix
        cy_row_complement_basis(ldpc::gf2sparse::GF2Sparse<ldpc::gf2sparse::GF2Entry> *mat) {
            return row_complement_basis(*mat).to_csr_matrix();
        }


    }//end namespace gf2sparse_linalg
}//end namespace ldpc

// typedef gf2sparse_linalg::RowReduce<ldpc::gf2sparse::ldpc::gf2sparse::GF2Sparse<ldpc::gf2sparse::GF2Entry>> cy_row_reduce;

#endif