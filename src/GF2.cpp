/// Originally taken from the LDPC library by quantumgizmos released under the MIT license.
/// https://github.com/quantumgizmos/ldpc_v2/blob/cb85ee8601ee3fe59482985ab20840525c452a98/src_cpp/gf2dense.hpp

#include "GF2.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

void PluDecomposition::rref() {
    rows.reserve(rowCount);
    for (std::size_t i = 0; i < rowCount; i++) {
        rows.emplace_back(i);
    }

    const auto maxRank = std::min(rowCount, colCount);

    lower.resize(rowCount, std::vector<uint64_t>{});

    for (std::size_t colIdx = 0; colIdx < colCount; colIdx++) {
        eliminateColumn(colIdx);
        if (matrixRank == maxRank) {
            break;
        }
    }
}

bool PluDecomposition::eliminateColumn(std::size_t colIdx) {
    auto rrCol = std::vector<bool>(rowCount, false);

    for (auto rowIndex : cscMat[colIdx]) {
        rrCol[rowIndex] = true;
    }
    // apply previous operations to current column
    for (std::size_t i = 0; i < matrixRank; i++) {
        swap(rrCol[i], rrCol[swapRows[i]]);
        if (rrCol[i]) {
            // if row elem is one, do elimination for current column below the pivot
            // elimination operations to apply are stored in the `elimination_rows` attribute,
            for (std::size_t const rowIdx : eliminationRows[i]) {
                rrCol[rowIdx] = !rrCol[rowIdx];
            }
        }
    }
    bool pivotFound = false;
    for (auto i = matrixRank; i < rowCount; i++) {
        if (rrCol[i]) {
            pivotFound = true;
            swapRows.emplace_back(i);
            pivotCols.emplace_back(colIdx);
            break;
        }
    }

    // if no pivot was found, we go to next column
    if (!pivotFound) {
        return false;
    }

    swap(rrCol[matrixRank], rrCol[swapRows[matrixRank]]);
    std::swap(rows[matrixRank], rows[swapRows[matrixRank]]);
    eliminationRows.emplace_back();

    std::swap(lower[matrixRank], lower[swapRows[matrixRank]]);
    lower[matrixRank].emplace_back(matrixRank);

    for (auto i = matrixRank + 1; i < rowCount; i++) {
        if (rrCol[i]) {
            eliminationRows[matrixRank].emplace_back(i);
            lower[i].emplace_back(matrixRank);
        }
    }

    upper.emplace_back();
    for (std::size_t i = 0; i <= matrixRank; i++) {
        if (rrCol[i]) {
            upper[i].emplace_back(colIdx);
        }
    }

    matrixRank++;
    return true;
}

std::vector<uint64_t> PluDecomposition::luSolve(const std::vector<uint64_t>& y) {
    /*
    Equation: Ax = y

    We use LU decomposition to arrange the above into the form:
    LU(Qx) = PAQ^T(Qx)=Py

    We can then solve for x using forward-backward substitution:
    1. Forward substitution: Solve Lb = Py for b
    2. Backward substitution: Solve UQx = b for x
    */
    if (y.size() != rowCount) {
        throw std::invalid_argument("Input parameter `y` is of the incorrect length for luSolve.");
    }
    auto x = std::vector<uint64_t>(colCount, 0);
    auto b = std::vector<uint64_t>(matrixRank, 0);
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
    for (auto rowIndex = static_cast<int64_t>(matrixRank) - 1; rowIndex >= 0; rowIndex--) {
        std::size_t rowSum = 0;
        for (std::size_t const colIndex : upper[static_cast<size_t>(rowIndex)]) {
            rowSum ^= x[colIndex];
        }
        x[pivotCols[static_cast<size_t>(rowIndex)]] = rowSum ^ b[static_cast<size_t>(rowIndex)];
    }
    return x;
}
