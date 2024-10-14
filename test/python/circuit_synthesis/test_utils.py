"""Test utility functions for the circuit synthesis module."""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pytest
from ldpc import mod2

from mqt.qecc.circuit_synthesis.state_prep import final_matrix_constraint
from mqt.qecc.circuit_synthesis.synthesis_utils import (
    gaussian_elimination_min_column_ops,
    gaussian_elimination_min_parallel_eliminations,
)

if TYPE_CHECKING:
    import numpy.typing as npt


class MatrixTest(NamedTuple):
    """Test matrix and expected results."""

    matrix: npt.NDArray[np.int8]
    min_col_ops: int
    max_parallel_steps: int


def check_correct_elimination(
    matrix: npt.NDArray[np.int8], final_matrix: npt.NDArray[np.int8], col_ops: list[tuple[int, int]]
) -> bool:
    """Check if the matrix is correctly eliminated."""
    matrix = matrix.copy()
    for op in col_ops:
        matrix[:, op[1]] = (matrix[:, op[0]] + matrix[:, op[1]]) % 2
    print(matrix)
    return np.array_equal(matrix, final_matrix)


def get_n_parallel_layers(ops: list[tuple[int, int]]) -> int:
    """Get the number of parallel layers in the elimination."""
    used_cols: set[int] = set()
    layer = 0
    for op in ops:
        if op[0] in used_cols or op[1] in used_cols:
            layer += 1
            used_cols = set()
        used_cols.add(op[0])
        used_cols.add(op[1])
    return layer


@pytest.fixture
def identity_matrix() -> MatrixTest:
    """Return a 4x4 identity matrix."""
    return MatrixTest(np.eye(4, dtype=np.int8), 0, 0)


@pytest.fixture
def full_matrix() -> MatrixTest:
    """Return a 4x4 matrix with all ones."""
    return MatrixTest(np.array([[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]], dtype=np.int8), 3, 2)


@pytest.mark.parametrize("test_vals", ["identity_matrix", "full_matrix"])
def test_min_column_ops(test_vals: MatrixTest, request) -> None:  # type: ignore[no-untyped-def]
    """Check correct number of column operations are returned."""
    fixture = request.getfixturevalue(test_vals)
    matrix = fixture.matrix
    min_col_ops = fixture.min_col_ops
    rank = mod2.rank(matrix)
    res = gaussian_elimination_min_column_ops(
        matrix, lambda checks: final_matrix_constraint(checks, rank), max_eliminations=fixture.min_col_ops
    )
    assert res is not None
    reduced, ops = res
    assert len(ops) == min_col_ops
    assert check_correct_elimination(matrix, reduced, ops)


@pytest.mark.parametrize("test_vals", ["identity_matrix", "full_matrix"])
def test_min_parallel_eliminations(test_vals: MatrixTest, request) -> None:  # type: ignore[no-untyped-def]
    """Check correct number of parallel eliminations are returned."""
    fixture = request.getfixturevalue(test_vals)
    matrix = fixture.matrix
    rank = mod2.rank(matrix)
    max_parallel_steps = fixture.max_parallel_steps
    res = gaussian_elimination_min_parallel_eliminations(
        matrix, lambda checks: final_matrix_constraint(checks, rank), max_parallel_steps=fixture.max_parallel_steps
    )
    assert res is not None
    reduced, ops = res

    assert check_correct_elimination(matrix, reduced, ops)

    n_parallel_layers = get_n_parallel_layers(ops)
    assert n_parallel_layers <= max_parallel_steps
