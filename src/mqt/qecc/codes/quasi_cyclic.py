"""Construction of Quasi-Cyclic LDPC codes from https://arxiv.org/abs/2308.07915."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .css_code import CSSCode, InvalidCSSCodeError

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt


def construct_quasi_cyclic_code(n: int) -> CSSCode:
    """Construct the check matrices for a quasi-cyclic LDPC code.

    Args:
        n: The number of qubits. Currently supported values are 72, 90, 108, 144, and 288.
    """
    if n == 72:
        x, z = _cyclic_code_check_matrix(6, 6, ([3], [1, 2]), ([1, 2], [3]))
        d = 6
    elif n == 90:
        x, z = _cyclic_code_check_matrix(15, 3, ([9], [1, 2]), ([2, 7], [0]))
        d = 10
    elif n == 108:
        x, z = _cyclic_code_check_matrix(9, 6, ([3], [1, 2]), ([1, 2], [3]))
        d = 10
    elif n == 144:
        x, z = _cyclic_code_check_matrix(12, 6, ([3], [1, 2]), ([1, 2], [3]))
        d = 12
    elif n == 288:
        x, z = _cyclic_code_check_matrix(12, 12, ([3], [2, 7]), ([1, 2], [3]))
        d = 18
    else:
        msg = f"No quasi-cyclic code with n = {n}."
        raise InvalidCSSCodeError(msg)
    return CSSCode(d, x, z)


def _make_check_array(n: int) -> npt.NDArray[np.int8]:
    checks = []
    for _ in range(n // 2 - 1):
        check = np.random.choice([0, 1], size=(n,), p=[1 - 4 / n, 4 / n])
        while np.sum(check) == 0:
            check = np.random.choice([0, 1], size=(n,), p=[1 - 4 / n, 4 / n])
        checks.append(check)
    while np.any(np.sum(np.array(checks), axis=0) == 0):
        checks = []
        for _ in range(n // 2 - 1):
            check = np.random.choice([0, 1], size=(n,), p=[1 - 4 / n, 4 / n])
            while np.sum(check) == 0:
                check = np.random.choice([0, 1], size=(n,), p=[1 - 4 / n, 4 / n])
            checks.append(check)
    return np.array(checks, dtype=np.int8)


def _shift_matrix(l: int) -> npt.NDArray[np.int8]:
    s = np.zeros((l, l), dtype=np.int8)
    for i in range(l):
        s[i, (i + 1) % l] = 1
    return s


def _x_matrix(l: int, m: int) -> npt.NDArray[np.int8]:
    return np.kron(_shift_matrix(l), np.eye(m))


def _y_matrix(l: int, m: int) -> npt.NDArray[np.int8]:
    return np.kron(np.eye(l), _shift_matrix(m))


def _cyclic_matrix_sum(x_terms, y_terms, l, m) -> npt.NDArray[np.int8]:
    x = _x_matrix(l, m)
    y = _y_matrix(l, m)
    a = np.zeros(x.shape)
    for pow in x_terms:
        a += np.linalg.matrix_power(x, pow)
    for pow in y_terms:
        a += np.linalg.matrix_power(y, pow)

    return a


def _cyclic_code_check_matrix(
    m: int, l: int, a_terms: tuple[list[int], list[int]], b_terms: tuple[list[int], list[int]]
) -> tuple[npt.NDArray[np.int8], npt.NDArray[np.int8]]:
    a = _cyclic_matrix_sum(a_terms[0], a_terms[1], l, m) % 2
    b = _cyclic_matrix_sum(b_terms[0], b_terms[1], l, m) % 2
    return np.hstack((a, b)), np.hstack((b.T, a.T))
