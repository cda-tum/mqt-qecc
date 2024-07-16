"""Construction of Quasi-Cyclic LDPC codes from https://arxiv.org/abs/2308.07915."""

from __future__ import annotations

import numpy as np


def quasi_cyclic_check_matrices(n: int) -> tuple[np.array, np.array]:
    """Construct the check matrices for a quasi-cyclic LDPC code."""
    if n == 72:
        return _cyclic_code_check_matrix(6, 6, ([3], [1, 2]), ([1, 2], [3]))
    if n == 90:
        return _cyclic_code_check_matrix(15, 3, ([9], [1, 2]), ([2, 7], [0]))
    if n == 108:
        return _cyclic_code_check_matrix(9, 6, ([3], [1, 2]), ([1, 2], [3]))
    if n == 144:
        return _cyclic_code_check_matrix(12, 6, ([3], [1, 2]), ([1, 2], [3]))
    if n == 288:
        return _cyclic_code_check_matrix(12, 12, ([3], [2, 7]), ([1, 2], [3]))
    msg = f"No quasi-cyclic code with n = {n}."
    raise ValueError(msg)


def _make_check_array(n):
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
    return np.array(checks)


def _shift_matrix(l: int) -> np.array:
    s = np.zeros((l, l))
    for i in range(l):
        s[i, (i + 1) % l] = 1
    return s


def _x_matrix(l: int, m: int) -> np.array:
    return np.kron(_shift_matrix(l), np.eye(m))


def _y_matrix(l: int, m: int) -> np.array:
    return np.kron(np.eye(l), _shift_matrix(m))


def _cyclic_matrix_sum(x_terms, y_terms, l, m):
    x = _x_matrix(l, m)
    y = _y_matrix(l, m)
    a = np.zeros(x.shape)
    for pow in x_terms:
        a += np.linalg.matrix_power(x, pow)
    for pow in y_terms:
        a += np.linalg.matrix_power(y, pow)

    return a


def _cyclic_code_check_matrix(
    m: int, l: int, a_terms: tuple[list, list], b_terms: tuple[list, list]
) -> tuple[np.array, np.array]:
    a = _cyclic_matrix_sum(a_terms[0], a_terms[1], l, m) % 2
    b = _cyclic_matrix_sum(b_terms[0], b_terms[1], l, m) % 2
    return np.hstack((a, b)), np.hstack((b.T, a.T))
