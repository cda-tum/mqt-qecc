"""Constructions for (quantum) Hamming codes."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .css_code import CSSCode

if TYPE_CHECKING:
    import numpy.typing as npt


def hamming_code_checks(r: int) -> npt.NDArray[np.int8]:
    """Return the check matrix for the [2^r-1, 2^r-r-1, 3] Hamming code."""
    n = 2**r - 1
    h = np.zeros((r, n), dtype=int)
    # columns are all binary strings up to 2^r
    for i in range(1, n + 1):
        h[:, i - 1] = np.array([int(x) for x in f"{i:b}".zfill(r)])

    return h


def construct_quantum_hamming_code(r: int) -> CSSCode:
    """Return the [2^r, 2^r-r-1, 3] quantum Hamming code."""
    h = hamming_code_checks(r)
    return CSSCode(3, h, h)
