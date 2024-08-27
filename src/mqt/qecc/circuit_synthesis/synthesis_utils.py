"""Utility functions for synthesizing circuits."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import multiprocess
import numpy as np
from ldpc import mod2

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Callable

    import numpy.typing as npt


logger = logging.getLogger(__name__)


def run_with_timeout(func: Callable[[Any], Any], *args: Any, timeout: int = 10) -> Any | str | None:  # noqa: ANN401
    """Run a function with a timeout.

    If the function does not complete within the timeout, return None.

    Args:
        func: The function to run.
        args: The arguments to pass to the function.
        timeout: The maximum time to allow the function to run for in seconds.
    """
    manager = multiprocess.Manager()
    return_list = manager.list()
    p = multiprocess.Process(target=lambda: return_list.append(func(*args)))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
        return "timeout"
    return return_list[0]


def heuristic_gaussian_elimination(
    matrix: npt.NDArray[np.int8], parallel_elimination: bool = True
) -> tuple[npt.NDArray[np.int8], list[tuple[int, int]]]:
    """Perform Gaussian elimination on the column space of a matrix using as few eliminations as possible.

    The algorithm utilizes a greedy heuristic to select the columns to eliminate in order to minimize the number of eliminations required.

    Args:
        matrix: The matrix to perform Gaussian elimination on.
        parallel_elimination: Whether to prioritize elimination steps that act on disjoint columns.

    returns:
        The reduced matrix and a list of the elimination steps taken. The elimination steps are represented as tuples of the form (i, j) where i is the column being eliminated with and j is the column being eliminated.
    """
    matrix = matrix.copy()
    rank = mod2.rank(matrix)

    def is_reduced() -> bool:
        return bool(len(np.where(np.all(matrix == 0, axis=0))[0]) == matrix.shape[1] - rank)

    costs = np.array([
        [np.sum((matrix[:, i] + matrix[:, j]) % 2) for j in range(matrix.shape[1])] for i in range(matrix.shape[1])
    ])

    costs -= np.sum(matrix, axis=0)
    np.fill_diagonal(costs, 1)

    used_columns = []  # type: list[np.int_]
    eliminations = []  # type: list[tuple[int, int]]
    while not is_reduced():
        m = np.zeros((matrix.shape[1], matrix.shape[1]), dtype=bool)  # type: npt.NDArray[np.bool_]
        m[used_columns, :] = True
        m[:, used_columns] = True

        costs_unused = np.ma.array(costs, mask=m)  # type: ignore[no-untyped-call]
        if np.all(costs_unused >= 0) or len(used_columns) == matrix.shape[1]:  # no more reductions possible
            if used_columns == []:  # local minimum => get out by making matrix triangular
                logging.warning("Local minimum reached. Making matrix triangular.")
                matrix = mod2.reduced_row_echelon(matrix)[0]
                costs = np.array([
                    [np.sum((matrix[:, i] + matrix[:, j]) % 2) for j in range(matrix.shape[1])]
                    for i in range(matrix.shape[1])
                ])
                costs -= np.sum(matrix, axis=0)
                np.fill_diagonal(costs, 1)
            else:  # try to move onto the next layer
                used_columns = []
            continue

        i, j = np.unravel_index(np.argmin(costs_unused), costs.shape)
        eliminations.append((int(i), int(j)))

        if parallel_elimination:
            used_columns.append(i)
            used_columns.append(j)

        # update matrix
        matrix[:, j] = (matrix[:, i] + matrix[:, j]) % 2
        # update costs
        new_weights = np.sum((matrix[:, j][:, np.newaxis] + matrix) % 2, axis=0)
        costs[j, :] = new_weights - np.sum(matrix, axis=0)
        costs[:, j] = new_weights - np.sum(matrix[:, j])
        np.fill_diagonal(costs, 1)

    return matrix, eliminations
