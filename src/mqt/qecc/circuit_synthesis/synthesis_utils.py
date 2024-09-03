"""Utility functions for synthesizing circuits."""

from __future__ import annotations

import functools
import logging
from typing import TYPE_CHECKING, Any

import multiprocess
import numpy as np
import z3
from ldpc import mod2
from qiskit import QuantumCircuit

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


def iterative_search_with_timeout(
    fun: Callable[[int], QuantumCircuit],
    min_param: int,
    max_param: int,
    min_timeout: int,
    max_timeout: int,
    param_factor: float = 2,
    timeout_factor: float = 2,
) -> None | tuple[None | QuantumCircuit, int]:
    """Geometrically increases the parameter and timeout until a result is found or the maximum timeout is reached.

    Args:
        fun: function to run with increasing parameters and timeouts
        min_param: minimum parameter to start with
        max_param: maximum parameter to reach
        min_timeout: minimum timeout to start with
        max_timeout: maximum timeout to reach
        param_factor: factor to increase the parameter by at each iteration
        timeout_factor: factor to increase the timeout by at each iteration
    """
    curr_timeout = min_timeout
    curr_param = min_param
    while curr_timeout <= max_timeout:
        while curr_param <= max_param:
            logging.info(f"Running iterative search with param={curr_param} and timeout={curr_timeout}")
            res = run_with_timeout(fun, curr_param, timeout=curr_timeout)
            if res is not None and (not isinstance(res, str) or res != "timeout"):
                return res, curr_param
            if curr_param == max_param:
                break

            curr_param = int(curr_param * param_factor)
            curr_param = min(curr_param, max_param)

        curr_timeout = int(curr_timeout * timeout_factor)
        curr_param = min_param
    return None, max_param


def heuristic_gaussian_elimination(
    matrix: npt.NDArray[np.int8], parallel_elimination: bool = True
) -> tuple[npt.NDArray[np.int8], list[tuple[int, int]]]:
    """Perform Gaussian elimination on the column space of a matrix using as few eliminations as possible.

    The algorithm utilizes a greedy heuristic to select the columns to eliminate in order to minimize the number of eliminations required.

    The matrix is reduced until there are exactly rnk(matrix) columns with non-zero entries.

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


def gaussian_elimination_min_column_ops(
    matrix: npt.NDArray[np.int8],
    termination_criteria: Callable[[Any], z3.BoolRef],
    max_eliminations: int,
) -> tuple[npt.NDArray[np.int8], list[tuple[int, int]]] | None:
    """Perform Gaussian elimination on the column space of a matrix using at most `max_eliminations` eliminations.

    The algorithm encodes the elimination into an SMT problem and uses Z3 to find the optimal solution.

    Args:
        matrix: The matrix to perform Gaussian elimination on.
        termination_criteria: A function that takes a boolean matrix as input and returns a Z3 boolean expression that is true if the matrix is considered reduced.
        max_eliminations: The maximum number of eliminations to perform.

    returns:
        The reduced matrix and a list of the elimination steps taken. The elimination steps are represented as tuples of the form (i, j) where i is the column being eliminated with and j is the column being eliminated.
    """
    n = matrix.shape[1]
    columns = np.array([
        [[z3.Bool(f"x_{d}_{i}_{j}") for j in range(n)] for i in range(matrix.shape[0])]
        for d in range(max_eliminations + 1)
    ])

    n_bits = int(np.ceil(np.log2(n)))
    targets = [z3.BitVec(f"target_{d}", n_bits) for d in range(max_eliminations)]
    controls = [z3.BitVec(f"control_{d}", n_bits) for d in range(max_eliminations)]
    s = z3.Solver()

    additions = np.array([
        [[z3.And(controls[d] == col_1, targets[d] == col_2) for col_2 in range(n)] for col_1 in range(n)]
        for d in range(max_eliminations)
    ])

    # create initial matrix
    columns[0, :, :] = matrix.astype(bool)
    s.add(_column_addition_constraint(columns, additions))

    for d in range(1, max_eliminations + 1):
        # two columns cannot be in two elimination steps at the same time
        s.add(controls[d - 1] != targets[d - 1])

        # control and target must be valid qubits

        if n and (n - 1) != 0 and not ((n & (n - 1) == 0) and n != 0):  # check if n is a power of 2 or 1 or 0
            s.add(z3.ULT(controls[d - 1], n))
            s.add(z3.ULT(targets[d - 1], n))

    # if column is not involved in any addition at certain depth, it is the same as the previous column
    for d in range(1, max_eliminations + 1):
        for col in range(n):
            s.add(z3.Implies(targets[d - 1] != col, symbolic_vector_eq(columns[d, :, col], columns[d - 1, :, col])))

    # assert that final check matrix has n-checks.shape[0] zero columns
    s.add(termination_criteria(columns))

    if s.check() == z3.sat:
        m = s.model()
        eliminations = [(m[controls[d]].as_long(), m[targets[d]].as_long()) for d in range(max_eliminations)]
        reduced = np.array([
            [bool(m[columns[max_eliminations][i][j]]) for j in range(n)] for i in range(matrix.shape[0])
        ]).astype(np.int8)  # type: npt.NDArray[np.int8]
        return reduced, eliminations

    return None


def gaussian_elimination_min_parallel_eliminations(
    matrix: npt.NDArray[np.int8], termination_criteria: Callable[[Any], z3.BoolRef], max_parallel_steps: int
) -> tuple[npt.NDArray[np.int8], list[tuple[int, int]]] | None:
    """Perform Gaussian elimination on the column space of a matrix using at most `max_parallel_steps` parallel column elimination steps.

    The algorithm encodes the elimination into a SAT problem and uses Z3 to find the optimal solution.

    Args:
        matrix: The matrix to perform Gaussian elimination on.
        termination_criteria: A function that takes a boolean matrix as input and returns a Z3 boolean expression that is true if the matrix is considered reduced.
        max_parallel_steps: The maximum number of parallel elimination steps to perform.

    returns:
        The reduced matrix and a list of the elimination steps taken. The elimination steps are represented as tuples of the form (i, j) where i is the column being eliminated with and j is the column being eliminated.
    """
    assert max_parallel_steps > 0, "max_parallel_steps should be greater than 0"
    columns = np.array([
        [[z3.Bool(f"x_{d}_{i}_{j}") for j in range(matrix.shape[1])] for i in range(matrix.shape[0])]
        for d in range(max_parallel_steps + 1)
    ])

    additions = np.array([
        [[z3.Bool(f"add_{d}_{i}_{j}") for j in range(matrix.shape[1])] for i in range(matrix.shape[1])]
        for d in range(max_parallel_steps)
    ])
    n_cols = matrix.shape[1]
    s = z3.Solver()

    # create initial matrix
    columns[0, :, :] = matrix.astype(bool)

    s.add(_column_addition_constraint(columns, additions))

    # qubit can be involved in at most one addition at each depth
    for d in range(max_parallel_steps):
        for col in range(n_cols):
            s.add(
                z3.PbLe(
                    [(additions[d, col_1, col], 1) for col_1 in range(n_cols) if col != col_1]
                    + [(additions[d, col, col_2], 1) for col_2 in range(n_cols) if col != col_2],
                    1,
                )
            )

    # if column is not involved in any addition at certain depth, it is the same as the previous column
    for d in range(1, max_parallel_steps + 1):
        for col in range(n_cols):
            s.add(
                z3.Implies(
                    z3.Not(
                        z3.Or(
                            list(np.delete(additions[d - 1, :, col], [col]))
                            + list(np.delete(additions[d - 1, col, :], [col]))
                        )
                    ),
                    symbolic_vector_eq(columns[d, :, col], columns[d - 1, :, col]),
                )
            )

    s.add(termination_criteria(columns))

    if s.check() == z3.sat:
        m = s.model()
        eliminations = [
            (i, j)
            for d in range(max_parallel_steps)
            for j in range(matrix.shape[1])
            for i in range(matrix.shape[1])
            if m[additions[d, i, j]]
        ]
        reduced = np.array([
            [bool(m[columns[max_parallel_steps, i, j]]) for j in range(matrix.shape[1])] for i in range(matrix.shape[0])
        ]).astype(np.int8)  # type: npt.NDArray[np.int8]
        return reduced, eliminations

    return None


def build_css_circuit_from_cnot_list(n: int, cnots: list[tuple[int, int]], hadamards: list[int]) -> QuantumCircuit:
    """Build a quantum circuit consisting of Hadamards followed by a layer of CNOTs from a list of CNOTs and a list of checks.

    Args:
        n: Number of qubits in the circuit.
        cnots: List of CNOTs to apply. Each CNOT is a tuple of the form (control, target).
        hadamards: List of qubits to apply Hadamards to.

    Returns:
        The quantum circuit.
    """
    circ = QuantumCircuit(n)
    circ.h(hadamards)
    for i, j in cnots:
        circ.cx(i, j)
    return circ


def _column_addition_constraint(
    columns: npt.NDArray[z3.BoolRef | bool],  # type: ignore[type-var]
    col_add_vars: npt.NDArray[z3.BoolRef],
) -> z3.BoolRef:
    assert len(columns.shape) == 3
    max_parallel_steps = col_add_vars.shape[0]  # type: ignore[unreachable]
    n_cols = col_add_vars.shape[2]

    constraints = []
    for d in range(1, max_parallel_steps + 1):
        for col_1 in range(n_cols):
            for col_2 in range(col_1 + 1, n_cols):
                col_sum = symbolic_vector_add(columns[d - 1, :, col_1], columns[d - 1, :, col_2])

                # encode col_2 += col_1
                add_col1_to_col2 = z3.Implies(
                    col_add_vars[d - 1, col_1, col_2],
                    z3.And(
                        symbolic_vector_eq(columns[d, :, col_2], col_sum),
                        symbolic_vector_eq(columns[d, :, col_1], columns[d - 1, :, col_1]),
                    ),
                )

                # encode col_1 += col_2
                add_col2_to_col1 = z3.Implies(
                    col_add_vars[d - 1, col_2, col_1],
                    z3.And(
                        symbolic_vector_eq(columns[d, :, col_1], col_sum),
                        symbolic_vector_eq(columns[d, :, col_2], columns[d - 1, :, col_2]),
                    ),
                )

                constraints.extend([add_col1_to_col2, add_col2_to_col1])

    return z3.And(constraints)


def symbolic_vector_eq(v1: npt.NDArray[z3.BoolRef | bool], v2: npt.NDArray[z3.BoolRef | bool]) -> z3.BoolRef:
    """Return assertion that two symbolic vectors should be equal."""
    constraints = [False for _ in v1]
    for i in range(len(v1)):
        # If one of the elements is a bool, we can simplify the expression
        v1_i_is_bool = isinstance(v1[i], (bool, np.bool_))
        v2_i_is_bool = isinstance(v2[i], (bool, np.bool_))
        if v1_i_is_bool:
            v1[i] = bool(v1[i])
            if v1[i]:
                constraints[i] = v2[i]
            else:
                constraints[i] = z3.Not(v2[i]) if not v2_i_is_bool else not v2[i]

        elif v2_i_is_bool:
            v2[i] = bool(v2[i])
            if v2[i]:
                constraints[i] = v1[i]
            else:
                constraints[i] = z3.Not(v1[i])
        else:
            constraints[i] = v1[i] == v2[i]
    return z3.And(constraints)


def odd_overlap(v_sym: npt.NDArray[z3.BoolRef | bool], v_con: npt.NDArray[np.int8]) -> z3.BoolRef:
    """Return True if the overlap of symbolic vector with constant vector is odd."""
    return z3.PbEq([(v_sym[i], 1) for i, c in enumerate(v_con) if c == 1], 1)


def symbolic_scalar_mult(v: npt.NDArray[np.int8], a: z3.BoolRef | bool) -> npt.NDArray[z3.BoolRef]:
    """Multiply a concrete vector by a symbolic scalar."""
    return np.array([a if s == 1 else False for s in v])


def symbolic_vector_add(
    v1: npt.NDArray[z3.BoolRef | bool], v2: npt.NDArray[z3.BoolRef | bool]
) -> npt.NDArray[z3.BoolRef | bool]:
    """Add two symbolic vectors."""
    v_new = [False for _ in range(len(v1))]
    for i in range(len(v1)):
        # If one of the elements is a bool, we can simplify the expression
        v1_i_is_bool = isinstance(v1[i], (bool, np.bool_))
        v2_i_is_bool = isinstance(v2[i], (bool, np.bool_))
        if v1_i_is_bool:
            v1[i] = bool(v1[i])
            if v1[i]:
                v_new[i] = z3.Not(v2[i]) if not v2_i_is_bool else not v2[i]
            else:
                v_new[i] = v2[i]

        elif v2_i_is_bool:
            v2[i] = bool(v2[i])
            if v2[i]:
                v_new[i] = z3.Not(v1[i])
            else:
                v_new[i] = v1[i]

        else:
            v_new[i] = z3.Xor(v1[i], v2[i])

    return np.array(v_new)


def optimal_elimination(
    matrix: npt.NDArray[np.int8],
    termination_criteria: Callable[[Any], z3.BoolRef],
    optimization_metric: str = "column_ops",
    min_param: int = 1,
    max_param: int = 10,
    min_timeout: int = 1,
    max_timeout: int = 3600,
) -> tuple[npt.NDArray[np.int8], list[tuple[int, int]]] | None:
    """Synthesize a state preparation circuit for a CSS code that minimizes the circuit w.r.t. some metric param according to prep_func.

    Args:
        matrix: The stabilizer matrix of the CSS code.
        termination_criteria: The termination criteria for when the matrix is considered reduced.
        optimization_metric: The metric to optimize the circuit w.r.t. to. Can be either "column_ops" or "parallel_ops".
        zero_state: Whether to start from the zero state.
        min_param: The minimum value of the metric parameter.
        max_param: The maximum value of the metric parameter.
        min_timeout: The minimum time to run one search iteration for.
        max_timeout: The maximum time to run one search iteration for.
    """
    if optimization_metric not in {"column_ops", "parallel_ops"}:
        msg = "Invalid optimization metric"
        raise ValueError(msg)

    opt_fun = {
        "column_ops": gaussian_elimination_min_column_ops,
        "parallel_ops": gaussian_elimination_min_parallel_eliminations,
    }[optimization_metric]

    fun = functools.partial(
        opt_fun,
        matrix,
        termination_criteria,
    )

    res = iterative_search_with_timeout(
        fun,
        min_param,
        max_param,
        min_timeout,
        max_timeout,
    )

    if res is None:
        return None
    reduced = res[0]
    if reduced is None:
        return None
    reduced, eliminations = reduced
    curr_param = res[1]

    logging.info(f"Solution found with param {curr_param}")
    # Solving a SAT instance is much faster than proving unsat in this case
    # so we iterate backwards until we find an unsat instance or hit a timeout
    logging.info("Trying to minimize param")
    while True:
        logging.info(f"Trying param {curr_param - 1}")
        opt_res = run_with_timeout(fun, curr_param - 1, timeout=max_timeout)
        if opt_res is None or (isinstance(opt_res, str) and opt_res == "timeout"):
            break
        assert not isinstance(opt_res, str)
        reduced, eliminations = opt_res
        curr_param -= 1

    logging.info(f"Optimal param: {curr_param}")
    return reduced, eliminations
