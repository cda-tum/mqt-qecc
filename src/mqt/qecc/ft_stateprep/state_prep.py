"""Synthesizing state preparation circuits for CSS codes."""

from __future__ import annotations

from ldpc import mod2
# from code import CSSCode
import numpy as np
from qiskit import QuantumCircuit
import z3

import multiprocessing
import time

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt


def _build_circuit_from_list_and_checks(cnots: list[tuple], checks: npt.NDArray[np.int_], zero_state=True) -> QuantumCircuit:
    # Build circuit
    n = checks.shape[1]
    circ = QuantumCircuit(n)

    controls = [i for i in range(n) if np.sum(checks[:, i]) >= 1]
    if zero_state:
        for control in controls:
            circ.h(control)
    else:
        for i in range(n):
            if i not in controls:
                circ.h(i)

    for i, j in reversed(cnots):
        if not zero_state:
            i, j = j, i
        circ.cx(i, j)

    return circ


def heuristic_prep_circuit(code: CSSCode, optimize_depth: bool=True, zero_state: bool=True):
    """Return a circuit that prepares the +1 eigenstate of the code w.r.t. the Z or X basis.

    Args:
        code: The CSS code to prepare the state for.
        zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
    """
    checks = code.Hx.copy() if zero_state else code.Hz.copy()
    rank = mod2.rank(checks)
    
    def is_reduced():
        return len(np.where(np.all(checks==0, axis=0))[0]) == checks.shape[1]-rank
    
    costs = np.array([[np.sum((checks[:, i] + checks[:, j]) % 2) for j in range(checks.shape[1])] for i in range(checks.shape[1])])
    costs -= np.sum(checks, axis=0)
    np.fill_diagonal(costs, 1)

    used_qubits = []
    cnots = []
    while not is_reduced():
        m = np.zeros((checks.shape[1], checks.shape[1]), dtype=bool)
        m[used_qubits, :] = True
        m[:, used_qubits] = True

        costs_unused = np.ma.array(costs, mask=m)
        if np.all(costs_unused >= 0):  # no more reductions possible
            if used_qubits == []:  # local minimum => get out by making matrix triangular
                costs = np.array([[np.sum((checks[:, i] + checks[:, j]) % 2)
                                   for j in range(checks.shape[1])]
                                  for i in range(checks.shape[1])])
                costs -= np.sum(checks, axis=0)
                np.fill_diagonal(costs, 1)
                break
                checks = mod2.reduced_row_echelon(checks)[0]
            else:  # try to move onto the next layer
                used_qubits = []
            continue

        i, j = np.unravel_index(np.argmin(costs_unused), costs.shape)
        cnots.append((i, j))

        if optimize_depth:
            used_qubits.append(i)
            used_qubits.append(j)

        # update checks
        checks[:, j] = (checks[:, i] + checks[:, j]) % 2

        # update costs
        new_weights = np.sum((checks[:, j][:, np.newaxis] + checks) % 2, axis=0)
        costs[:, j] = new_weights - np.sum(checks, axis=0)
        costs[j, :] = new_weights - np.sum(checks[:, j])
        np.fill_diagonal(costs, 1)

    return _build_circuit_from_list_and_checks(cnots, checks, zero_state)


def _run_with_timeout(func, *args, timeout: int=10):
    """Run a function with a timeout. If the function does not complete within the timeout, return None.

    Args:
        func: The function to run.
        args: The arguments to pass to the function.
        timeout: The maximum time to allow the function to run for in seconds."""
    
    manager = multiprocessing.Manager()
    return_list = manager.list()
    p = multiprocessing.Process(target=lambda: return_list.append(func(*args)))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
        return "timeout"
    return return_list[0]


def _symbolic_scalar_mult(v: npt.NDArray[np.int_], a: z3.BoolRef | bool):
    """Multiply a concrete vector by a symbolic scalar."""
    return [a if s == 1 else False for s in v]


SymOrBool = z3.BoolRef | bool
SymVec = list[SymOrBool]


def _symbolic_vector_add(v1: SymVec, v2: SymVec):
    """Add two symbolic vectors."""
    if v1 is None:
        return v2
    if v2 is None:
        return v1

    v_new = [False for _ in range(len(v1))]
    for i in range(len(v1)):
        if isinstance(v1[i], bool):
            if v1[i]:
                v_new[i] = z3.Not(v2[i])
            else:
                v_new[i] = v2[i]

        elif isinstance(v2[i], bool):
            if v2[i]:
                v_new[i] = z3.Not(v1[i])
            else:
                v_new[i] = v1[i]

        else:
            v_new[i] = z3.Xor(v1[i], v2[i])

    return v_new


def _odd_overlap(v_sym: SymVec, v_con: npt.NDArray[np.int_]):
    """Return True if the overlap of symbolic vector with constant vector is odd."""
    return z3.PbEq([(v_sym[i], 1) for i, c in enumerate(v_con) if c == 1], 1)


def _generate_circ_with_bounded_depth(checks: npt.NDArray, max_depth) -> np.array:
    columns = np.array([[[z3.Bool(f'x_{d}_{i}_{j}')
                          for j in range(checks.shape[1])]
                         for i in range(checks.shape[0])]
                        for d in range(max_depth+1)])
    
    additions = np.array([[[z3.Bool(f'add_{d}_{i}_{j}')
                            for j in range(checks.shape[1])]
                           for i in range(checks.shape[1])]
                          for d in range(max_depth)])
    s = z3.Solver()

    # create initial matrix
    columns[0, :, :] = checks.astype(bool)
    
    # encode all possible column additions
    for d in range(1, max_depth+1):
        for col_1 in range(checks.shape[1]):
            for col_2 in range(col_1+1, checks.shape[1]):
                col_sum = _symbolic_vector_add(columns[d-1, :, col_1], columns[d-1, :, col_2])

                # encode col_2 += col_1
                s.add(z3.Implies(additions[d-1, col_1, col_2],
                                 z3.And([columns[d, i, col_2] == col_sum[i]
                                         for i in range(checks.shape[0])] +
                                        [columns[d, i, col_1] == columns[d-1, i, col_1]
                                         for i in range(checks.shape[0])])))
                # encode col_1 += col_2
                s.add(z3.Implies(additions[d-1, col_2, col_1],
                                 z3.And([columns[d, i, col_1] == col_sum[i]
                                         for i in range(checks.shape[0])] +
                                        [columns[d, i, col_2] == columns[d-1, i, col_2]
                                         for i in range(checks.shape[0])])))

    # at most one addition per column
    for d in range(max_depth):
        for col in range(checks.shape[1]):
            s.add(z3.PbLe([(additions[d, col_1, col], 1)
                           for col_1 in range(checks.shape[1])
                           if col != col_1] +
                          [(additions[d, col, col_2], 1)
                           for col_2 in range(checks.shape[1])
                           if col != col_2], 1
                          )
                  )

    # if column is not involved in any addition at certain depth, it is the same as the previous column
    for d in range(1, max_depth+1):
        for col in range(checks.shape[1]):
            s.add(z3.Implies(
                z3.Not(z3.Or([additions[d-1, col_1, col]
                              for col_1 in range(checks.shape[1]) if col != col_1] +
                             [additions[d-1, col, col_1]
                              for col_1 in range(checks.shape[1]) if col != col_1])),
                z3.And([columns[d, i, col] == columns[d-1, i, col]
                        for i in range(checks.shape[0])])))

    # assert that final check matrix has checks.shape[1]-checks.shape[0] zero columns
    s.add(z3.PbEq([(z3.Not(z3.Or([columns[max_depth][i][col]
                                  for i in range(checks.shape[0])])),
                    1) for col in range(checks.shape[1])],
                  checks.shape[1]-checks.shape[0]
                  )
          )

    if s.check() == z3.sat:
        m = s.model()
        additions = [(i, j) for d in range(max_depth)
                     for j in range(checks.shape[1])
                     for i in range(checks.shape[1]) if m[additions[d, i, j]]]

        checks = np.array([[bool(m[columns[max_depth, i, j]]) for j in range(checks.shape[1])]
                           for i in range(checks.shape[0])])

        return additions, checks.astype(int)

    return False


def iterative_search_with_timeout(fun, min_param, max_param, min_timeout, max_timeout, param_factor=2, timeout_factor=2):
    """Geometrically increases the parameter and timeout until a result is found or the maximum timeout is reached.

    Args:
        fun: function to run with increasing parameters and timeouts
        min_param: minimum parameter to start with
        max_param: maximum parameter to reach
        min_timeout: minimum timeout to start with
        max_timeout: maximum timeout to reach
    """
    curr_timeout = min_timeout
    curr_param = min_param
    param_type = type(min_param)
    found = False
    while curr_timeout <= max_timeout:
        while curr_param <= max_param:
            res = _run_with_timeout(fun, curr_param, timeout=curr_timeout)
            if res and res != "timeout":
                return res, curr_param
            curr_param = param_type(curr_param*param_factor)
        curr_timeout *= 2
        curr_param = min_param
    return None, max_param


def depth_optimal_prep_circuit(code: CSSCode, zero_state: bool=True, min_depth=1, max_depth=10, min_timeout=1, max_timeout=3600) -> QuantumCircuit:
    """Synthesize a state preparation circuit for a CSS code that minimizes the depth of the circuit.

    Args:
        code: The CSS code to prepare the state for.
        zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
        starting_depth: The depth of the circuit to start with.
        depth_limit: The maximum depth of the circuit to search for.
        max_cnots: The maximum number of CNOT gates to allow in the circuit. If None, no limit is imposed.
    """
    # first try to find any circuit by exponentially increasing depth
    checks = code.Hx if zero_state else code.Hz
    
    curr_timeout = min_timeout
    curr_depth = min_depth
    circ = None
    
    # while curr_timeout <= max_timeout:
    #     while curr_depth <= max_depth:
    #         res = _run_with_timeout(_generate_circ_with_bounded_depth, checks, curr_depth, timeout=curr_timeout)
    #         if res and res != "timeout":
    #             cnots, reduced_checks = res
    #             circ = _build_circuit_from_list_and_checks(cnots, reduced_checks)
    #             break
    #         curr_depth *= 2
    #     if circ is not None:
    #         break
    #     curr_timeout *= 2
    #     curr_depth = min_depth

    # if circ is None:
    #     return None
    res, curr_depth = iterative_search_with_timeout(lambda depth: _generate_circ_with_bounded_depth(checks, depth), min_depth, max_depth, min_timeout, max_timeout)
    if res:
        cnots, reduced_checks = res
        circ = _build_circuit_from_list_and_checks(cnots, reduced_checks)
    else:
        return None

    # Solving a SAT instance is much faster than proving unsat in this case
    # so we iterate backwards until we find an unsat instance or hit a timeout
    curr_depth -= 1
    while True:
        res = _run_with_timeout(_generate_circ_with_bounded_depth, checks, curr_depth, timeout=max_timeout)
        if res and res != "timeout":
            cnots, reduced_checks = res
            circ = _build_circuit_from_list_and_checks(cnots, reduced_checks)
        else:
            break
        curr_depth -= 1

    return circ
    
