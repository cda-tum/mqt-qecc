"""Synthesizing state preparation circuits for CSS codes."""

from __future__ import annotations

import logging
from collections import deque
from typing import TYPE_CHECKING, Any

import multiprocess
import numpy as np
import z3
from ldpc import mod2
from qiskit import AncillaRegister, ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.converters import circuit_to_dag
from qiskit.dagcircuit import DAGOutNode

logger = logging.getLogger(__name__)

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Callable

    import numpy.typing as npt
    from qiskit import DagCircuit, DAGNode
    from qiskit.quantum_info import PauliList

    from ..code import CSSCode


class StatePrepCircuit:
    """Represents a state preparation circuit for a CSS code."""

    def __init__(self, circ: QuantumCircuit, code: CSSCode, zero_state: bool = True) -> None:
        """Initialize a state preparation circuit.

        Args:
            circ: The state preparation circuit.
            code: The CSS code to prepare the state for.
            zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
        """
        self.circ = circ
        self.code = code
        self.zero_state = zero_state
        self.x_checks = code.Hx.copy() if zero_state else np.vstack((code.Lx.copy(), code.Hx.copy()))
        self.z_checks = code.Hz.copy() if not zero_state else np.vstack((code.Lz.copy(), code.Hz.copy()))
        self.num_qubits = circ.num_qubits
        self.x_fault_sets = [None for _ in range((code.distance - 1) // 2 + 1)]  # type: list[npt.NDArray[np.int8] | None]
        self.z_fault_sets = [None for _ in range((code.distance - 1) // 2 + 1)]  # type: list[npt.NDArray[np.int8] | None]
        self.max_x_measurements = len(self.x_checks)
        self.max_z_measurements = len(self.z_checks)

    def compute_fault_set(
        self, num_errors: int = 1, x_errors: bool = True, reduce: bool = True
    ) -> npt.NDArray[np.int8]:
        """Compute the fault set of the state.

        Args:
            state: The stabilizer state to compute the fault set for.
            num_errors: The number of independent errors to propagate through the circuit.
            x_errors: If True, compute the fault set for X errors. If False, compute the fault set for Z errors.
            reduce: If True, reduce the fault set by the stabilizers of the code to reduce weights.

        Returns:
            The fault set of the state.
        """
        faults = self.x_fault_sets[num_errors] if x_errors else self.z_fault_sets[num_errors]  # type: npt.NDArray[np.int8] | None
        if faults is not None:
            return faults

        if num_errors == 1:
            logging.info("Computing fault set for 1 error.")
            dag = circuit_to_dag(self.circ)
            for node in dag.front_layer():  # remove hadamards
                dag.remove_op_node(node)
            fault_list = []
            # propagate every error before a control
            for node in dag.topological_op_nodes():
                error = _propagate_error(dag, node, x_errors=x_errors)
                fault_list.append(error)
            faults = np.array(fault_list, dtype=np.int8)
            faults = np.unique(faults, axis=0)
        else:
            logging.info(f"Computing fault set for {num_errors} errors.")
            faults = self.compute_fault_set(num_errors - 1, x_errors, reduce=reduce)
            assert faults is not None
            single_faults = self.compute_fault_set(1, x_errors, reduce=reduce)
            new_faults = (faults[:, np.newaxis, :] + single_faults).reshape(-1, self.num_qubits) % 2
            non_propagated_single_errors = np.eye(self.num_qubits, dtype=np.int8)  # type: npt.NDArray[np.int8]
            faults = (new_faults[:, np.newaxis, :] + non_propagated_single_errors).reshape(-1, self.num_qubits) % 2
            # remove duplicates
            faults = np.unique(faults, axis=0)

        # reduce faults by stabilizer
        logging.info("Reducing fault set.")
        reduced = True

        stabs = self.x_checks if x_errors else self.z_checks
        for i, fault in enumerate(faults):
            while reduced:
                reduced = False
                prev_fault = fault
                for stab in stabs:
                    reduced_fault = (prev_fault + stab) % 2
                    if np.sum(prev_fault) > np.sum(reduced_fault):
                        faults[i] = reduced_fault
                        prev_fault = reduced_fault
                        reduced = True
                        break

        # remove trivial faults
        logging.info("Removing trivial faults.")
        for i, fault in enumerate(faults):
            faults[i] = _coset_leader(fault, stabs)
        faults = faults[np.where(np.sum(faults, axis=1) > num_errors)[0]]

        # remove stabilizer equivalent faults
        if reduce:
            logging.info("Removing stabilizer equivalent faults.")
            faults = _remove_stabilizer_equivalent_faults(faults, stabs)
        if x_errors:
            self.x_fault_sets[num_errors] = faults
        else:
            self.z_fault_sets[num_errors] = faults
        return faults


def heuristic_prep_circuit(code: CSSCode, optimize_depth: bool = True, zero_state: bool = True) -> StatePrepCircuit:
    """Return a circuit that prepares the +1 eigenstate of the code w.r.t. the Z or X basis.

    Args:
        code: The CSS code to prepare the state for.
        optimize_depth: If True, optimize the depth of the circuit. This may lead to a higher number of CNOTs.
        zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
    """
    logging.info("Starting heuristic state preparation.")
    checks = code.Hx.copy() if zero_state else code.Hz.copy()
    rank = mod2.rank(checks)

    def is_reduced() -> bool:
        return bool(len(np.where(np.all(checks == 0, axis=0))[0]) == checks.shape[1] - rank)

    costs = np.array([
        [np.sum((checks[:, i] + checks[:, j]) % 2) for j in range(checks.shape[1])] for i in range(checks.shape[1])
    ])
    costs -= np.sum(checks, axis=0)
    np.fill_diagonal(costs, 1)

    used_qubits = []  # type: list[np.int_]
    cnots = []  # type: list[tuple[int, int]]
    while not is_reduced():
        m = np.zeros((checks.shape[1], checks.shape[1]), dtype=bool)  # type: npt.NDArray[np.bool_]
        m[used_qubits, :] = True
        m[:, used_qubits] = True

        costs_unused = np.ma.array(costs, mask=m)  # type: ignore[no-untyped-call]
        if np.all(costs_unused >= 0):  # no more reductions possible
            if used_qubits == []:  # local minimum => get out by making matrix triangular
                logging.warning("Local minimum reached. Making matrix triangular.")
                checks = mod2.reduced_row_echelon(checks)[0]
                costs = np.array([
                    [np.sum((checks[:, i] + checks[:, j]) % 2) for j in range(checks.shape[1])]
                    for i in range(checks.shape[1])
                ])
                costs -= np.sum(checks, axis=0)
                np.fill_diagonal(costs, 1)
            else:  # try to move onto the next layer
                used_qubits = []
            continue

        i, j = np.unravel_index(np.argmin(costs_unused), costs.shape)
        cnots.append((int(i), int(j)))

        if optimize_depth:
            used_qubits.append(i)
            used_qubits.append(j)

        # update checks
        checks[:, j] = (checks[:, i] + checks[:, j]) % 2
        # update costs
        new_weights = np.sum((checks[:, j][:, np.newaxis] + checks) % 2, axis=0)
        costs[j, :] = new_weights - np.sum(checks, axis=0)
        costs[:, j] = new_weights - np.sum(checks[:, j])
        np.fill_diagonal(costs, 1)

    circ = _build_circuit_from_list_and_checks(cnots, checks, zero_state)
    return StatePrepCircuit(circ, code, zero_state)


def _generate_circ_with_bounded_depth(
    checks: npt.NDArray[np.int8], max_depth: int, zero_state: bool = True
) -> QuantumCircuit | None:
    assert max_depth > 0, "max_depth should be greater than 0"
    columns = np.array([
        [[z3.Bool(f"x_{d}_{i}_{j}") for j in range(checks.shape[1])] for i in range(checks.shape[0])]
        for d in range(max_depth + 1)
    ])

    additions = np.array([
        [[z3.Bool(f"add_{d}_{i}_{j}") for j in range(checks.shape[1])] for i in range(checks.shape[1])]
        for d in range(max_depth)
    ])
    n_cols = checks.shape[1]
    s = z3.Solver()

    # create initial matrix
    columns[0, :, :] = checks.astype(bool)

    s.add(_column_addition_contraint(columns, additions))

    # qubit can be involved in at most one addition at each depth
    for d in range(max_depth):
        for col in range(n_cols):
            s.add(
                z3.PbLe(
                    [(additions[d, col_1, col], 1) for col_1 in range(n_cols) if col != col_1]
                    + [(additions[d, col, col_2], 1) for col_2 in range(n_cols) if col != col_2],
                    1,
                )
            )

    # if column is not involved in any addition at certain depth, it is the same as the previous column
    for d in range(1, max_depth + 1):
        for col in range(n_cols):
            s.add(
                z3.Implies(
                    z3.Not(
                        z3.Or(
                            list(np.delete(additions[d - 1, :, col], [col]))
                            + list(np.delete(additions[d - 1, col, :], [col]))
                        )
                    ),
                    _symbolic_vector_eq(columns[d, :, col], columns[d - 1, :, col]),
                )
            )

    s.add(_final_matrix_constraint(columns))

    if s.check() == z3.sat:
        m = s.model()
        cnots = [
            (i, j)
            for d in range(max_depth)
            for j in range(checks.shape[1])
            for i in range(checks.shape[1])
            if m[additions[d, i, j]]
        ]

        checks = np.array([
            [bool(m[columns[max_depth, i, j]]) for j in range(checks.shape[1])] for i in range(checks.shape[0])
        ])

        return _build_circuit_from_list_and_checks(cnots, checks, zero_state=zero_state)

    return None


def _generate_circ_with_bounded_gates(
    checks: npt.NDArray[np.int8], max_cnots: int, zero_state: bool = True
) -> QuantumCircuit:
    """Find the gate optimal circuit for a given check matrix and maximum depth."""
    columns = np.array([
        [[z3.Bool(f"x_{d}_{i}_{j}") for j in range(checks.shape[1])] for i in range(checks.shape[0])]
        for d in range(max_cnots + 1)
    ])
    n_bits = int(np.ceil(np.log2(checks.shape[1])))
    targets = [z3.BitVec(f"target_{d}", n_bits) for d in range(max_cnots)]
    controls = [z3.BitVec(f"control_{d}", n_bits) for d in range(max_cnots)]
    s = z3.Solver()

    additions = np.array([
        [
            [z3.And(controls[d] == col_1, targets[d] == col_2) for col_2 in range(checks.shape[1])]
            for col_1 in range(checks.shape[1])
        ]
        for d in range(max_cnots)
    ])

    # create initial matrix
    columns[0, :, :] = checks.astype(bool)
    s.add(_column_addition_contraint(columns, additions))

    for d in range(1, max_cnots + 1):
        # qubit cannot be control and target at the same time
        s.add(controls[d - 1] != targets[d - 1])

        # control and target must be valid qubits
        if checks.shape[1] and (checks.shape[1] - 1) != 0:
            s.add(z3.ULT(controls[d - 1], checks.shape[1]))
            s.add(z3.ULT(targets[d - 1], checks.shape[1]))

    # if column is not involved in any addition at certain depth, it is the same as the previous column
    for d in range(1, max_cnots + 1):
        for col in range(checks.shape[1]):
            s.add(z3.Implies(targets[d - 1] != col, _symbolic_vector_eq(columns[d, :, col], columns[d - 1, :, col])))

    # assert that final check matrix has checks.shape[1]-checks.shape[0] zero columns
    s.add(_final_matrix_constraint(columns))

    if s.check() == z3.sat:
        m = s.model()
        cnots = [(m[controls[d]].as_long(), m[targets[d]].as_long()) for d in range(max_cnots)]
        checks = np.array([
            [bool(m[columns[max_cnots][i][j]]) for j in range(checks.shape[1])] for i in range(checks.shape[0])
        ]).astype(int)
        return _build_circuit_from_list_and_checks(cnots, checks, zero_state=zero_state)

    return None


def _optimal_circuit(
    code: CSSCode,
    prep_func: Callable[[npt.NDArray[np.int8], int, bool], QuantumCircuit | None],
    zero_state: bool = True,
    min_param: int = 1,
    max_param: int = 10,
    min_timeout: int = 1,
    max_timeout: int = 3600,
) -> StatePrepCircuit | None:
    """Synthesize a state preparation circuit for a CSS code that minimizes the circuit w.r.t. some metric param according to prep_func.

    Args:
        code: The CSS code to prepare the state for.
        zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
        prep_func: The function to optimize the circuit with respect to.
        min_param: minimum parameter to start with
        max_param: maximum parameter to reach
        min_timeout: minimum timeout to start with
        max_timeout: maximum timeout to reach
    """
    checks = code.Hx if zero_state else code.Hz

    def fun(param: int) -> QuantumCircuit | None:
        return prep_func(checks, param, zero_state)

    res = iterative_search_with_timeout(
        fun,
        min_param,
        max_param,
        min_timeout,
        max_timeout,
    )

    if res is None:
        return None
    circ, curr_param = res
    if circ is None:
        return None

    logging.info(f"Solution found with param {curr_param}")
    # Solving a SAT instance is much faster than proving unsat in this case
    # so we iterate backwards until we find an unsat instance or hit a timeout
    logging.info("Trying to minimize param")
    while True:
        logging.info(f"Trying param {curr_param - 1}")
        opt_res = _run_with_timeout(prep_func, checks, curr_param - 1, timeout=max_timeout)  # type: ignore[arg-type]
        if opt_res is None or (isinstance(opt_res, str) and opt_res == "timeout"):
            break
        circ = opt_res
        curr_param -= 1

    logging.info(f"Optimal param: {curr_param}")
    return StatePrepCircuit(circ, code, zero_state)


def depth_optimal_prep_circuit(
    code: CSSCode,
    zero_state: bool = True,
    min_depth: int = 1,
    max_depth: int = 10,
    min_timeout: int = 1,
    max_timeout: int = 3600,
) -> StatePrepCircuit | None:
    """Synthesize a state preparation circuit for a CSS code that minimizes the circuit depth.

    Args:
        code: The CSS code to prepare the state for.
        zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
        min_depth: minimum depth to start with
        max_depth: maximum depth to reach
        min_timeout: minimum timeout to start with
        max_timeout: maximum timeout to reach
    """
    return _optimal_circuit(
        code, _generate_circ_with_bounded_depth, zero_state, min_depth, max_depth, min_timeout, max_timeout
    )


def gate_optimal_prep_circuit(
    code: CSSCode,
    zero_state: bool = True,
    min_gates: int = 1,
    max_gates: int = 10,
    min_timeout: int = 1,
    max_timeout: int = 3600,
) -> StatePrepCircuit | None:
    """Synthesize a state preparation circuit for a CSS code that minimizes the number of gates.

    Args:
        code: The CSS code to prepare the state for.
        zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
        min_gates: minimum number of gates to start with
        max_gates: maximum number of gates to reach
        min_timeout: minimum timeout to start with
        max_timeout: maximum timeout to reach
    """
    return _optimal_circuit(
        code, _generate_circ_with_bounded_gates, zero_state, min_gates, max_gates, min_timeout, max_timeout
    )


def _build_circuit_from_list_and_checks(
    cnots: list[tuple[int, int]], checks: npt.NDArray[np.int8], zero_state: bool = True
) -> QuantumCircuit:
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
        if zero_state:
            ctrl, tar = i, j
        else:
            ctrl, tar = j, i
        circ.cx(ctrl, tar)
    return circ


def _run_with_timeout(func: Callable[[Any], Any], *args: Any, timeout: int = 10) -> Any | str | None:  # noqa: ANN401
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
            res = _run_with_timeout(fun, curr_param, timeout=curr_timeout)
            if res is not None and (not isinstance(res, str) or res != "timeout"):
                return res, curr_param
            if curr_param == max_param:
                break

            curr_param = int(curr_param * param_factor)
            curr_param = min(curr_param, max_param)

        curr_timeout = int(curr_timeout * timeout_factor)
        curr_param = min_param
    return None, max_param


def gate_optimal_verification_stabilizers(
    sp_circ: StatePrepCircuit,
    x_errors: bool = True,
    min_timeout: int = 1,
    max_timeout: int = 3600,
    max_ancillas: int | None = None,
) -> list[list[npt.NDArray[np.int8]]]:
    """Return verification stabilizers for the state preparation circuit.

    The method uses an iterative search to find the optimal set of stabilizers by repeatedly computing the optimal circuit for each number of ancillas and cnots. This is repeated for each number of independent correctable errors in the state preparation circuit. Thus the verification circuit is constructed of multiple "layers" of stabilizers, each layer corresponding to a fault set it verifies.

    Args:
        sp_circ: The state preparation circuit to verify.
        x_errors: If True, verify the X errors. If False, verify the Z errors.
        min_timeout: The minimum time to allow each search to run for.
        max_timeout: The maximum time to allow each search to run for.
        max_ancillas: The maximum number of ancillas to allow in each layer verification circuit.

    Returns:
        A list of stabilizers to verify the state preparation circuit.
    """
    max_errors = (sp_circ.code.distance - 1) // 2
    layers = [[] for _ in range(max_errors)]  # type: list[list[npt.NDArray[np.int8]]]
    if max_ancillas is None:
        max_ancillas = sp_circ.max_z_measurements if x_errors else sp_circ.max_x_measurements
    # Find the optimal circuit for every number of errors in the preparation circuit
    for num_errors in range(1, max_errors + 1):
        logging.info(f"Finding verification stabilizers for {num_errors} errors")
        faults = sp_circ.compute_fault_set(num_errors, x_errors)
        if len(faults) == 0:
            logging.info(f"No non-trivial faults for {num_errors} errors")
            layers[num_errors - 1] = []
            continue
        # Start with maximal number of ancillas
        # Minimal CNOT solution must be achievable with these
        num_anc = max_ancillas
        checks = sp_circ.z_checks if x_errors else sp_circ.x_checks
        min_cnots = np.min(np.sum(checks, axis=1))
        max_cnots = np.sum(checks)

        logging.info(
            f"Finding verification stabilizers for {num_errors} errors with {min_cnots} to {max_cnots} CNOTs using {num_anc} ancillas"
        )

        def fun(num_cnots: int) -> list[npt.NDArray[np.int8]] | None:
            return verification_stabilizers(sp_circ, num_anc, num_cnots, num_errors, x_errors=x_errors)  # noqa: B023

        res = iterative_search_with_timeout(
            fun,
            min_cnots,
            max_cnots,
            min_timeout,
            max_timeout,
        )

        if res is None:
            logging.info(f"No verification stabilizers found for {num_errors} errors")
            layers[num_errors - 1] = []
            continue
        measurements, num_cnots = res
        if measurements is None or (isinstance(measurements, str) and measurements == "timeout"):
            logging.info(f"No verification stabilizers found for {num_errors} errors")
            return []  # No solution found

        logging.info(f"Found verification stabilizers for {num_errors} errors with {num_cnots} CNOTs")
        # If any measurements are unused we can reduce the number of ancillas at least by that
        num_anc = np.sum([np.any(m) for m in measurements])
        measurements = [m for m in measurements if np.any(m)]

        # Iterate backwards to find the minimal number of cnots
        logging.info(f"Finding minimal number of CNOTs for {num_errors} errors")

        def search_cnots(num_cnots: int) -> list[npt.NDArray[np.int8]] | None:
            return verification_stabilizers(sp_circ, num_anc, num_cnots, num_errors, x_errors=x_errors)  # noqa: B023

        while num_cnots - 1 > 0:
            logging.info(f"Trying {num_cnots - 1} CNOTs")

            cnot_opt = _run_with_timeout(
                search_cnots,
                num_cnots - 1,
                timeout=max_timeout,
            )
            if cnot_opt is None or (isinstance(cnot_opt, str) and cnot_opt == "timeout"):
                break
            num_cnots -= 1
            measurements = cnot_opt
        logging.info(f"Minimal number of CNOTs for {num_errors} errors is: {num_cnots}")

        # If the number of CNOTs is minimal, we can reduce the number of ancillas
        logging.info(f"Finding minimal number of ancillas for {num_errors} errors")
        while num_anc - 1 > 0:
            logging.info(f"Trying {num_anc - 1} ancillas")

            def search_anc(num_anc: int) -> list[npt.NDArray[np.int8]] | None:
                return verification_stabilizers(sp_circ, num_anc, num_cnots, num_errors, x_errors=x_errors)  # noqa: B023

            anc_opt = _run_with_timeout(
                search_anc,
                num_cnots - 1,
                timeout=max_timeout,
            )
            if anc_opt is None or (isinstance(anc_opt, str) and anc_opt == "timeout"):
                break
            num_anc -= 1
            measurements = anc_opt
        logging.info(f"Minimal number of ancillas for {num_errors} errors is: {num_anc}")
        layers[num_errors - 1] = measurements

    return layers


def gate_optimal_verification_circuit(
    sp_circ: StatePrepCircuit, min_timeout: int = 1, max_timeout: int = 3600, max_ancillas: int | None = None
) -> QuantumCircuit:
    """Return a verified state preparation circuit.

    The verification circuit is a set of stabilizers such that each propagated error in sp_circ anticommutes with some verification stabilizer.

    The method uses an iterative search to find the optimal set of stabilizers by repeatedly computing the optimal circuit for each number of ancillas and cnots. This is repeated for each number of independent correctable errors in the state preparation circuit. Thus the verification circuit is constructed of multiple "layers" of stabilizers, each layer corresponding to a fault set it verifies.

    Args:
        sp_circ: The state preparation circuit to verify.
        min_timeout: The minimum time to allow each search to run for.
        max_timeout: The maximum time to allow each search to run for.
        max_ancillas: The maximum number of ancillas to allow in each layer verification circuit.
    """
    logging.info("Finding optimal verification stabilizers for X errors")
    x_layers = gate_optimal_verification_stabilizers(sp_circ, True, min_timeout, max_timeout, max_ancillas)

    z_layers = gate_optimal_verification_stabilizers(sp_circ, False, min_timeout, max_timeout, max_ancillas)

    z_measurements = [measurement for layer in x_layers for measurement in layer]
    x_measurements = [measurement for layer in z_layers for measurement in layer]
    return _measure_stabs(sp_circ.circ, x_measurements, z_measurements)


def heuristic_verification_circuit(
    sp_circ: StatePrepCircuit, max_covering_sets: int = 10000, find_coset_leaders: bool = True
) -> QuantumCircuit:
    """Return a verified state preparation circuit.

    The method uses a greedy set covering heuristic to find a small set of stabilizers that verifies the state preparation circuit. The heuristic is not guaranteed to find the optimal set of stabilizers.

    Args:
        sp_circ: The state preparation circuit to verify.
        max_covering_sets: The maximum number of covering sets to consider.
        find_coset_leaders: Whether to find coset leaders for the found measurements. This is done using SAT solvers so it can be slow.
    """
    x_layers = heuristic_verification_stabilizers(sp_circ, True, max_covering_sets, find_coset_leaders)
    z_layers = heuristic_verification_stabilizers(sp_circ, False, max_covering_sets, find_coset_leaders)

    z_measurements = [measurement for layer in x_layers for measurement in layer]
    x_measurements = [measurement for layer in z_layers for measurement in layer]
    return _measure_stabs(sp_circ.circ, x_measurements, z_measurements)


def heuristic_verification_stabilizers(
    sp_circ: StatePrepCircuit, x_errors: bool = True, max_covering_sets: int = 10000, find_coset_leaders: bool = True
) -> list[list[npt.NDArray[np.int8]]]:
    """Return verification stabilizers for the preparation circuit.

    Args:
        sp_circ: The state preparation circuit to verify.
        x_errors: Whether to find verification stabilizers for X errors. If False, find for Z errors.
        max_covering_sets: The maximum number of covering sets to consider.
        find_coset_leaders: Whether to find coset leaders for the found measurements. This is done using SAT solvers so it can be slow.
    """
    logging.info("Finding verification stabilizers using heuristic method")
    max_errors = (sp_circ.code.distance - 1) // 2
    layers = [[] for _ in range(max_errors)]  # type: list[list[npt.NDArray[np.int8]]]
    for num_errors in range(1, max_errors + 1):
        logging.info(f"Finding verification stabilizers for {num_errors} errors")
        faults = sp_circ.compute_fault_set(num_errors, x_errors)
        logging.info(f"There are {len(faults)} faults")
        if len(faults) == 0:
            layers[num_errors - 1] = []
            continue

        orthogonal_checks = sp_circ.z_checks if x_errors else sp_circ.x_checks
        syndromes = orthogonal_checks @ faults.T % 2
        candidates = np.where(np.any(syndromes != 0, axis=1))[0]
        non_candidates = np.where(np.all(syndromes == 0, axis=1))[0]
        candidate_checks = orthogonal_checks[candidates]
        non_candidate_checks = orthogonal_checks[non_candidates]

        def covers(s: npt.NDArray[np.int8], faults: npt.NDArray[np.int8]) -> frozenset[int]:
            return frozenset(np.where(s @ faults.T % 2 != 0)[0])

        logging.info("Converting Stabilizer Checks to covering sets")
        candidate_sets_ordered = [(covers(s, faults), s, i) for i, s in enumerate(candidate_checks)]
        candidate_sets_ordered.sort(key=lambda x: -np.sum(x[1]))
        mapping = {
            cand: _coset_leader(candidate_checks[i], non_candidate_checks) for cand, _, i in candidate_sets_ordered
        }
        candidate_sets = {cand for cand, _, _ in candidate_sets_ordered}

        def set_cover(
            n: int, cands: set[frozenset[int]], mapping: dict[frozenset[int], npt.NDArray[np.int8]]
        ) -> list[frozenset[int]]:
            universe = set(range(n))
            cover = []

            def sort_key(stab: frozenset[int], universe: set[int] = universe) -> tuple[int, np.int_]:
                return (len(stab & universe), -np.sum(mapping[stab]))  # noqa: B023

            while universe:
                best = max(cands, key=sort_key)
                cover.append(best)
                universe -= best
            return cover

        improved = True
        logging.info("Finding initial set cover")
        cover = set_cover(len(faults), candidate_sets, mapping)
        logging.info(f"Initial set cover has {len(cover)} sets")
        cost1 = len(cover)
        cost2 = sum(np.sum(mapping[stab]) for stab in cover)
        prev_candidates = candidate_sets.copy()
        while improved and len(candidate_sets) < max_covering_sets:
            improved = False
            # add all symmetric differences to candidates
            to_remove = set()  # type: set[frozenset[int]]
            to_add = set()  # type: set[frozenset[int]]
            for c1 in candidate_sets:
                for c2 in candidate_sets:
                    if len(to_add) >= max_covering_sets:
                        break
                    comb = c1 ^ c2
                    if c1 == c2 or comb in candidate_sets or comb in to_add or comb in to_remove:
                        continue

                    mapping[comb] = (mapping[c1] + mapping[c2]) % 2
                    if len(c1 & c2) == 0:
                        to_remove.add(c1)
                        to_remove.add(c2)
                    to_add.add(comb)
            candidate_sets = candidate_sets.union(to_add)
            new_cover = set_cover(len(faults), candidate_sets, mapping)
            logging.info(f"New Covering set has {len(new_cover)} sets")
            new_cost1 = len(new_cover)
            new_cost2 = sum(np.sum(mapping[stab]) for stab in new_cover)
            if new_cost1 < cost1 or (new_cost1 == cost1 and new_cost2 < cost2):
                cover = new_cover
                cost1 = new_cost1
                cost2 = new_cost2
                improved = True
            elif candidate_sets == prev_candidates:
                break
            prev_candidates = candidate_sets

        logging.info(f"Found covering set of size {len(cover)} for {num_errors} errors")
        measurements = [mapping[c] for c in cover]
        if find_coset_leaders and len(non_candidates) > 0:
            logging.info(f"Finding coset leaders for {num_errors} errors")
            measurements = [_coset_leader(m, non_candidate_checks) for m in measurements]
        logging.info(f"Found {np.sum(measurements)} CNOTS for {num_errors} errors")
        layers[num_errors - 1] = measurements

    return layers


def _measure_stabs(
    circ: QuantumCircuit, x_measurements: list[npt.NDArray[np.int8]], z_measurements: list[npt.NDArray[np.int8]]
) -> QuantumCircuit:
    # Create the verification circuit
    num_z_anc = len(z_measurements)
    num_x_anc = len(x_measurements)
    q = QuantumRegister(circ.num_qubits, "q")
    z_anc = AncillaRegister(num_z_anc, "z_anc")
    x_anc = AncillaRegister(num_x_anc, "x_anc")
    z_c = ClassicalRegister(num_z_anc, "z_c")
    x_c = ClassicalRegister(num_x_anc, "x_c")
    measured_circ = QuantumCircuit(q, x_anc, z_anc, x_c, z_c)

    measured_circ.compose(circ, inplace=True)
    current_anc = 0
    for measurement in z_measurements:
        for qubit in np.where(measurement == 1)[0]:
            measured_circ.cx(q[qubit], z_anc[current_anc])
        measured_circ.measure(z_anc[current_anc], z_c[current_anc])

    current_anc = 0
    for measurement in x_measurements:
        measured_circ.h(x_anc[current_anc])
        for qubit in np.where(measurement == 1)[0]:
            measured_circ.cx(x_anc[current_anc], q[qubit])
        measured_circ.h(x_anc[current_anc])
        measured_circ.measure(x_anc[current_anc], x_c[current_anc])
        current_anc += 1

    return measured_circ


def _vars_to_stab(
    measurement: list[z3.BoolRef | bool], generators: npt.NDArray[np.int8]
) -> npt.NDArray[z3.BoolRef | bool]:  # type: ignore[type-var]
    measurement_stab = _symbolic_scalar_mult(generators[0], measurement[0])
    for i, scalar in enumerate(measurement[1:]):
        measurement_stab = _symbolic_vector_add(measurement_stab, _symbolic_scalar_mult(generators[i + 1], scalar))
    return measurement_stab


def verification_stabilizers(
    sp_circ: StatePrepCircuit, num_anc: int, num_cnots: int, num_errors: int, x_errors: bool = True
) -> list[npt.NDArray[np.int8]] | None:
    """Return verification stabilizers for num_errors independent errors in the state preparation circuit using z3.

    Args:
        sp_circ: The state preparation circuit.
        num_anc: The maximum number of ancilla qubits to use.
        num_cnots: The maximumg number of CNOT gates to use.
        num_errors: The number of errors occur in the state prep circuit.
        x_errors: If True, the errors are X errors. Otherwise, the errors are Z errors.
    """
    # Measurements are written as sums of generators
    # The variables indicate which generators are non-zero in the sum
    gens = sp_circ.z_checks if x_errors else sp_circ.x_checks
    n_gens = gens.shape[0]

    measurement_vars = [[z3.Bool(f"m_{anc}_{i}") for i in range(n_gens)] for anc in range(num_anc)]
    solver = z3.Solver()

    measurement_stabs = [_vars_to_stab(vars_, gens) for vars_ in measurement_vars]

    # assert that each error is detected
    errors = sp_circ.compute_fault_set(num_errors, x_errors)
    solver.add(
        z3.And([
            z3.PbGe([(_odd_overlap(measurement, error), 1) for measurement in measurement_stabs], 1) for error in errors
        ])
    )

    # assert that not too many CNOTs are used
    solver.add(
        z3.PbLe(
            [(measurement[q], 1) for measurement in measurement_stabs for q in range(sp_circ.num_qubits)], num_cnots
        )
    )

    if solver.check() == z3.sat:
        model = solver.model()
        # Extract stabilizer measurements from model
        actual_measurements = []
        for m in measurement_vars:
            v = np.zeros(sp_circ.num_qubits, dtype=np.int8)  # type: npt.NDArray[np.int8]
            for g in range(n_gens):
                if model[m[g]]:
                    v += gens[g]
            actual_measurements.append(v % 2)

        return actual_measurements
    return None


def _coset_leader(error: npt.NDArray[np.int8], generators: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
    if len(generators) == 0:
        return error
    s = z3.Optimize()
    leader = [z3.Bool(f"e_{i}") for i in range(len(error))]
    coeff = [z3.Bool(f"c_{i}") for i in range(len(generators))]

    g = _vars_to_stab(coeff, generators)

    s.add(_symbolic_vector_eq(np.array(leader), _symbolic_vector_add(error.astype(bool), g)))
    s.minimize(z3.Sum(leader))

    s.check()  # always SAT
    m = s.model()
    return np.array([bool(m[leader[i]]) for i in range(len(error))]).astype(int)


def _symbolic_scalar_mult(v: npt.NDArray[np.int8], a: z3.BoolRef | bool) -> npt.NDArray[z3.BoolRef]:
    """Multiply a concrete vector by a symbolic scalar."""
    return np.array([a if s == 1 else False for s in v])


def _symbolic_vector_add(
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


def _odd_overlap(v_sym: npt.NDArray[z3.BoolRef | bool], v_con: npt.NDArray[np.int8]) -> z3.BoolRef:
    """Return True if the overlap of symbolic vector with constant vector is odd."""
    return z3.PbEq([(v_sym[i], 1) for i, c in enumerate(v_con) if c == 1], 1)


def _symbolic_vector_eq(v1: npt.NDArray[z3.BoolRef | bool], v2: npt.NDArray[z3.BoolRef | bool]) -> z3.BoolRef:
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


def _column_addition_contraint(
    columns: npt.NDArray[z3.BoolRef | bool],
    col_add_vars: npt.NDArray[z3.BoolRef],
) -> z3.BoolRef:
    assert len(columns.shape) == 3
    max_depth = col_add_vars.shape[0]  # type: ignore[unreachable]
    n_cols = col_add_vars.shape[2]

    constraints = []
    for d in range(1, max_depth + 1):
        for col_1 in range(n_cols):
            for col_2 in range(col_1 + 1, n_cols):
                col_sum = _symbolic_vector_add(columns[d - 1, :, col_1], columns[d - 1, :, col_2])

                # encode col_2 += col_1
                add_col1_to_col2 = z3.Implies(
                    col_add_vars[d - 1, col_1, col_2],
                    z3.And(
                        _symbolic_vector_eq(columns[d, :, col_2], col_sum),
                        _symbolic_vector_eq(columns[d, :, col_1], columns[d - 1, :, col_1]),
                    ),
                )

                # encode col_1 += col_2
                add_col2_to_col1 = z3.Implies(
                    col_add_vars[d - 1, col_2, col_1],
                    z3.And(
                        _symbolic_vector_eq(columns[d, :, col_1], col_sum),
                        _symbolic_vector_eq(columns[d, :, col_2], columns[d - 1, :, col_2]),
                    ),
                )

                constraints.extend([add_col1_to_col2, add_col2_to_col1])

    return z3.And(constraints)


def _final_matrix_constraint(columns: npt.NDArray[z3.BoolRef | bool]) -> z3.BoolRef:
    assert len(columns.shape) == 3
    return z3.PbEq(  # type: ignore[unreachable]
        [(z3.Not(z3.Or(list(columns[-1, :, col]))), 1) for col in range(columns.shape[2])],
        columns.shape[2] - columns.shape[1],
    )


def _propagate_error(dag: DagCircuit, node: DAGNode, x_errors: bool = True) -> PauliList:
    """Propagates a Pauli error through a circuit beginning from control of node."""
    control = node.qargs[0]._index  # noqa: SLF001
    error = np.array([0] * dag.num_qubits(), dtype=np.int8)  # type: npt.NDArray[np.int8]
    error[control] = 1
    # propagate error through circuit via bfs
    q = deque([node])
    visited = set()  # type: set[DAGNode]
    while q:
        node = q.popleft()
        if node in visited or isinstance(node, DAGOutNode):
            continue
        control = node.qargs[0]._index  # noqa: SLF001
        target = node.qargs[1]._index  # noqa: SLF001
        if x_errors:
            error[target] = (error[target] + error[control]) % 2
        else:
            error[control] = (error[target] + error[control]) % 2
        for succ in dag.successors(node):
            q.append(succ)
    return error


def _remove_stabilizer_equivalent_faults(
    faults: npt.NDArray[np.int8], stabilizers: npt.NDArray[np.int8]
) -> npt.NDArray[np.int8]:
    """Remove stabilizer equivalent faults from a list of faults."""
    faults = faults.copy()
    stabilizers = stabilizers.copy()
    removed = set()

    logging.debug(f"Removing stabilizer equivalent faults from {len(faults)} faults.")
    for i, f1 in enumerate(faults):
        if i in removed:
            continue
        stabs_ext1 = np.vstack((stabilizers, f1))
        if mod2.rank(stabs_ext1) == mod2.rank(stabilizers):
            removed.add(i)
            continue

        for j, f2 in enumerate(faults[i + 1 :]):
            if j + i + 1 in removed:
                continue
            stabs_ext2 = np.vstack((stabs_ext1, f2))

            if mod2.rank(stabs_ext2) == mod2.rank(stabs_ext1):
                removed.add(j + i + 1)

    logging.debug(f"Removed {len(removed)} stabilizer equivalent faults.")
    indices = list(set(range(len(faults))) - removed)
    if len(indices) == 0:
        return np.array([])

    return faults[indices]
