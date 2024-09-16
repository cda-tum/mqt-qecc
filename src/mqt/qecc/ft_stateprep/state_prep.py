"""Synthesizing state preparation circuits for CSS codes."""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import TYPE_CHECKING, Any

import multiprocess
import numpy as np
import z3
from ldpc import mod2
from qiskit import AncillaRegister, ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.converters import circuit_to_dag

from ..codes import InvalidCSSCodeError

logger = logging.getLogger(__name__)

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Callable

    import numpy.typing as npt
    from qiskit import AncillaQubit, ClBit, DAGNode, Qubit
    from qiskit.quantum_info import PauliList

    from ..codes import CSSCode


class StatePrepCircuit:
    """Represents a state preparation circuit for a CSS code."""

    def __init__(
        self, circ: QuantumCircuit, code: CSSCode, zero_state: bool = True, error_detection_code: bool = False
    ) -> None:
        """Initialize a state preparation circuit.

        Args:
            circ: The state preparation circuit.
            code: The CSS code to prepare the state for.
            zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
            error_detection_code: If True, prepare the state for error detection. This ensures that when computing the fault set of the circuit, up to d//2 errors errors can occur in the circuit.
        """
        self.circ = circ
        self.code = code
        self.zero_state = zero_state

        if code.Hx is None or code.Hz is None:
            msg = "The CSS code must have both X and Z checks."
            raise InvalidCSSCodeError(msg)

        self.x_checks = code.Hx.copy() if zero_state else np.vstack((code.Lx.copy(), code.Hx.copy()))
        self.z_checks = code.Hz.copy() if not zero_state else np.vstack((code.Lz.copy(), code.Hz.copy()))

        self.num_qubits = circ.num_qubits

        self.error_detection_code = error_detection_code
        self._set_max_errors()

        self.max_x_measurements = len(self.x_checks)
        self.max_z_measurements = len(self.z_checks)

    def set_error_detection(self, error_detection: bool) -> None:
        """Set whether the state preparation circuit is for error detection."""
        self.error_detection_code = error_detection
        self._set_max_errors()

    def compute_fault_sets(self, reduce: bool = True) -> None:
        """Compute the fault sets for the state preparation circuit."""
        self.compute_fault_set(self.max_errors, x_errors=True, reduce=reduce)
        self.compute_fault_set(self.max_errors, x_errors=False, reduce=reduce)

    def compute_fault_set(
        self, num_errors: int = 1, x_errors: bool = True, reduce: bool = True
    ) -> npt.NDArray[np.int8]:
        """Compute the fault set of the state.

        Args:
            num_errors: The number of independent errors to propagate through the circuit.
            x_errors: If True, compute the fault set for X errors. If False, compute the fault set for Z errors.
            reduce: If True, reduce the fault set by the stabilizers of the code to reduce weights.

        Returns:
            The fault set of the state.
        """
        faults: npt.NDArray[np.int8] | None = (
            self.x_fault_sets[num_errors] if x_errors else self.z_fault_sets[num_errors]
        )
        if faults is not None:
            return faults

        if num_errors == 1:
            logging.info("Computing fault set for 1 error.")
            dag = circuit_to_dag(self.circ)
            for node in dag.front_layer():  # remove hadamards
                dag.remove_op_node(node)
            fault_list = []
            # propagate every error before a control
            nodes = list(dag.topological_op_nodes())
            for i in range(len(nodes)):
                error = _propagate_error(nodes[i:], dag.num_qubits(), x_errors=x_errors)
                fault_list.append(error)
            faults = np.array(fault_list, dtype=np.int8)
            faults = np.unique(faults, axis=0)

            if x_errors and self.x_fault_sets_unreduced[1] is None:
                non_propagated_single_errors = np.eye(self.num_qubits, dtype=np.int8)
                self.x_fault_sets_unreduced[1] = np.vstack((faults, non_propagated_single_errors))
            elif not x_errors and self.z_fault_sets[1] is None:
                non_propagated_single_errors = np.eye(self.num_qubits, dtype=np.int8)
                self.z_fault_sets_unreduced[1] = np.vstack((faults, non_propagated_single_errors))
        else:
            logging.info(f"Computing fault set for {num_errors} errors.")
            self.compute_fault_set(num_errors - 1, x_errors, reduce=reduce)
            if x_errors:
                faults = self.x_fault_sets_unreduced[num_errors - 1]
                single_faults = self.x_fault_sets_unreduced[1]
            else:
                faults = self.z_fault_sets_unreduced[num_errors - 1]
                single_faults = self.z_fault_sets_unreduced[1]

            assert faults is not None
            assert single_faults is not None

            new_faults = (faults[:, np.newaxis, :] + single_faults).reshape(-1, self.num_qubits) % 2
            # remove duplicates
            faults = np.unique(new_faults, axis=0)
            if x_errors:
                self.x_fault_sets_unreduced[num_errors] = faults.copy()
            else:
                self.z_fault_sets_unreduced[num_errors] = faults.copy()

        # reduce faults by stabilizer
        stabs = self.x_checks if x_errors else self.z_checks
        faults = _remove_trivial_faults(faults, stabs, self.code, x_errors, num_errors)

        # remove stabilizer equivalent faults
        if reduce:
            logging.info("Removing stabilizer equivalent faults.")
            faults = _remove_stabilizer_equivalent_faults(faults, stabs)
        if x_errors:
            self.x_fault_sets[num_errors] = faults
        else:
            self.z_fault_sets[num_errors] = faults
        return faults

    def combine_faults(
        self, additional_faults: npt.NDArray[np.int8], x_errors: bool = True
    ) -> list[npt.NDArray[np.int8] | None]:
        """Combine fault sets of circuit with additional independent faults.

        Args:
            additional_faults: The additional faults to combine with the fault set of the circuit.
            x_errors: If True, combine the fault sets for X errors. If False, combine the fault sets for Z errors.
        """
        self.compute_fault_sets()
        if len(additional_faults) == 0:
            return self.x_fault_sets if x_errors else self.z_fault_sets

        fault_sets_unreduced = self.x_fault_sets_unreduced.copy() if x_errors else self.z_fault_sets_unreduced.copy()
        assert fault_sets_unreduced[1] is not None
        fault_sets_unreduced[1] = np.vstack((fault_sets_unreduced[1], additional_faults))

        for i in range(1, self.max_errors):
            uncombined = fault_sets_unreduced[i]
            assert uncombined is not None
            combined = (uncombined[:, np.newaxis, :] + additional_faults).reshape(-1, self.num_qubits) % 2
            next_faults = fault_sets_unreduced[i + 1]
            assert next_faults is not None
            fault_sets_unreduced[i + 1] = np.vstack((next_faults, combined))
        fault_sets: list[npt.NDArray[np.int8] | None] = [None for _ in range(self.max_errors + 1)]
        stabs = self.x_checks if x_errors else self.z_checks
        for num_errors in range(1, self.max_errors + 1):
            fs = fault_sets_unreduced[num_errors]
            assert fs is not None
            fault_sets[num_errors] = _remove_trivial_faults(fs, stabs, self.code, x_errors, num_errors)
        return fault_sets

    def _set_max_errors(self) -> None:
        if self.code.distance == 2:
            logging.warning("Code distance is 2, assuming error detection code.")
            self.error_detection_code = True

        self.max_errors = (self.code.distance - 1) // 2 if not self.error_detection_code else self.code.distance // 2
        self.max_x_errors = (
            (self.code.x_distance - 1) // 2 if not self.error_detection_code else self.code.x_distance // 2
        )
        self.max_z_errors = (
            (self.code.z_distance - 1) // 2 if not self.error_detection_code else self.code.z_distance // 2
        )
        self.x_fault_sets: list[npt.NDArray[np.int8] | None] = [None for _ in range(self.max_errors + 1)]
        self.z_fault_sets: list[npt.NDArray[np.int8] | None] = [None for _ in range(self.max_errors + 1)]
        self.x_fault_sets_unreduced: list[npt.NDArray[np.int8] | None] = [None for _ in range(self.max_errors + 1)]
        self.z_fault_sets_unreduced: list[npt.NDArray[np.int8] | None] = [None for _ in range(self.max_errors + 1)]


def heuristic_prep_circuit(code: CSSCode, optimize_depth: bool = True, zero_state: bool = True) -> StatePrepCircuit:
    """Return a circuit that prepares the +1 eigenstate of the code w.r.t. the Z or X basis.

    Args:
        code: The CSS code to prepare the state for.
        optimize_depth: If True, optimize the depth of the circuit. This may lead to a higher number of CNOTs.
        zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
    """
    logging.info("Starting heuristic state preparation.")
    if code.Hx is None or code.Hz is None:
        msg = "The code must have both X and Z stabilizers defined."
        raise InvalidCSSCodeError(msg)
    checks = code.Hx.copy() if zero_state else code.Hz.copy()
    rank = mod2.rank(checks)

    def is_reduced() -> bool:
        return bool(len(np.where(np.all(checks == 0, axis=0))[0]) == checks.shape[1] - rank)

    costs = np.array([
        [np.sum((checks[:, i] + checks[:, j]) % 2) for j in range(checks.shape[1])] for i in range(checks.shape[1])
    ])
    costs -= np.sum(checks, axis=0)
    np.fill_diagonal(costs, 1)

    used_qubits: list[np.int_] = []
    cnots: list[tuple[int, int]] = []
    while not is_reduced():
        m = np.zeros((checks.shape[1], checks.shape[1]), dtype=bool)
        m[used_qubits, :] = True
        m[:, used_qubits] = True

        costs_unused = np.ma.array(costs, mask=m)  # type: ignore[no-untyped-call]
        if np.all(costs_unused >= 0) or len(used_qubits) == checks.shape[1]:  # no more reductions possible
            if not used_qubits:  # local minimum => get out by making matrix triangular
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

    s.add(_column_addition_constraint(columns, additions))

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
) -> QuantumCircuit | None:
    """Find the gate optimal circuit for a given check matrix and maximum depth."""
    n = checks.shape[1]
    columns = np.array([
        [[z3.Bool(f"x_{d}_{i}_{j}") for j in range(n)] for i in range(checks.shape[0])] for d in range(max_cnots + 1)
    ])
    n_bits = int(np.ceil(np.log2(n)))
    targets = [z3.BitVec(f"target_{d}", n_bits) for d in range(max_cnots)]
    controls = [z3.BitVec(f"control_{d}", n_bits) for d in range(max_cnots)]
    s = z3.Solver()

    additions = np.array([
        [[z3.And(controls[d] == col_1, targets[d] == col_2) for col_2 in range(n)] for col_1 in range(n)]
        for d in range(max_cnots)
    ])

    # create initial matrix
    columns[0, :, :] = checks.astype(bool)
    s.add(_column_addition_constraint(columns, additions))

    for d in range(1, max_cnots + 1):
        # qubit cannot be control and target at the same time
        s.add(controls[d - 1] != targets[d - 1])

        # control and target must be valid qubits

        if n and (n - 1) != 0 and not ((n & (n - 1) == 0) and n != 0):  # check if n is a power of 2 or 1 or 0
            s.add(z3.ULT(controls[d - 1], n))
            s.add(z3.ULT(targets[d - 1], n))

    # if column is not involved in any addition at certain depth, it is the same as the previous column
    for d in range(1, max_cnots + 1):
        for col in range(n):
            s.add(z3.Implies(targets[d - 1] != col, _symbolic_vector_eq(columns[d, :, col], columns[d - 1, :, col])))

    # assert that final check matrix has n-checks.shape[0] zero columns
    s.add(_final_matrix_constraint(columns))

    if s.check() == z3.sat:
        m = s.model()
        cnots = [(m[controls[d]].as_long(), m[targets[d]].as_long()) for d in range(max_cnots)]
        checks = np.array([
            [bool(m[columns[max_cnots][i][j]]) for j in range(n)] for i in range(checks.shape[0])
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
    if code.Hx is None or code.Hz is None:
        msg = "Code must have both X and Z stabilizers defined."
        raise ValueError(msg)
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
        opt_res = _run_with_timeout(fun, curr_param - 1, timeout=max_timeout)
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
    additional_faults: npt.NDArray[np.int8] | None = None,
) -> list[list[npt.NDArray[np.int8]]]:
    """Return verification stabilizers for the state preparation circuit.

    The method uses an iterative search to find the optimal set of stabilizers by repeatedly computing the optimal circuit for each number of ancillas and cnots. This is repeated for each number of independent correctable errors in the state preparation circuit. Thus the verification circuit is constructed of multiple "layers" of stabilizers, each layer corresponding to a fault set it verifies.

    Args:
        sp_circ: The state preparation circuit to verify.
        x_errors: If True, verify the X errors. If False, verify the Z errors.
        min_timeout: The minimum time to allow each search to run for.
        max_timeout: The maximum time to allow each search to run for.
        max_ancillas: The maximum number of ancillas to allow in each layer verification circuit.
        additional_faults: Faults to verify in addition to the faults propagating in the state preparation circuit.

    Returns:
        A list of stabilizers to verify the state preparation circuit.
    """
    max_errors = sp_circ.max_errors
    layers: list[list[npt.NDArray[np.int8]]] = [[] for _ in range(max_errors)]
    if max_ancillas is None:
        max_ancillas = sp_circ.max_z_measurements if x_errors else sp_circ.max_x_measurements

    sp_circ.compute_fault_sets()
    fault_sets = (
        sp_circ.combine_faults(additional_faults, x_errors)
        if additional_faults is not None
        else sp_circ.x_fault_sets
        if x_errors
        else sp_circ.z_fault_sets
    )

    # Find the optimal circuit for every number of errors in the preparation circuit
    for num_errors in range(1, max_errors + 1):
        logging.info(f"Finding verification stabilizers for {num_errors} errors")
        faults = fault_sets[num_errors]
        assert faults is not None

        if len(faults) == 0:
            logging.info(f"No non-trivial faults for {num_errors} errors")
            layers[num_errors - 1] = []
            continue
        # Start with maximal number of ancillas
        # Minimal CNOT solution must be achievable with these
        num_anc = max_ancillas
        checks = sp_circ.z_checks if x_errors else sp_circ.x_checks
        min_cnots: int = np.min(np.sum(checks, axis=1))
        max_cnots: int = np.sum(checks)

        logging.info(
            f"Finding verification stabilizers for {num_errors} errors with {min_cnots} to {max_cnots} CNOTs using {num_anc} ancillas"
        )

        def fun(num_cnots: int) -> list[npt.NDArray[np.int8]] | None:
            return verification_stabilizers(sp_circ, faults, num_anc, num_cnots, x_errors=x_errors)  # noqa: B023

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
        measurements = [m for m in measurements if np.any(m)]
        num_anc = len(measurements)
        # Iterate backwards to find the minimal number of cnots
        logging.info(f"Finding minimal number of CNOTs for {num_errors} errors")

        def search_cnots(num_cnots: int) -> list[npt.NDArray[np.int8]] | None:
            return verification_stabilizers(sp_circ, faults, num_anc, num_cnots, x_errors=x_errors)  # noqa: B023

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
                return verification_stabilizers(sp_circ, faults, num_anc, num_cnots, x_errors=x_errors)  # noqa: B023

            anc_opt = _run_with_timeout(
                search_anc,
                num_anc - 1,
                timeout=max_timeout,
            )
            if anc_opt is None or (isinstance(anc_opt, str) and anc_opt == "timeout"):
                break
            num_anc -= 1
            measurements = anc_opt
        logging.info(f"Minimal number of ancillas for {num_errors} errors is: {num_anc}")
        layers[num_errors - 1] = measurements

    return layers


def _verification_circuit(
    sp_circ: StatePrepCircuit,
    verification_stabs_fun: Callable[
        [StatePrepCircuit, bool, npt.NDArray[np.int8] | None], list[list[npt.NDArray[np.int8]]]
    ],
    full_fault_tolerance: bool = True,
) -> QuantumCircuit:
    logging.info("Finding verification stabilizers for the state preparation circuit")
    layers_1 = verification_stabs_fun(sp_circ, sp_circ.zero_state, None)
    measurements_1 = [measurement for layer in layers_1 for measurement in layer]

    if full_fault_tolerance:
        additional_errors = _hook_errors(measurements_1)
        layers_2 = verification_stabs_fun(sp_circ, not sp_circ.zero_state, additional_errors)
        measurements_2 = [measurement for layer in layers_2 for measurement in layer]
    else:
        measurements_2 = []

    if sp_circ.zero_state:
        return _measure_ft_stabs(sp_circ, measurements_2, measurements_1, full_fault_tolerance=full_fault_tolerance)
    return _measure_ft_stabs(sp_circ, measurements_1, measurements_2, full_fault_tolerance=full_fault_tolerance)


def gate_optimal_verification_circuit(
    sp_circ: StatePrepCircuit,
    min_timeout: int = 1,
    max_timeout: int = 3600,
    max_ancillas: int | None = None,
    full_fault_tolerance: bool = True,
) -> QuantumCircuit:
    """Return a verified state preparation circuit.

    The verification circuit is a set of stabilizers such that each propagated error in sp_circ anticommutes with some verification stabilizer.

    The method uses an iterative search to find the optimal set of stabilizers by repeatedly computing the optimal circuit for each number of ancillas and cnots. This is repeated for each number of independent correctable errors in the state preparation circuit. Thus the verification circuit is constructed of multiple "layers" of stabilizers, each layer corresponding to a fault set it verifies.

    Args:
        sp_circ: The state preparation circuit to verify.
        min_timeout: The minimum time to allow each search to run for.
        max_timeout: The maximum time to allow each search to run for.
        max_ancillas: The maximum number of ancillas to allow in each layer verification circuit.
        full_fault_tolerance: If True, the verification circuit will be constructed to be fault tolerant to all errors in the state preparation circuit. If False, the verification circuit will be constructed to be fault tolerant only to the type of errors that can cause a logical error. For a logical |0> state preparation circuit, this means the verification circuit will be fault tolerant to X errors but not for Z errors. For a logical |+> state preparation circuit, this means the verification circuit will be fault tolerant to Z errors but not for X errors.
    """

    def verification_stabs_fun(
        sp_circ: StatePrepCircuit,
        zero_state: bool,
        additional_errors: npt.NDArray[np.int8] | None = None,
    ) -> list[list[npt.NDArray[np.int8]]]:
        return gate_optimal_verification_stabilizers(
            sp_circ, zero_state, min_timeout, max_timeout, max_ancillas, additional_errors
        )

    return _verification_circuit(sp_circ, verification_stabs_fun, full_fault_tolerance=full_fault_tolerance)


def heuristic_verification_circuit(
    sp_circ: StatePrepCircuit,
    max_covering_sets: int = 10000,
    find_coset_leaders: bool = True,
    full_fault_tolerance: bool = True,
) -> QuantumCircuit:
    """Return a verified state preparation circuit.

    The method uses a greedy set covering heuristic to find a small set of stabilizers that verifies the state preparation circuit. The heuristic is not guaranteed to find the optimal set of stabilizers.

    Args:
        sp_circ: The state preparation circuit to verify.
        max_covering_sets: The maximum number of covering sets to consider.
        find_coset_leaders: Whether to find coset leaders for the found measurements. This is done using SAT solvers so it can be slow.
        full_fault_tolerance: If True, the verification circuit will be constructed to be fault tolerant to all errors in the state preparation circuit. If False, the verification circuit will be constructed to be fault tolerant only to the type of errors that can cause a logical error. For a logical |0> state preparation circuit, this means the verification circuit will be fault tolerant to X errors but not for Z errors. For a logical |+> state preparation circuit, this means the verification circuit will be fault tolerant to Z errors but not for X errors.
    """

    def verification_stabs_fun(
        sp_circ: StatePrepCircuit, zero_state: bool, additional_errors: npt.NDArray[np.int8] | None = None
    ) -> list[list[npt.NDArray[np.int8]]]:
        return heuristic_verification_stabilizers(
            sp_circ, zero_state, max_covering_sets, find_coset_leaders, additional_errors
        )

    return _verification_circuit(sp_circ, verification_stabs_fun, full_fault_tolerance=full_fault_tolerance)


def heuristic_verification_stabilizers(
    sp_circ: StatePrepCircuit,
    x_errors: bool = True,
    max_covering_sets: int = 10000,
    find_coset_leaders: bool = True,
    additional_faults: npt.NDArray[np.int8] | None = None,
) -> list[list[npt.NDArray[np.int8]]]:
    """Return verification stabilizers for the preparation circuit.

    Args:
        sp_circ: The state preparation circuit to verify.
        x_errors: Whether to find verification stabilizers for X errors. If False, find for Z errors.
        max_covering_sets: The maximum number of covering sets to consider.
        find_coset_leaders: Whether to find coset leaders for the found measurements. This is done using SAT solvers so it can be slow.
        additional_faults: Faults to verify in addition to the faults propagating in the state preparation circuit.
    """
    logging.info("Finding verification stabilizers using heuristic method")
    max_errors = sp_circ.max_errors
    layers: list[list[npt.NDArray[np.int8]]] = [[] for _ in range(max_errors)]
    sp_circ.compute_fault_sets()
    fault_sets = (
        sp_circ.combine_faults(additional_faults, x_errors)
        if additional_faults is not None
        else sp_circ.x_fault_sets
        if x_errors
        else sp_circ.z_fault_sets
    )
    orthogonal_checks = sp_circ.z_checks if x_errors else sp_circ.x_checks
    for num_errors in range(1, max_errors + 1):
        logging.info(f"Finding verification stabilizers for {num_errors} errors")
        faults = fault_sets[num_errors]
        assert faults is not None
        logging.info(f"There are {len(faults)} faults")
        if len(faults) == 0:
            layers[num_errors - 1] = []
            continue

        layers[num_errors - 1] = _heuristic_layer(faults, orthogonal_checks, find_coset_leaders, max_covering_sets)

    return layers


def _covers(s: npt.NDArray[np.int8], faults: npt.NDArray[np.int8]) -> frozenset[int]:
    return frozenset(np.where(s @ faults.T % 2 != 0)[0])


def _set_cover(
    n: int, cands: set[frozenset[int]], mapping: dict[frozenset[int], list[npt.NDArray[np.int8]]]
) -> set[frozenset[int]]:
    universe = set(range(n))
    cover: set[frozenset[int]] = set()

    while universe:
        best = max(cands, key=lambda stab: (len(stab & universe), -np.sum(mapping[stab])))  # type: ignore[operator]
        cover.add(best)
        universe -= best
    return cover


def _extend_covering_sets(
    candidate_sets: set[frozenset[int]], size_limit: int, mapping: dict[frozenset[int], list[npt.NDArray[np.int8]]]
) -> set[frozenset[int]]:
    to_remove: set[frozenset[int]] = set()
    to_add: set[frozenset[int]] = set()
    for c1 in candidate_sets:
        for c2 in candidate_sets:
            if len(to_add) >= size_limit:
                break

            comb = c1 ^ c2
            if c1 == c2 or comb in candidate_sets or comb in to_add or comb in to_remove:
                continue

            mapping[comb].extend([(s1 + s2) % 2 for s1 in mapping[c1] for s2 in mapping[c2]])

            if len(c1 & c2) == 0:
                to_remove.add(c1)
                to_remove.add(c2)
            to_add.add(c1 ^ c2)

    return candidate_sets.union(to_add)


def _heuristic_layer(
    faults: npt.NDArray[np.int8], checks: npt.NDArray[np.int8], find_coset_leaders: bool, max_covering_sets: int
) -> list[npt.NDArray[np.int8]]:
    syndromes = checks @ faults.T % 2
    candidates = np.where(np.any(syndromes != 0, axis=1))[0]
    non_candidates = np.where(np.all(syndromes == 0, axis=1))[0]
    candidate_checks = checks[candidates]
    non_candidate_checks = checks[non_candidates]

    logging.info("Converting Stabilizer Checks to covering sets")
    candidate_sets_ordered = [(_covers(s, faults), s, i) for i, s in enumerate(candidate_checks)]
    mapping = defaultdict(list)
    for cand, _, i in candidate_sets_ordered:
        mapping[cand].append(candidate_checks[i])
    candidate_sets = {cand for cand, _, _ in candidate_sets_ordered}

    logging.info("Finding initial set cover")
    cover = _set_cover(len(faults), candidate_sets, mapping)
    logging.info(f"Initial set cover has {len(cover)} sets")

    def cost(cover: set[frozenset[int]]) -> tuple[int, int]:
        cost1 = len(cover)
        cost2 = sum(np.sum(mapping[stab]) for stab in cover)
        return cost1, cost2

    cost1, cost2 = cost(cover)
    prev_candidates = candidate_sets.copy()

    # find good cover
    improved = True
    while improved and len(candidate_sets) < max_covering_sets:
        improved = False
        # add all symmetric differences to candidates
        candidate_sets = _extend_covering_sets(candidate_sets, max_covering_sets, mapping)
        new_cover = _set_cover(len(faults), candidate_sets, mapping)
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

    # reduce stabilizers in cover
    logging.info(f"Found covering set of size {len(cover)}.")
    if find_coset_leaders and len(non_candidates) > 0:
        logging.info("Finding coset leaders.")
        measurements = []
        for c in cover:
            leaders = [_coset_leader(m, non_candidate_checks) for m in mapping[c]]
            leaders.sort(key=np.sum)
            measurements.append(leaders[0])
    else:
        measurements = [min(mapping[c], key=np.sum) for c in cover]

    return measurements


def _measure_ft_x(qc: QuantumCircuit, x_measurements: list[npt.NDArray[np.int8]], t: int, flags: bool = False) -> None:
    if len(x_measurements) == 0:
        return
    num_x_anc = len(x_measurements)
    x_anc = AncillaRegister(num_x_anc, "x_anc")
    x_c = ClassicalRegister(num_x_anc, "x_c")
    qc.add_register(x_anc)
    qc.add_register(x_c)

    for i, m in enumerate(x_measurements):
        stab = np.where(m != 0)[0]
        if flags:
            measure_flagged(qc, stab, x_anc[i], x_c[i], z_measurement=False, t=t)
        else:
            qc.h(x_anc[i])
            qc.cx([x_anc[i]] * len(stab), stab)
            qc.h(x_anc[i])
            qc.measure(x_anc[i], x_c[i])


def _measure_ft_z(qc: QuantumCircuit, z_measurements: list[npt.NDArray[np.int8]], t: int, flags: bool = False) -> None:
    if len(z_measurements) == 0:
        return
    num_z_anc = len(z_measurements)
    z_anc = AncillaRegister(num_z_anc, "z_anc")
    z_c = ClassicalRegister(num_z_anc, "z_c")
    qc.add_register(z_anc)
    qc.add_register(z_c)

    for i, m in enumerate(z_measurements):
        stab = np.where(m != 0)[0]
        if flags:
            measure_flagged(qc, stab, z_anc[i], z_c[i], z_measurement=True, t=t)
        else:
            qc.cx(stab, [z_anc[i]] * len(stab))
    qc.measure(z_anc, z_c)


def _measure_ft_stabs(
    sp_circ: StatePrepCircuit,
    x_measurements: list[npt.NDArray[np.int8]],
    z_measurements: list[npt.NDArray[np.int8]],
    full_fault_tolerance: bool = True,
) -> QuantumCircuit:
    # Create the verification circuit
    q = QuantumRegister(sp_circ.num_qubits, "q")
    measured_circ = QuantumCircuit(q)
    measured_circ.compose(sp_circ.circ, inplace=True)

    if sp_circ.zero_state:
        _measure_ft_z(measured_circ, z_measurements, t=sp_circ.max_x_errors)
        if full_fault_tolerance:
            _measure_ft_x(measured_circ, x_measurements, flags=True, t=sp_circ.max_x_errors)
    else:
        _measure_ft_x(measured_circ, x_measurements, t=sp_circ.max_z_errors)
        if full_fault_tolerance:
            _measure_ft_z(measured_circ, z_measurements, flags=True, t=sp_circ.max_z_errors)

    return measured_circ


def _vars_to_stab(
    measurement: list[z3.BoolRef | bool], generators: npt.NDArray[np.int8]
) -> npt.NDArray[z3.BoolRef | bool]:
    measurement_stab = _symbolic_scalar_mult(generators[0], measurement[0])
    for i, scalar in enumerate(measurement[1:]):
        measurement_stab = _symbolic_vector_add(measurement_stab, _symbolic_scalar_mult(generators[i + 1], scalar))
    return measurement_stab


def verification_stabilizers(
    sp_circ: StatePrepCircuit,
    fault_set: npt.NDArray[np.int8],
    num_anc: int,
    num_cnots: int,
    x_errors: bool = True,
) -> list[npt.NDArray[np.int8]] | None:
    """Return verification stabilizers for num_errors independent errors in the state preparation circuit using z3.

    Args:
        sp_circ: The state preparation circuit.
        fault_set: The set of errors to verify.
        num_anc: The maximum number of ancilla qubits to use.
        num_cnots: The maximumg number of CNOT gates to use.
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
    solver.add(
        z3.And([
            z3.PbGe([(_odd_overlap(measurement, error), 1) for measurement in measurement_stabs], 1)
            for error in fault_set
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
            v = np.zeros(sp_circ.num_qubits, dtype=np.int8)
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


def _column_addition_constraint(
    columns: npt.NDArray[z3.BoolRef | bool],
    col_add_vars: npt.NDArray[z3.BoolRef],
) -> z3.BoolRef:
    assert len(columns.shape) == 3
    max_depth = col_add_vars.shape[0]
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
    return z3.PbEq(
        [(z3.Not(z3.Or(list(columns[-1, :, col]))), 1) for col in range(columns.shape[2])],
        columns.shape[2] - columns.shape[1],
    )


def _propagate_error(nodes: list[DAGNode], n_qubits: int, x_errors: bool = True) -> PauliList:
    """Propagates a Pauli error through a circuit beginning from first node.

    Args:
        nodes: List of nodes in the circuit in topological order.
        n_qubits: Number of qubits in the circuit.
        x_errors: If True, propagate X errors. Otherwise, propagate Z errors.
    """
    start = nodes[0]
    control = start.qargs[0]._index if x_errors else start.qargs[1]._index  # noqa: SLF001
    error: npt.NDArray[np.int8] = np.array([0] * n_qubits, dtype=np.int8)
    error[control] = 1
    # propagate error through circuit via bfs
    for node in nodes[1:]:
        control = node.qargs[0]._index  # noqa: SLF001
        target = node.qargs[1]._index  # noqa: SLF001
        if x_errors:
            error[target] = (error[target] + error[control]) % 2
        else:
            error[control] = (error[target] + error[control]) % 2
    return error


def _remove_trivial_faults(
    faults: npt.NDArray[np.int8], stabs: npt.NDArray[np.int8], code: CSSCode, x_errors: bool, num_errors: int
) -> npt.NDArray[np.int8]:
    faults = faults.copy()
    logging.info("Removing trivial faults.")
    d_error = code.x_distance if x_errors else code.z_distance
    t_error = max((d_error - 1) // 2, 1)
    t = max((code.distance - 1) // 2, 1)
    max_w = t_error // t
    for i, fault in enumerate(faults):
        faults[i] = _coset_leader(fault, stabs)
    faults = faults[np.where(np.sum(faults, axis=1) > max_w * num_errors)[0]]

    # unique faults
    return np.unique(faults, axis=0)


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


def naive_verification_circuit(sp_circ: StatePrepCircuit) -> QuantumCircuit:
    """Naive verification circuit for a state preparation circuit."""
    if sp_circ.code.Hx is None or sp_circ.code.Hz is None:
        msg = "Code must have stabilizers defined."
        raise ValueError(msg)

    z_measurements = list(sp_circ.code.Hx)
    x_measurements = list(sp_circ.code.Hz)
    reps = sp_circ.max_errors
    return _measure_ft_stabs(sp_circ, z_measurements * reps, x_measurements * reps)


def w_flag_pattern(w: int) -> list[int]:
    """Return the w-flag construction from https://arxiv.org/abs/1708.02246.

    Args:
        w: The number of w-flags to construct.

    Returns:
        The w-flag pattern.
    """
    s1 = [2 * j + 2 for j in reversed(range((w - 4) // 2))]
    s2 = [w - 3, 0]
    s3 = [2 * j + 1 for j in reversed(range((w - 4) // 2))]
    return s1 + s2 + s3 + [w - 2]


def _ancilla_cnot(qc: QuantumCircuit, qubit: Qubit | AncillaQubit, ancilla: AncillaQubit, z_measurement: bool) -> None:
    if z_measurement:
        qc.cx(qubit, ancilla)
    else:
        qc.cx(ancilla, qubit)


def _flag_measure(qc: QuantumCircuit, flag: AncillaQubit, meas_bit: ClBit, z_measurement: bool) -> None:
    if z_measurement:
        qc.h(flag)
    qc.measure(flag, meas_bit)


def _flag_reset(qc: QuantumCircuit, flag: AncillaQubit, z_measurement: bool) -> None:
    qc.reset(flag)
    if z_measurement:
        qc.h(flag)


def _flag_init(qc: QuantumCircuit, flag: AncillaQubit, z_measurement: bool) -> None:
    if z_measurement:
        qc.h(flag)


def _measure_stab_unflagged(
    qc: QuantumCircuit,
    stab: list[Qubit] | npt.NDArray[np.int_],
    ancilla: AncillaQubit,
    measurement_bit: ClBit,
    z_measurement: bool = True,
) -> None:
    if not z_measurement:
        qc.h(ancilla)
        qc.cx([ancilla] * len(stab), stab)
        qc.h(ancilla)
    else:
        qc.cx(stab, [ancilla] * len(stab))
    qc.measure(ancilla, measurement_bit)


def measure_flagged(
    qc: QuantumCircuit,
    stab: list[Qubit] | npt.NDArray[np.int_],
    ancilla: AncillaQubit,
    measurement_bit: ClBit,
    t: int,
    z_measurement: bool = True,
) -> None:
    """Measure a w-flagged stabilizer with the general scheme.

    The measurement is done in place.

    Args:
        qc: The quantum circuit to add the measurement to.
        stab: The qubits to measure.
        ancilla: The ancilla qubit to use for the measurement.
        measurement_bit: The classical bit to store the measurement result of the ancilla.
        t: The number of errors to protect from.
        z_measurement: Whether to measure an X (False) or Z (True) stabilizer.
    """
    w = len(stab)
    if t == 1:
        _measure_stab_one_flagged(qc, stab, ancilla, measurement_bit, z_measurement)
        return

    if w < 3 and t == 2:
        _measure_stab_unflagged(qc, stab, ancilla, measurement_bit, z_measurement)
        return

    if w == 4 and t == 2:
        measure_flagged_4(qc, stab, ancilla, measurement_bit, z_measurement)
        return

    if w == 6 and t == 2:
        measure_flagged_6(qc, stab, ancilla, measurement_bit, z_measurement)
        return

    if w == 8 and t == 2:
        measure_flagged_8(qc, stab, ancilla, measurement_bit, z_measurement)
        return

    if t == 2:
        measure_stab_two_flagged(qc, stab, ancilla, measurement_bit, z_measurement)
        return

    msg = f"Flagged measurement for w={w} and t={t} not implemented."
    raise NotImplementedError(msg)


def _measure_stab_one_flagged(
    qc: QuantumCircuit,
    stab: list[Qubit] | npt.NDArray[np.int_],
    ancilla: AncillaQubit,
    measurement_bit: ClBit,
    z_measurement: bool = True,
) -> None:
    """Measure a 1-flagged stabilizer using an optimized scheme."""
    flag_reg = AncillaRegister(1)
    meas_reg = ClassicalRegister(1)
    qc.add_register(flag_reg)
    qc.add_register(meas_reg)
    flag = flag_reg[0]
    flag_meas = meas_reg[0]
    if not z_measurement:
        qc.h(ancilla)

    _ancilla_cnot(qc, stab[0], ancilla, z_measurement)
    _flag_init(qc, flag, z_measurement)

    _ancilla_cnot(qc, flag, ancilla, z_measurement)

    for q in stab[1:-1]:
        _ancilla_cnot(qc, q, ancilla, z_measurement)

    _ancilla_cnot(qc, flag, ancilla, z_measurement)
    _flag_measure(qc, flag, flag_meas, z_measurement)

    _ancilla_cnot(qc, stab[-1], ancilla, z_measurement)

    if not z_measurement:
        qc.h(ancilla)
    qc.measure(ancilla, measurement_bit)


def measure_stab_two_flagged(
    qc: QuantumCircuit,
    stab: list[Qubit] | npt.NDArray[np.int_],
    ancilla: AncillaQubit,
    measurement_bit: ClBit,
    z_measurement: bool = True,
) -> None:
    """Measure a 2-flagged stabilizer using the scheme of https://arxiv.org/abs/1708.02246 (page 13)."""
    assert len(stab) > 4
    n_flags = (len(stab) + 1) // 2 - 1
    flag_reg = AncillaRegister(n_flags)
    meas_reg = ClassicalRegister(n_flags)

    qc.add_register(flag_reg)
    qc.add_register(meas_reg)

    if not z_measurement:
        qc.h(ancilla)

    _ancilla_cnot(qc, stab[0], ancilla, z_measurement)

    _flag_init(qc, flag_reg[0], z_measurement)
    _ancilla_cnot(qc, flag_reg[0], ancilla, z_measurement)

    _ancilla_cnot(qc, stab[1], ancilla, z_measurement)
    _flag_init(qc, flag_reg[1], z_measurement)
    _ancilla_cnot(qc, flag_reg[1], ancilla, z_measurement)

    cnots = 2
    flags = 2
    for q in stab[2:-2]:
        _ancilla_cnot(qc, q, ancilla, z_measurement)
        cnots += 1
        if cnots % 2 == 0 and cnots < len(stab) - 2:
            _flag_init(qc, flag_reg[flags], z_measurement)
            _ancilla_cnot(qc, flag_reg[flags], ancilla, z_measurement)
        if cnots >= 7 and cnots % 2 == 1:
            _ancilla_cnot(qc, flag_reg[flags - 1], ancilla, z_measurement)
            _flag_measure(qc, flag_reg[flags - 1], meas_reg[flags - 1], z_measurement)
        flags += 1

    _ancilla_cnot(qc, flag_reg[0], ancilla, z_measurement)
    _flag_measure(qc, flag_reg[0], meas_reg[0], z_measurement)

    _ancilla_cnot(qc, stab[-2], ancilla, z_measurement)

    cnots += 1
    if cnots >= 7 and cnots % 2 == 1:
        _ancilla_cnot(qc, flag_reg[flags - 1], ancilla, z_measurement)
        _flag_measure(qc, flag_reg[flags - 1], meas_reg[flags - 1], z_measurement)

    _ancilla_cnot(qc, flag_reg[1], ancilla, z_measurement)
    _flag_measure(qc, flag_reg[1], meas_reg[1], z_measurement)

    _ancilla_cnot(qc, stab[-1], ancilla, z_measurement)

    cnots += 1
    if cnots >= 7 and cnots % 2 == 1:
        _ancilla_cnot(qc, flag_reg[flags - 1], ancilla, z_measurement)
        _flag_measure(qc, flag_reg[flags - 1], meas_reg[flags - 1], z_measurement)
    if not z_measurement:
        qc.h(ancilla)

    qc.measure(ancilla, measurement_bit)


def measure_flagged_4(
    qc: QuantumCircuit,
    stab: list[Qubit] | npt.NDArray[np.int_],
    ancilla: AncillaQubit,
    measurement_bit: ClBit,
    z_measurement: bool = True,
) -> None:
    """Measure a 4-flagged stabilizer using an optimized scheme."""
    assert len(stab) == 4
    flag_reg = AncillaRegister(1)
    meas_reg = ClassicalRegister(1)
    qc.add_register(flag_reg)
    qc.add_register(meas_reg)
    flag = flag_reg[0]
    flag_meas = meas_reg[0]

    if not z_measurement:
        qc.h(ancilla)

    _ancilla_cnot(qc, stab[0], ancilla, z_measurement)
    _flag_init(qc, flag, z_measurement)

    _ancilla_cnot(qc, flag, ancilla, z_measurement)

    _ancilla_cnot(qc, stab[1], ancilla, z_measurement)
    _ancilla_cnot(qc, stab[2], ancilla, z_measurement)

    _ancilla_cnot(qc, flag, ancilla, z_measurement)
    _flag_measure(qc, flag, flag_meas, z_measurement)

    _ancilla_cnot(qc, stab[3], ancilla, z_measurement)

    if not z_measurement:
        qc.h(ancilla)
    qc.measure(ancilla, measurement_bit)


def measure_flagged_6(
    qc: QuantumCircuit,
    stab: list[Qubit] | npt.NDArray[np.int_],
    ancilla: AncillaQubit,
    measurement_bit: ClBit,
    z_measurement: bool = True,
) -> None:
    """Measure a 6-flagged stabilizer using an optimized scheme."""
    assert len(stab) == 6
    flag = AncillaRegister(2)
    meas = ClassicalRegister(2)

    qc.add_register(flag)
    qc.add_register(meas)

    if not z_measurement:
        qc.h(ancilla)

    _ancilla_cnot(qc, stab[0], ancilla, z_measurement)

    _flag_init(qc, flag[0], z_measurement)
    _ancilla_cnot(qc, flag[0], ancilla, z_measurement)

    _ancilla_cnot(qc, stab[1], ancilla, z_measurement)

    _flag_init(qc, flag[1], z_measurement)
    _ancilla_cnot(qc, flag[1], ancilla, z_measurement)

    _ancilla_cnot(qc, stab[2], ancilla, z_measurement)
    _ancilla_cnot(qc, stab[3], ancilla, z_measurement)

    _ancilla_cnot(qc, flag[0], ancilla, z_measurement)
    _flag_measure(qc, flag[0], meas[0], z_measurement)

    _ancilla_cnot(qc, stab[4], ancilla, z_measurement)

    _ancilla_cnot(qc, flag[1], ancilla, z_measurement)
    _flag_measure(qc, flag[1], meas[1], z_measurement)

    _ancilla_cnot(qc, stab[5], ancilla, z_measurement)

    if not z_measurement:
        qc.h(ancilla)
    qc.measure(ancilla, measurement_bit)


def measure_flagged_8(
    qc: QuantumCircuit,
    stab: list[Qubit] | npt.NDArray[np.int_],
    ancilla: AncillaQubit,
    measurement_bit: ClBit,
    z_measurement: bool = True,
) -> None:
    """Measure an 8-flagged stabilizer using an optimized scheme."""
    assert len(stab) == 8
    flag = AncillaRegister(3)
    meas = ClassicalRegister(3)
    qc.add_register(flag)
    qc.add_register(meas)

    if not z_measurement:
        qc.h(ancilla)

    _ancilla_cnot(qc, stab[0], ancilla, z_measurement)

    _flag_init(qc, flag[0], z_measurement)
    _ancilla_cnot(qc, flag[0], ancilla, z_measurement)

    _ancilla_cnot(qc, stab[1], ancilla, z_measurement)

    _flag_init(qc, flag[1], z_measurement)
    _ancilla_cnot(qc, flag[1], ancilla, z_measurement)

    _ancilla_cnot(qc, stab[2], ancilla, z_measurement)
    _ancilla_cnot(qc, stab[3], ancilla, z_measurement)

    _flag_init(qc, flag[2], z_measurement)
    _ancilla_cnot(qc, flag[2], ancilla, z_measurement)

    _ancilla_cnot(qc, stab[4], ancilla, z_measurement)
    _ancilla_cnot(qc, stab[5], ancilla, z_measurement)

    _ancilla_cnot(qc, flag[0], ancilla, z_measurement)
    _flag_measure(qc, flag[0], meas[0], z_measurement)

    _ancilla_cnot(qc, stab[6], ancilla, z_measurement)

    _ancilla_cnot(qc, flag[2], ancilla, z_measurement)
    _flag_measure(qc, flag[2], meas[2], z_measurement)

    _ancilla_cnot(qc, flag[1], ancilla, z_measurement)
    _flag_measure(qc, flag[1], meas[1], z_measurement)

    _ancilla_cnot(qc, stab[7], ancilla, z_measurement)

    if not z_measurement:
        qc.h(ancilla)
    qc.measure(ancilla, measurement_bit)


def _hook_errors(stabs: list[npt.NDArray[np.int8]]) -> npt.NDArray[np.int8]:
    """Assuming CNOTs are executed in ascending order of qubit index, this function gives all the hook errors of the given stabilizer measurements."""
    errors = []
    for stab in stabs:
        error = stab.copy()
        for i in range(len(stab)):
            if stab[i] == 1:
                error[i] = 0
                errors.append(error.copy())
                error[i] = 1
    return np.array(errors)
