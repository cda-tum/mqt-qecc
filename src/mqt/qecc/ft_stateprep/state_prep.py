"""Synthesizing state preparation circuits for CSS codes."""

from __future__ import annotations

from ldpc import mod2
from code import CSSCode
import numpy as np
from qiskit import AncillaRegister, ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.quantum_info import PauliList
from qiskit.converters import circuit_to_dag
from collections import deque
from qiskit.dagcircuit import DAGOutNode
import z3
import multiprocessing

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt
    SymOrBool = z3.BoolRef | bool
    SymVec = list[SymOrBool] | npt.NDArray[SymOrBool]

class StatePrepCircuit:
    """Represents a state preparation circuit for a CSS code."""

    def __init__(self, circ: QuantumCircuit, code: CSSCode, zero_state: bool = True):
        """Initialize a state preparation circuit.

        Args:
            circ: The state preparation circuit.
            code: The CSS code to prepare the state for.
            zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
        """
        self.circ = circ
        self.code = code
        self.zero_state = zero_state
        self.checks = code.Hx.copy() if zero_state else code.Hz.copy()
        self.orthogonal_checks = np.vstack((code.Hz.copy(), code.Lz.copy())) if zero_state else np.vstack((code.Hx.copy(), code.Lx.copy()))
        self.num_qubits = circ.num_qubits
        self.fault_sets = [None for _ in range((code.distance-1) // 2 + 1)]
        self.max_measurements = len(self.orthogonal_checks)

    def compute_fault_set(self, n_errors=1, reduce=True) -> npt.NDArray[np.bool_]:
        """Compute the fault set of the state.

        Args:
            state: The stabilizer state to compute the fault set for.

        Returns:
            The fault set of the state.
        """
        if self.fault_sets[n_errors] is not None:
            return self.fault_sets[n_errors]

        dag = circuit_to_dag(self.circ)
        for node in dag.front_layer(): # remove hadamards
            dag.remove_op_node(node)
        faults = []
        # propagate every error before a control
        for node in dag.topological_op_nodes():
            error = _propagate_error(dag, node, zero_state=self.zero_state)
            faults.append(error)
        faults = _remove_stabilizer_equivalent_faults(faults, self.code.Hx)

        # combine faults to generate higher weight errors
        for k in range(1, n_errors):
            new_faults = []
            # Append single_qubit errors
            for i, fault in enumerate(faults):
                for j, fault2 in enumerate(faults[i+1:]):
                    new_faults.append((fault+fault2) % 2)

                for j in range(self.num_qubits):
                    single_fault = np.array([0]*self.num_qubits)
                    single_fault[j] = 1
                    new_faults.append((fault+single_fault) % 2)
            faults = new_faults
        # reduce faults by stabilizer
        reduced = True
        while reduced:
            reduced = False
            for i, fault in enumerate(faults):
                for check in self.checks:
                    reduced_fault = (fault+check) % 2
                    if np.sum(fault) > np.sum(reduced_fault):
                        faults[i] = reduced_fault
                        reduced = True
                        break

        # remove trivial faults
        faults = np.array(faults)
        faults = faults[np.where(np.sum(faults, axis=1) > (self.code.distance-1) // 2)[0]]
        # remove stabilizer equivalent faults

        faults = _remove_stabilizer_equivalent_faults(faults, self.checks)
        self.fault_sets[n_errors] = np.array(faults)
        return faults


class NDFTStatePrepCircuit:
    """Non-deterministic fault-tolerant state preparation circuit for a CSS code."""

    def __init__(self, circ: QuantumCircuit, code: CSSCode, zero_state: bool = True):
        """Initialize a state preparation circuit.

        Args:
            circ: The state preparation circuit.
            code: The CSS code to prepare the state for.
            zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
        """
        self.circ = circ
        self.code = code
        self.zero_state = zero_state
        
        
def heuristic_prep_circuit(code: CSSCode, optimize_depth: bool = True, zero_state: bool = True):
    """Return a circuit that prepares the +1 eigenstate of the code w.r.t. the Z or X basis.

    Args:
        code: The CSS code to prepare the state for.
        zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
    """
    checks = code.Hx.copy() if zero_state else code.Hz.copy()
    rank = mod2.rank(checks)

    def is_reduced():
        return len(np.where(np.all(checks == 0, axis=0))[0]) == checks.shape[1]-rank

    costs = np.array([[np.sum((checks[:, i] + checks[:, j]) % 2)
                       for j in range(checks.shape[1])]
                      for i in range(checks.shape[1])])
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

    circ = _build_circuit_from_list_and_checks(cnots, checks, zero_state)
    return StatePrepCircuit(circ, code, zero_state)


def _generate_circ_with_bounded_depth(checks, max_depth, zero_state=True) -> np.array:
    assert max_depth > 0, "max_depth should be greater than 0"
    columns = np.array([[[z3.Bool(f'x_{d}_{i}_{j}')
                          for j in range(checks.shape[1])]
                         for i in range(checks.shape[0])]
                        for d in range(max_depth+1)])

    additions = np.array([[[z3.Bool(f'add_{d}_{i}_{j}')
                            for j in range(checks.shape[1])]
                           for i in range(checks.shape[1])]
                          for d in range(max_depth)])
    n_cols = checks.shape[1]
    s = z3.Solver()

    # create initial matrix
    columns[0, :, :] = checks.astype(bool)

    s.add(_column_addition_contraint(columns, additions))

    # qubit can be involved in at most one addition at each depth
    for d in range(max_depth):
        for col in range(n_cols):
            s.add(z3.PbLe([(additions[d, col_1, col], 1)
                           for col_1 in range(n_cols)
                           if col != col_1] +
                          [(additions[d, col, col_2], 1)
                           for col_2 in range(n_cols)
                           if col != col_2], 1
                          )
                  )

    # if column is not involved in any addition at certain depth, it is the same as the previous column
    for d in range(1, max_depth+1):
        for col in range(n_cols):
            s.add(
                z3.Implies(
                    z3.Not(z3.Or(list(np.delete(additions[d-1, :, col], [col])) +
                                 list(np.delete(additions[d-1, col, :], [col])))),
                    _symbolic_vector_eq(columns[d, :, col], columns[d-1, :, col])))

    s.add(_final_matrix_constraint(columns))

    if s.check() == z3.sat:
        m = s.model()
        additions = [(i, j) for d in range(max_depth)
                     for j in range(checks.shape[1])
                     for i in range(checks.shape[1]) if m[additions[d, i, j]]]

        checks = np.array([[bool(m[columns[max_depth, i, j]]) for j in range(checks.shape[1])]
                           for i in range(checks.shape[0])])

        return _build_circuit_from_list_and_checks(additions, checks, zero_state=zero_state)

    return None


def _generate_circ_with_bounded_gates(checks, max_cnots: int, zero_state=True):
    """Find the gate optimal circuit for a given check matrix and maximum depth."""
    columns = np.array([[[z3.Bool(f'x_{d}_{i}_{j}')
                          for j in range(checks.shape[1])]
                         for i in range(checks.shape[0])]
                        for d in range(max_cnots+1)])
    n_bits = int(np.ceil(np.log2(checks.shape[1])))
    targets = [z3.BitVec(f'target_{d}', n_bits) for d in range(max_cnots)]
    controls = [z3.BitVec(f'control_{d}', n_bits) for d in range(max_cnots)]
    s = z3.Solver()

    additions = np.array([[[z3.And(controls[d] == col_1, targets[d] == col_2)
                            for col_2 in range(checks.shape[1])]
                           for col_1 in range(checks.shape[1])]
                          for d in range(max_cnots)])

    # create initial matrix
    columns[0, :, :] = checks.astype(bool)
    s.add(_column_addition_contraint(columns, additions))

    for d in range(1, max_cnots+1):
        # qubit cannot be control and target at the same time
        s.add(controls[d-1] != targets[d-1])

        # control and target must be valid qubits
        if checks.shape[1] and (checks.shape[1]-1) != 0:
            s.add(z3.ULT(controls[d-1], checks.shape[1]))
            s.add(z3.ULT(targets[d-1], checks.shape[1]))



    # if column is not involved in any addition at certain depth, it is the same as the previous column
    for d in range(1, max_cnots+1):
            for col in range(checks.shape[1]):
                s.add(z3.Implies(
                    targets[d-1] != col,
                    _symbolic_vector_eq(columns[d, :, col], columns[d-1, :, col])))

    # assert that final check matrix has checks.shape[1]-checks.shape[0] zero columns
    s.add(_final_matrix_constraint(columns))

    if s.check() == z3.sat:
        m = s.model()
        additions = [(m[controls[d]].as_long(), m[targets[d]].as_long()) for d in range(max_cnots)]
        checks = np.array([[bool(m[columns[max_cnots][i][j]]) for j in range(checks.shape[1])] for i in range(checks.shape[0])]).astype(int)
        return _build_circuit_from_list_and_checks(additions, checks, zero_state=zero_state)

    return None


def _optimal_circuit(code: CSSCode, prep_func, zero_state: bool=True, min_param=1, max_param=10, min_timeout=1, max_timeout=3600) -> QuantumCircuit:
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
    circ = None

    circ, curr_param = iterative_search_with_timeout(lambda param: prep_func(checks, param), min_param, max_param, min_timeout, max_timeout)

    if not circ:
        return None

    # Solving a SAT instance is much faster than proving unsat in this case
    # so we iterate backwards until we find an unsat instance or hit a timeout
    curr_param -= 1
    while True:
        res = _run_with_timeout(prep_func, checks, curr_param, timeout=max_timeout)

        if res is not None or res == "timeout":
            break
        circ = res
        curr_param -= 1

    return StatePrepCircuit(circ, code, zero_state)


def depth_optimal_prep_circuit(code: CSSCode, zero_state: bool=True, min_depth=1, max_depth=10, min_timeout=1, max_timeout=3600) -> QuantumCircuit:
    """Synthesize a state preparation circuit for a CSS code that minimizes the circuit depth.

    Args:
        code: The CSS code to prepare the state for.
        zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
        min_depth: minimum depth to start with
        max_depth: maximum depth to reach
        min_timeout: minimum timeout to start with
        max_timeout: maximum timeout to reach
    """
    return _optimal_circuit(code, _generate_circ_with_bounded_depth, zero_state, min_depth, max_depth, min_timeout, max_timeout)


def gate_optimal_prep_circuit(code: CSSCode, zero_state: bool=True, min_gates=1, max_gates=10, min_timeout=1, max_timeout=3600) -> QuantumCircuit:
    """Synthesize a state preparation circuit for a CSS code that minimizes the number of gates.

    Args:
        code: The CSS code to prepare the state for.
        zero_state: If True, prepare the +1 eigenstate of the Z basis. If False, prepare the +1 eigenstate of the X basis.
        min_gates: minimum number of gates to start with
        max_gates: maximum number of gates to reach
        min_timeout: minimum timeout to start with
        max_timeout: maximum timeout to reach
    """
    return _optimal_circuit(code, _generate_circ_with_bounded_gates, zero_state, min_gates, max_gates, min_timeout, max_timeout)


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


def _run_with_timeout(func, *args, timeout: int=10):
    """Run a function with a timeout.

    If the function does not complete within the timeout, return None.

    Args:
        func: The function to run.
        args: The arguments to pass to the function.
        timeout: The maximum time to allow the function to run for in seconds.
    """
    manager = multiprocessing.Manager()
    return_list = manager.list()
    p = multiprocessing.Process(target=lambda: return_list.append(func(*args)))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
        return "timeout"
    return return_list[0]


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
    while curr_timeout <= max_timeout:
        while curr_param <= max_param:
            res = _run_with_timeout(fun, curr_param, timeout=curr_timeout)
            if res is not None or (isinstance(res, str) and res != "timeout"):
                return res, curr_param
            curr_param = param_type(curr_param*param_factor)
        curr_timeout *= timeout_factor
        curr_param = min_param
    return None, max_param


def gate_optimal_verification_stabilizers(sp_circ: StatePrepCircuit, n_errors: int = 1, min_timeout=10, max_timeout=3600) -> list[npt.NDArray[int]]:
    """Return a verification circuit for the state preparation circuit.

    Args:
        n_errors: The number of errors to detect.
        reduce: If True, reduce the fault set before computing the verification circuit.

    Returns:
        The verification circuit.
    """
    layers = [None for _ in range(n_errors)]
    # Find the optimal circuit for every number of errors in the preparation circuit
    for num_errors in range(1, (sp_circ.code.distance-1) // 2 + 1):
        # Start with maximal number of ancillas
        # Minimal CNOT solution must be achievable with these
        num_anc = sp_circ.max_measurements
        min_cnots = np.min(np.sum(sp_circ.orthogonal_checks, axis=1))
        max_cnots = np.sum(sp_circ.orthogonal_checks)

        measurements, num_cnots = iterative_search_with_timeout(lambda num_cnots: verification_stabilizers(sp_circ, num_anc, num_cnots, num_errors), min_cnots, max_cnots, min_timeout, max_timeout)

        if measurements is None or (isinstance(measurements, str) and measurements == "timeout"):
            return None  # No solution found

        # If any measurements are unused we can reduce the number of ancillas at least by that
        num_anc = np.sum([np.any(m) for m in measurements])
        measurements = [m for m in measurements if np.any(m)]

        # Iterate backwards to find the minimal number of cnots
        num_cnots -= 1
        while num_cnots > 0:
            res = verification_stabilizers(sp_circ, num_anc, num_cnots, num_errors)
            if res is None or res == "timeout":
                break
            num_cnots -= 1
            measurements = res

        # If the number of CNOTs is minimal, we can reduce the number of ancillas
        num_anc -= 1
        while num_anc > 0:
            res = verification_stabilizers(sp_circ, num_anc, num_cnots, num_errors)
            if res is None or res == "timeout":
                break
            num_anc -= 1
            measurements = res

        layers[num_errors-1] = measurements

    return layers

def gate_optimal_verification_circuit(sp_circ: StatePrepCircuit, n_errors: int = 1, min_timeout: int=10, max_timeout: int=3600) -> QuantumCircuit:
    """Return a verified state preparation circuit."""
    layers = gate_optimal_verification_stabilizers(sp_circ, n_errors, min_timeout, max_timeout)
    if layers is None:
        return None

    measured_circ = _measure_stabs(sp_circ.circ, [measurement for layer in layers for measurement in layer])
    return NDFTStatePrepCircuit(sp_circ.code, measured_circ, zero_state=sp_circ.zero_state)


def _measure_stabs(circ: QuantumCircuit, measurements: list(npt.NDArray[np.int_])) -> QuantumCircuit:
    # Create the verification circuit
    num_anc = len(measurements)
    q = QuantumRegister(circ.num_qubits, "q")
    anc = AncillaRegister(num_anc, "anc")
    c = ClassicalRegister(num_anc, "c")
    measured_circ = QuantumCircuit(q, anc, c)

    measured_circ.compose(circ, inplace=True)
    current_anc = 0
    for measurement in measurements:
        if not self.zero_state:
            measured_circ.h(q)
        for qubit in np.where(measurement == 1)[0]:
            if self.zero_state:
                measured_circ.cx(q[qubit], anc[current_anc])
            else:
                measured_circ.cx(anc[current_anc], q[qubit])
        if not self.zero_state:
            measured_circ.h(q)
        current_anc += 1
    return measured_circ

def verification_stabilizers(self, num_anc, num_cnots, num_errors):
    """Return a verification circuit for the state preparation circuit.

    Args:
        num_anc: The maximum number of ancilla qubits to use.
        num_cnots: The maximumg number of CNOT gates to use.
        num_errors: The number of errors occur in the state prep circuit.
    """
    # Measurements are written as sums of generators
    # The variables indicate which generators are non-zero in the sum
    measurement_vars = [[z3.Bool("m_{}_{}".format(anc, i))
                         for i in range(self.max_measurements)]
                        for anc in range(num_anc)]
    solver = z3.Solver()

    def vars_to_stab(measurement):
        measurement_stab = _symbolic_scalar_mult(self.orthogonal_checks[0], measurement[0])
        for i, scalar in enumerate(measurement[1:]):
            measurement_stab = _symbolic_vector_add(measurement_stab, _symbolic_scalar_mult(self.orthogonal_checks[i+1], scalar))
        return measurement_stab

    measurement_stabs = [vars_to_stab(vars_) for vars_ in measurement_vars]

    # assert that each error is detected
    errors = self._compute_fault_set(num_errors)
    solver.add(z3.And([z3.PbGe([(_odd_overlap(measurement, error), 1)
                                for measurement in measurement_stabs], 1)
                       for error in errors]))

    # assert that not too many CNOTs are used
    solver.add(z3.PbLe([(measurement[q], 1)
                        for measurement in measurement_stabs
                        for q in range(self.num_qubits)],
                       num_cnots))

    if solver.check() == z3.sat:
        model = solver.model()
        # Extract stabilizer measurements from model
        actual_measurements = []
        for m in measurement_vars:
            v = np.zeros(self.num_qubits, dtype=int)
            for g in range(self.max_measurements):
                if model[m[g]]:
                    v += self.orthogonal_checks[g]
            actual_measurements.append(v % 2)

        return np.array(actual_measurements)

    return None

    
def _symbolic_scalar_mult(v: npt.NDArray[np.int_], a: z3.BoolRef | bool):
    """Multiply a concrete vector by a symbolic scalar."""
    return [a if s == 1 else False for s in v]


def _symbolic_vector_add(v1: SymVec, v2: SymVec):
    """Add two symbolic vectors."""
    if v1 is None:
        return v2
    if v2 is None:
        return v1

    v_new = [False for _ in range(len(v1))]
    for i in range(len(v1)):
        # If one of the elements is a bool, we can simplify the expression
        v1_i_is_bool = isinstance(v1[i], bool)
        v2_i_is_bool = isinstance(v2[i], bool)
        if v1_i_is_bool:
            if v1[i]:
                v_new[i] = z3.Not(v2[i]) if not v2_i_is_bool else not v2[i]
            else:
                v_new[i] = v2[i]

        elif v2_i_is_bool:
            if v2[i]:
                v_new[i] = z3.Not(v1[i]) if not v1_i_is_bool else not v1[i]
            else:
                v_new[i] = v1[i]

        else:
            v_new[i] = z3.Xor(v1[i], v2[i])

    return v_new


def _odd_overlap(v_sym: SymVec, v_con: npt.NDArray[np.int_]):
    """Return True if the overlap of symbolic vector with constant vector is odd."""
    return z3.PbEq([(v_sym[i], 1) for i, c in enumerate(v_con) if c == 1], 1)


def _symbolic_vector_eq(v1: SymVec, v2: SymVec):
    """Return assertion that two symbolic vectors should be equal."""
    return z3.And([v1[i] == v2[i] for i in range(len(v1))])


def _column_addition_contraint(columns, col_add_vars):
    max_depth = col_add_vars.shape[0]
    n_cols = col_add_vars.shape[2]

    constraints = []
    for d in range(1, max_depth+1):
        for col_1 in range(n_cols):
            for col_2 in range(col_1+1, n_cols):
                col_sum = _symbolic_vector_add(columns[d-1, :, col_1], columns[d-1, :, col_2])

                # encode col_2 += col_1
                constraints.append(z3.Implies(col_add_vars[d-1, col_1, col_2],
                                              z3.And(_symbolic_vector_eq(columns[d, :, col_2], col_sum),
                                                     _symbolic_vector_eq(columns[d, :, col_1], columns[d-1, :, col_1]))))

                # encode col_1 += col_2
                constraints.append(z3.Implies(col_add_vars[d-1, col_2, col_1],
                                              z3.And(_symbolic_vector_eq(columns[d, :, col_1], col_sum),
                                                     _symbolic_vector_eq(columns[d, :, col_2], columns[d-1, :, col_2]))))

    return z3.And(constraints)


def _final_matrix_constraint(columns):
    return z3.PbEq([(z3.Not(z3.Or(list(columns[-1, :, col]))), 1) for col in range(columns.shape[2])], columns.shape[2]-columns.shape[1])


def _propagate_error(dag, node, zero_state=True) -> PauliList:
    """Propagates a Pauli error through a circuit beginning from control of node."""
    control = node.qargs[0]._index
    error = np.array([0]*dag.num_qubits())
    error[control] = 1
    # propagate error through circuit via bfs
    q = deque([node])
    visited = set()
    while q:
        node = q.popleft()
        if node in visited or isinstance(node, DAGOutNode):
            continue
        control = node.qargs[0]._index
        target = node.qargs[1]._index
        if zero_state:
            error[target] = (error[target] + error[control])%2
        else:
            error[control] = (error[target] + error[control])%2
        for succ in dag.successors(node):
            q.append(succ)
    return error


def _remove_stabilizer_equivalent_faults(faults: np.array, stabilizers: np.array) -> np.array:
    """Remove stabilizer equivalent faults from a list of faults."""
    faults = faults.copy()
    stabilizers = stabilizers.copy()
    removed = set()

    for i, f1 in enumerate(faults):
        if i in removed:
            continue
        A = np.vstack((stabilizers, f1))
        A = mod2.reduced_row_echelon(A)[0]
        if np.all(A[-1] == 0):
            removed.add(i)
            continue

        for j, f2 in enumerate(faults[i+1:]):
            if j+i+1 in removed:
                continue
            S = np.vstack((A, f2))
            S = mod2.reduced_row_echelon(S)[0]

            if np.all(S[-1] == 0):
                removed.add(j+i+1)

    indices = np.array(list(set(range(len(faults)))-removed))
    if len(indices) == 0:
        return np.array([])
    return np.array(faults)[indices]
