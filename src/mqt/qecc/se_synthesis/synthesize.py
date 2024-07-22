from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
import z3
from qiskit import AncillaRegister, ClassicalRegister, QuantumCircuit, QuantumRegister

if TYPE_CHECKING:
    import numpy.typing as npt

    from ..codes import CSSCode


class OptimalSyndromeExtractionEncoder:
    """Encoder instance for optimal synthesis of decoding circuits for a CSS code.

    Attributes:
        solver (z3.Solver): The z3 solver instance.
        qubit_variables (list): A list of lists of lists of z3 variables. The first index corresponds to the qubit index, the second to the time step and the third to the variable index.
        n_xchecks (int): The number of x-checks.
        n_zchecks (int): The number of z-checks.
        x_checks (np.array): The x-check matrix of the code.
        z_checks (np.array): The z-check matrix of the code.
        T (int): The maximal (CNOT) depth of the encoder circuit.
        x_vars (list): Variables for the x-checks. x_vars[check][qubit][t] is a z3 variable that is true if and only if the CNOT on this qubit of this x-check is performed at time t.]
        z_vars (list): Variables for the z-checks. z_vars[check][qubit][t] is a z3 variable that is true if and only if the CNOT on this qubit of this z-check is performed at time t.
        x_lt_z_vars (list): Variables for the x-checks. x_lt_z_vars[check][qubit][t] is a z3 variable that is true if and only if the CNOT on this qubit of this x-check is performed before the CNOT on this qubit of this z-check.
        z_lt_x_vars (list): Variables for the z-checks. z_lt_x_vars[check][qubit][t] is a z3 variable that is true if and only if the CNOT on this qubit of this z-check is performed before the CNOT on this qubit of this x-check.
        overlaps (dict): Dictionary of the form overlaps[(i, j)] = [q1, q2, ...] where q1, q2, ... are the qubits that are involved in both the i-th x-check and the j-th z-check.
    """

    def __init__(
        self,
        code: CSSCode | tuple[npt.NDArray[np.int8], npt.NDArray[np.int8]],
        T: int,
        dependency_graph: npt.NDArray[np.int8] = None,
    ) -> None:
        """Construct encoder instance for optimal synthesis of decoding circuits for an ldpc code.

        Args:
            code (CSSCode | tuple): The CSS code for which the syndrome extraction circuit is to be synthesized. If a tuple is given, it is assumed to be of the form (x_checks, z_checks) where x_checks and z_checks are the x-check and z-check matrices of the code, respectively.
            T (int): The maximal (CNOT) depth of the decoder circuit.
            dependency_graph (np.array): Graph where nodes are the qubits and an edge between two qubits indicates that they can't be measured at the same time.
        """
        self.solver = z3.Solver()

        if isinstance(code, tuple):
            x_checks, z_checks = code
        else:
            x_checks = code.Hx
            z_checks = code.Hz

        self.dependency_graph = dependency_graph

        self.n_qubits = x_checks.shape[1]
        self.qubit_variables = [[[] for t in range(T)] for _ in range(x_checks.shape[1])]
        self.n_xchecks = x_checks.shape[0]
        self.n_zchecks = z_checks.shape[0]
        self.x_checks = x_checks.copy()
        self.z_checks = z_checks.copy()
        self.T = T
        self.x_vars = []
        self.z_vars = []
        self.x_lt_z_vars = []
        self.z_lt_x_vars = []
        self.overlaps = {}
        self._cnot_order_constraint_x = defaultdict(dict)
        self._cnot_order_constraint_z = defaultdict(dict)
        self._circuit = None

        for i in range(x_checks.shape[0]):
            check = {}
            for j in range(x_checks.shape[1]):
                cnots = []
                if x_checks[i, j] == 1:
                    for t in range(T):
                        var = z3.Bool(f"x_{i}_{j}_{t}")
                        cnots.append(var)
                        self.qubit_variables[j][t].append(var)
                if len(cnots) > 0:
                    check[j] = cnots

            self.x_vars.append(check)

        for i in range(z_checks.shape[0]):
            check = {}
            for j in range(z_checks.shape[1]):
                cnots = []
                if z_checks[i, j] == 1:
                    for t in range(T):
                        var = z3.Bool(f"z_{i}_{j}_{t}")
                        cnots.append(var)
                        self.qubit_variables[j][t].append(var)
                if len(cnots) > 0:
                    check[j] = cnots
            self.z_vars.append(check)

        for i, x_check in enumerate(x_checks):
            for j, z_check in enumerate(z_checks):
                qubits = []
                for k in range(len(x_check)):
                    if x_check[k] == 1 and z_check[k] == 1:
                        qubits.append(k)
                self.overlaps[(i, j)] = qubits

    def _cnot_exists_constraint(self, check: list, cnot: int):
        """Create a constraint that the CNOT is performed exactly once in the circuit."""
        return z3.PbEq([(check[cnot][t], 1) for t in range(self.T)], 1)
        # return z3.PbEq([(x, 1) for x in [check[cnot][t] for t in range(self.T)]], 1)

    def _cnot_comes_before_cnot(self, check1: list[z3.Bool], check2: list[z3.Bool], cnot):
        """Create a constraint that the CNOT of check1 comes before the CNOT of check2."""
        return z3.Or([z3.And(check1[cnot][k], z3.Or(check2[cnot][k + 1 :])) for k in range(self.T)])

    def _assert_check_order_constraint(self, i: int, j: int) -> None:
        """Assert that CNOTs of x and z checks alternate in pairs."""
        x_constraints = list(self._cnot_order_constraint_x[(i, j)].values())
        if len(x_constraints) == 0:
            return

        # If an odd number of x-checks appear before the z-checks the commutation relation is violated.
        formula = x_constraints[0]
        for constr in x_constraints[1:]:
            formula = z3.Xor(formula, constr)
        self.solver.add(z3.Not(formula))

        z_constraints = list(self._cnot_order_constraint_z[(j, i)].values())
        if len(z_constraints) == 0:
            return

        # If an odd number of z-checks appear before the x-checks the commutation relation is violated.
        formula = z_constraints[0]
        for constr in z_constraints[1:]:
            formula = z3.Xor(formula, constr)
        self.solver.add(z3.Not(formula))

    def _cnots_not_overlapping_constraint(self, qubit, t):
        """Create a constraint that ensures a qubit is not involved in two CNOTs at a time."""
        return z3.PbLe([(x, 1) for x in self.qubit_variables[qubit][t]], 1)

    def _assert_cnot_on_check_constraint(self, check) -> None:
        """Assert that the CNOTs of a check are performed at different times."""
        for t in range(self.T):
            self.solver.add(z3.PbLe([(x, 1) for x in [check[qubit][t] for qubit in check]], 1))

    def _assert_dependency_constraints(self) -> None:
        for i, j in np.argwhere(self.dependency_graph):
            for t in range(self.T):
                for var_i in self.qubit_variables[i][t]:
                    for var_j in self.qubit_variables[j][t]:
                        self.solver.add(z3.Not(z3.And(var_i, var_j)))

    def _encode_constraints(self) -> None:
        # create variables for cnot order constraints
        for i, j in self.overlaps:
            for cnot in self.overlaps[(i, j)]:
                self._cnot_order_constraint_x[(i, j)][cnot] = z3.Bool(f"xltz_{i}_{j}_{cnot}")
                self._cnot_order_constraint_z[(j, i)][cnot] = z3.Bool(f"zltx_{j}_{i}_{cnot}")
                self.solver.add(
                    self._cnot_order_constraint_x[(i, j)][cnot]
                    == self._cnot_comes_before_cnot(self.x_vars[i], self.z_vars[j], cnot)
                )
                self.solver.add(
                    self._cnot_order_constraint_z[(j, i)][cnot]
                    == self._cnot_comes_before_cnot(self.z_vars[j], self.x_vars[i], cnot)
                )

        for check in self.x_vars:
            self._assert_cnot_on_check_constraint(check)
            for cnot in check:
                self.solver.add(self._cnot_exists_constraint(check, cnot))
        for check in self.z_vars:
            self._assert_cnot_on_check_constraint(check)
            for cnot in check:
                self.solver.add(self._cnot_exists_constraint(check, cnot))

        # assert that non-commuting cnots are ordered correctly
        for i, j in self.overlaps:
            self._assert_check_order_constraint(i, j)

        # assert that each qubit is used at most once per timestep
        for i, _qubit in enumerate(self.qubit_variables):
            for t in range(self.T):
                self.solver.add(self._cnots_not_overlapping_constraint(i, t))

        if self.dependency_graph is not None:
            self._assert_dependency_constraints()

    def solve(self):
        """Solve the problem."""
        self._encode_constraints()

        result = self.solver.check()
        if result != z3.sat:
            return result

        self._circuit = self._extract_circuit()
        return str(result)

    def get_schedule(self) -> list[list[tuple[int, int, str]]]:
        """Return the schedule.

        Returns:
            list[list[tuple[int, str]]]: The schedule given as a list of lists of tuples. Each list corresponds to a
                timestep and each tuple corresponds to a CNOT gate. The first element of the tuple is the data qubit on which the CNOT is applied and the second element is the index of the check qubit. The third element is either 'X' or 'Z' depending on whether the CNOT is an X-check or a Z-check.
        """
        schedule = []
        for t in range(self.T):
            timestep = []
            for qubit in range(self.x_checks.shape[1]):
                for i, check in enumerate(self.x_checks):
                    if check[qubit] == 1 and self.solver.model()[self.x_vars[i][qubit][t]]:
                        timestep.append((qubit, i, "X"))

                for i, check in enumerate(self.z_checks):
                    if check[qubit] == 1 and self.solver.model()[self.z_vars[i][qubit][t]]:
                        timestep.append((qubit, i, "Z"))

            schedule.append(timestep)
        return schedule

    def get_circuit(self):
        """Return the circuit."""
        if self._circuit is None:
            msg = "No circuit available. Call solve() first."
            raise ValueError(msg)
        return self._circuit

    def _extract_circuit(self):
        q = QuantumRegister(len(self.qubit_variables), "q")
        x_c = ClassicalRegister(self.n_xchecks, "x_c")
        z_c = ClassicalRegister(self.n_zchecks, "z_c")
        x_anc = AncillaRegister(self.n_xchecks, "x")
        z_anc = AncillaRegister(self.n_zchecks, "z")
        circuit = QuantumCircuit(q, x_c, z_c, x_anc, z_anc)

        # hadamard gates for the x-checks
        for anc in x_anc:
            circuit.h(anc)

        for t in range(self.T):
            for qubit in range(self.x_checks.shape[1]):
                for i, check in enumerate(self.x_checks):
                    if check[qubit] == 1 and self.solver.model()[self.x_vars[i][qubit][t]]:
                        circuit.cx(x_anc[i], q[qubit])

                for i, check in enumerate(self.z_checks):
                    if check[qubit] == 1 and self.solver.model()[self.z_vars[i][qubit][t]]:
                        circuit.cx(q[qubit], z_anc[i])

        # measurements
        for anc, c in zip(x_anc, x_c):
            circuit.h(anc)
            circuit.measure(anc, c)

        for anc, c in zip(z_anc, z_c):
            circuit.measure(anc, c)

        return circuit
