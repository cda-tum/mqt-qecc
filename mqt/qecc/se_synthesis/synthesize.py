from collections import defaultdict

import numpy as np
import z3
from qiskit import (AncillaRegister, ClassicalRegister, QuantumCircuit,
                    QuantumRegister)


class SyndromeExtractionEncoder:
    """Encoder instance for optimal synthesis of decoding circuits for an ldpc code.

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

    def __init__(self, x_checks: np.array, z_checks: np.array, T: int):
        """Construct encoder instance for optimal synthesis of decoding circuits for an ldpc code.

        Args:
            x_checks (np.array): The x-check matrix of the code.
            z_checks (np.array): The z-check matrix of the code.
            T (int): The maximal (CNOT) depth of the encoder circuit.
        """
        self.solver = z3.Solver()
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
                        var = z3.Bool(f'x_{i}_{j}_{t}')
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
                        var = z3.Bool(f'z_{i}_{j}_{t}')
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
        return z3.Or([check[cnot][t] for t in range(self.T)])
        # return z3.PbEq([(x, 1) for x in [check[cnot][t] for t in range(self.T)]], 1)

    def _cnot_comes_before_cnot(self, x_check: list[z3.Bool], z_check: list[z3.Bool], cnot):
        return z3.Or([z3.And(x_check[cnot][k], z3.Or(z_check[cnot][k+1:])) for k in range(self.T)])

    def _assert_check_order_constraint(self, i, j):
        self.solver.add(z3.Or(z3.And(list(self._cnot_order_constraint_x[(i, j)].values())), z3.And(list(self._cnot_order_constraint_z[(j, i)].values()))))

    def _cnots_not_overlapping_constraint(self, qubit, t):
        return z3.PbLe([(x, 1) for x in self.qubit_variables[qubit][t]], 1)

    def _assert_cnot_on_check_constraint(self, check):
        for t in range(self.T):
            self.solver.add(z3.PbLe([(x, 1) for x in [check[qubit][t] for qubit in check.keys()]], 1))

    def _encode_constraints(self):
        # create variables for cnot order constraints
        for i, j in self.overlaps.keys():
            for cnot in self.overlaps[(i, j)]:
                self._cnot_order_constraint_x[(i, j)][cnot] = z3.Bool(f'xltz_{i}_{j}_{cnot}')
                self._cnot_order_constraint_z[(j, i)][cnot] = z3.Bool(f'zltx_{j}_{i}_{cnot}')
                self.solver.add(self._cnot_order_constraint_x[(i, j)][cnot] ==
                                      self._cnot_comes_before_cnot(self.x_vars[i], self.z_vars[j], cnot))
                self.solver.add(self._cnot_order_constraint_z[(j, i)][cnot] ==
                                      self._cnot_comes_before_cnot(self.z_vars[j], self.x_vars[i], cnot))

        for check in self.x_vars:
            self._assert_cnot_on_check_constraint(check)
            for cnot in check.keys():
                self.solver.add(self._cnot_exists_constraint(check, cnot))
        for check in self.z_vars:
            self._assert_cnot_on_check_constraint(check)
            for cnot in check.keys():
                self.solver.add(self._cnot_exists_constraint(check, cnot))

        # assert that non-commuting cnots are ordered correctly
        for i, j in self.overlaps.keys():
            self._assert_check_order_constraint(i, j)

        # assert that each qubit is used at most once per timestep
        for i, qubit in enumerate(self.qubit_variables):
            for t in range(self.T):
                self.solver.add(self._cnots_not_overlapping_constraint(i, t))

    def solve(self):
        """Solve the problem."""
        self._encode_constraints()

        print("finished encoding constraints")
        result = self.solver.check()
        if result != z3.sat:
            return result

        self._circuit = self._extract_circuit()
        return str(result)

    def get_circuit(self):
        """Return the circuit."""
        if self._circuit is None:
            raise ValueError("No circuit available. Call solve() first.")
        return self._circuit

    def _extract_circuit(self):
        q = QuantumRegister(len(self.qubit_variables), 'q')
        x_c = ClassicalRegister(self.n_xchecks, 'x_c')
        z_c = ClassicalRegister(self.n_zchecks, 'z_c')
        x_anc = AncillaRegister(self.n_xchecks, 'x')
        z_anc = AncillaRegister(self.n_zchecks, 'z')
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
