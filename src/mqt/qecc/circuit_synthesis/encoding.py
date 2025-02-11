"""Methods for synthesizing encoding circuits for CSS codes."""

from __future__ import annotations

import functools
import logging
import operator
from typing import TYPE_CHECKING

import numpy as np
import z3
from ldpc import mod2
from qiskit import QuantumCircuit

from ..codes import InvalidCSSCodeError
from .synthesis_utils import (
    build_css_circuit_from_cnot_list,
    heuristic_gaussian_elimination,
    optimal_elimination,
    symbolic_vector_add,
    symbolic_vector_eq,
)

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt

    from ..codes import CSSCode, StabilizerCode

logger = logging.getLogger(__name__)


def heuristic_encoding_circuit(
    code: CSSCode, optimize_depth: bool = True, balance_checks: bool = False
) -> QuantumCircuit:
    """Synthesize an encoding circuit for the given CSS code using a heuristic greedy search.

    Args:
        code: The CSS code to synthesize the encoding circuit for.
        optimize_depth: Whether to optimize the depth of the circuit.
        balance_checks: Whether to balance the entries of the stabilizer matrix via row operations.

    Returns:
        The synthesized encoding circuit and the qubits that are used to encode the logical qubits.
    """
    logger.info("Starting encoding circuit synthesis.")

    checks, logicals, use_x_checks = _get_matrix_with_fewest_checks(code)
    n_checks = checks.shape[0]

    if balance_checks:
        _balance_matrix(logicals)

    checks, cnots = heuristic_gaussian_elimination(
        np.vstack((checks, logicals)),
        parallel_elimination=optimize_depth,
    )

    # after reduction there still might be some overlap between initialized qubits and encoding qubits, we simply perform CNOTs to correct this
    encoding_qubits = np.where(checks[n_checks:, :].sum(axis=0) != 0)[0]
    initialization_qubits = np.where(checks[:n_checks, :].sum(axis=0) != 0)[0]
    # remove encoding qubits from initialization qubits
    initialization_qubits = np.setdiff1d(initialization_qubits, encoding_qubits)
    rows = []  # type: list[int]
    qubit_to_row = {}
    for qubit in initialization_qubits:
        cand = np.where(checks[:n_checks, qubit] == 1)[0]
        np.setdiff1d(cand, np.array(rows))
        rows.append(cand[0])
        qubit_to_row[qubit] = cand[0]

    for init_qubit in initialization_qubits:
        for encoding_qubit in encoding_qubits:
            row = qubit_to_row[init_qubit]
            if checks[row, encoding_qubit] == 1:
                cnots.append((init_qubit, encoding_qubit))
                checks[row, encoding_qubit] = 0

    cnots = cnots[::-1]

    return _build_css_encoder_from_cnot_list(n_checks, checks, cnots, use_x_checks)


def depth_optimal_encoding_circuit_non_css(
    code: StabilizerCode,
    min_depth: int = 1,
    max_depth: int = 10,
    min_timeout: int = 1,
    max_timeout: int = 3600,
) -> QuantumCircuit:
    n = code.n
    stabs = code.symplectic
    z_logicals = code.z_logicals.tableau.matrix
    x_logicals = code.x_logicals.tableau.matrix
    bit_width = 3

    # Gate constants
    IDENTITY = z3.BitVecVal(0, bit_width)
    HADAMARD = z3.BitVecVal(1, bit_width)
    SGATE = z3.BitVecVal(2, bit_width)
    SQRTX = z3.BitVecVal(3, bit_width)
    CXCTRL = z3.BitVecVal(4, bit_width)
    CXTAR = z3.BitVecVal(5, bit_width)
    CZ = z3.BitVecVal(6, bit_width)
    CZ2 = z3.BitVecVal(7, bit_width)

    # create variables
    sqgs = [
        [z3.BitVec(f"sqg_{t}_{q}", 3) for q in range(n)] for t in range(max_depth)
    ]  # I=0, H=1, S=2, SqrtX = 3, CXCtrl=4, CXTar=5, CZ=6|7
    czs = [[{q2: z3.Bool(f"cz_{t}_{q1}_{q2}") for q2 in range(q1 + 1, n)} for q1 in range(n)] for t in range(max_depth)]
    cxs = [[[z3.Bool(f"cx_{t}_{q1}_{q2}") for q2 in range(n)] for q1 in range(n)] for t in range(max_depth)]

    stab_x = [
        np.array([[z3.Bool(f"stab_x_{t}_{stab}_{q}") for q in range(n)] for stab in range(code.n - code.k)])
        for t in range(max_depth + 1)
    ]
    stab_z = [
        np.array([[z3.Bool(f"stab_z_{t}_{stab}_{q}") for q in range(n)] for stab in range(code.n - code.k)])
        for t in range(max_depth + 1)
    ]
    log_x_x = [
        np.array([[z3.Bool(f"log_x_x_{t}_{stab}_{q}") for q in range(n)] for stab in range(code.k)])
        for t in range(max_depth + 1)
    ]
    log_x_z = [
        np.array([[z3.Bool(f"log_x_z_{t}_{stab}_{q}") for q in range(n)] for stab in range(code.k)])
        for t in range(max_depth + 1)
    ]
    log_z_x = [
        np.array([[z3.Bool(f"log_z_x_{t}_{stab}_{q}") for q in range(n)] for stab in range(code.k)])
        for t in range(max_depth + 1)
    ]
    log_z_z = [
        np.array([[z3.Bool(f"log_z_z_{t}_{stab}_{q}") for q in range(n)] for stab in range(code.k)])
        for t in range(max_depth + 1)
    ]

    stab_x[0][:, :] = stabs[:, :n].astype(bool)
    stab_z[0][:, :] = stabs[:, n:].astype(bool)
    log_x_x[0][:, :] = x_logicals[:, :n].astype(bool)
    log_x_z[0][:, :] = x_logicals[:, n:].astype(bool)
    log_z_x[0][:, :] = z_logicals[:, :n].astype(bool)
    log_z_z[0][:, :] = z_logicals[:, n:].astype(bool)

    s = z3.Solver()
    # consistency
    for t in range(max_depth):
        for q1 in range(n):
            tqgs_q1_cz = list(czs[t][q1].values()) + [czs[t][q2][q1] for q2 in range(q1)]
            tqgs_q1_ctrl = cxs[t][q1]
            tqgs_q1_tar = [cxs[t][q2][q1] for q2 in range(n)]
            is_ctrl = z3.Or(tqgs_q1_ctrl)
            is_tar = z3.Or(tqgs_q1_tar)
            tqgs_q1_ctrl + tqgs_q1_tar
            tqgs_q1_ctrl + tqgs_q1_tar + tqgs_q1_cz
            # # if a single_qubit_gate is applied, then no two_qubit_gate can be applied
            # s.add(z3.Implies(z3.Or(tqgs_q1), sqgs[t][q1] == z3.BitVecVal(0, 2)))

            # if a gate is not a single-qubit gate, then the appropriate variables must be set
            s.add(z3.Implies(sqgs[t][q1] == CXCTRL, z3.Or(tqgs_q1_ctrl)))
            s.add(z3.Implies(sqgs[t][q1] == CXTAR, z3.Or(tqgs_q1_tar)))
            s.add(z3.Implies(sqgs[t][q1] == CZ | sqgs[t][q1] == CZ2, z3.Or(tqgs_q1_cz)))
            # if a two-qubit gate is applied it must either be a CZ or a CX
            # tqgs_q1 = z3.Not(z3.And(z3.Or(tqgs_q1_cz), z3.Or(is_tar, is_ctrl)))
            # a qubit can be either control or target, not both
            s.add(z3.Not(z3.And(is_ctrl, is_tar)))
            # a qubit can be the control of at most one CNOT
            s.add(z3.PbLe([(cxs[t][q1][q2], 1) for q2 in range(n)], 1))
            # a qubit can be the target of at most one CNOT
            s.add(z3.PbLe([(cxs[t][q2][q1], 1) for q2 in range(n)], 1))
            # a qubit can be involved in at most one CZ
            s.add(z3.PbLe([(cz, 1) for cz in tqgs_q1_cz], 1))

    # gate constraints
    tableaus_x = [np.vstack((stab_x[t], log_x_x[t], log_z_x[t])) for t in range(max_depth + 1)]
    tableaus_z = [np.vstack((stab_z[t], log_x_z[t], log_z_z[t])) for t in range(max_depth + 1)]
    for t in range(1, max_depth + 1):
        for q1 in range(n):
            # single-qubit gates
            s.add(
                z3.Implies(
                    sqgs[t - 1][q1] == IDENTITY,
                    z3.And(
                        symbolic_vector_eq(tableaus_x[t][:, q1], tableaus_x[t - 1][:, q1]),
                        symbolic_vector_eq(tableaus_z[t][:, q1], tableaus_z[t - 1][:, q1]),
                    ),
                )
            )  # I
            s.add(
                z3.Implies(
                    sqgs[t - 1][q1] == HADAMARD,
                    z3.And(
                        symbolic_vector_eq(tableaus_x[t][:, q1], tableaus_z[t - 1][:, q1]),
                        symbolic_vector_eq(tableaus_z[t][:, q1], tableaus_x[t - 1][:, q1]),
                    ),
                )
            )  # H
            s.add(
                z3.Implies(
                    sqgs[t - 1][q1] == SGATE,
                    z3.And(
                        symbolic_vector_eq(
                            symbolic_vector_add(tableaus_x[t - 1][:, q1], tableaus_z[t - 1][:, q1]),
                            tableaus_z[t][:, q1],
                        ),
                        symbolic_vector_eq(tableaus_x[t][:, q1], tableaus_x[t - 1][:, q1]),
                    ),
                )
            )  # S
            s.add(
                z3.Implies(
                    sqgs[t - 1][q1] == SQRTX,
                    z3.And(
                        symbolic_vector_eq(
                            symbolic_vector_add(tableaus_x[t - 1][:, q1], tableaus_z[t - 1][:, q1]),
                            tableaus_x[t][:, q1],
                        ),
                        symbolic_vector_eq(tableaus_z[t][:, q1], tableaus_z[t - 1][:, q1]),
                    ),
                )
            )  # SqrtX

            # two-qubit gates
            for q2 in range(n):
                s.add(
                    z3.Implies(
                        cxs[t - 1][q1][q2],
                        z3.And(
                            symbolic_vector_eq(tableaus_x[t][:, q1], tableaus_x[t - 1][:, q1]),
                            symbolic_vector_eq(
                                tableaus_x[t][:, q2],
                                symbolic_vector_add(tableaus_x[t - 1][:, q1], tableaus_x[t - 1][:, q2]),
                            ),
                            symbolic_vector_eq(tableaus_z[t][:, q2], tableaus_z[t - 1][:, q2]),
                            symbolic_vector_eq(
                                tableaus_z[t][:, q1],
                                symbolic_vector_add(tableaus_z[t - 1][:, q1], tableaus_z[t - 1][:, q2]),
                            ),
                        ),
                    )
                )  # CX

            for q2 in range(q1 + 1, n):
                s.add(
                    z3.Implies(
                        czs[t - 1][q1][q2],
                        z3.And(
                            symbolic_vector_eq(tableaus_x[t][:, q1], tableaus_x[t - 1][:, q1]),
                            symbolic_vector_eq(tableaus_x[t][:, q2], tableaus_x[t - 1][:, q2]),
                            symbolic_vector_eq(
                                tableaus_z[t][:, q1],
                                symbolic_vector_add(tableaus_z[t - 1][:, q1], tableaus_x[t - 1][:, q2]),
                            ),
                            symbolic_vector_eq(
                                tableaus_z[t][:, q2],
                                symbolic_vector_add(tableaus_x[t - 1][:, q1], tableaus_z[t - 1][:, q2]),
                            ),
                        ),
                    )
                )  # CZ
    # final matrix constraints

    # for the stabilizers, we only require that the final tableau is full rank
    # qubit_reduced = [z3.Not(z3.Or(
    #     z3.Or(
    #         [stab_x[-1][i, q] for i in range(code.n-code.k)]
    #     ),
    #     z3.Or(
    #         [stab_z[-1][i, q] for i in range(code.n-code.k)]
    #     )
    # ))
    #                  for q in range(n)
    #                  ]
    # s.add(z3.PbEq([(qubit_reduced[q], 1) for q in range(n)], code.k))
    s.add(
        z3.PbEq(
            [(z3.Or([stab_z[-1][i, q] for i in range(code.n - code.k)]), 1) for q in range(n)]
            + [(z3.Or([stab_x[-1][i, q] for i in range(code.n - code.k)]), 1) for q in range(n)],
            code.n - code.k,
        )
    )

    # s.add(z3.Not(z3.Or([stab_x[-1][i, q] for i in range(code.n-code.k) for q in range(n)])))

    # for the logicals we require that the final tableau is full rank and that the logicals are orthogonal to the stabilizers
    # this means that every row has exactly one x and one z entry for the same qubit
    s.add([
        z3.And(
            z3.PbEq([(x, 1) for x in log_x_x[-1][i]], 1),
            z3.PbEq([(z, 1) for z in log_z_z[-1][i]], 1),
            z3.Not(z3.Or(list(log_x_z[-1][i]) + list(log_z_x[-1][i]))),
        )
        for i in range(code.k)
    ])

    # solve
    if s.check() == z3.sat:
        m = s.model()
        # print last tableau
        for t in range(1, max_depth + 1):
            for _i in range(n - code.k):
                pass
            for _i in range(code.k):
                pass
        # extract circuit
        qc = QuantumCircuit(n)
        final_tableau_x = np.array([
            [bool(m[stab_x[-1][i, q]]) for q in range(n)] for i in range(code.n - code.k)
        ]).astype(int)
        [[m[stab_z[-1][i, q]] for i in range(code.n - code.k)] for q in range(n)]
        # check where hadamards need to be applied
        for t in range(max_depth):
            for q1 in range(n):
                if m[sqgs[t][q1]] == 1:
                    qc.h(q1)
                elif m[sqgs[t][q1]] == 2:
                    qc.sdg(q1)
                elif m[sqgs[t][q1]] == 3:
                    qc.sxdg(q1)
                for q2 in range(n):
                    if m[cxs[t][q1][q2]]:
                        qc.cx(q1, q2)
                    if q2 > q1 and m[czs[t][q1][q2]]:
                        qc.cz(q1, q2)
            qc.barrier()

        first_layer_hadamards = np.where(np.array(final_tableau_x).sum(axis=0) == 1)[0]
        if len(first_layer_hadamards) > 0:
            qc.h(first_layer_hadamards)
        return qc.inverse()
    return "UNSAT"


def gate_optimal_encoding_circuit(
    code: CSSCode,
    min_gates: int = 1,
    max_gates: int = 10,
    min_timeout: int = 1,
    max_timeout: int = 3600,
) -> QuantumCircuit:
    """Synthesize an encoding circuit for the given CSS code using the minimal number of gates.

    Args:
        code: The CSS code to synthesize the encoding circuit for.
        min_gates: The minimum number of gates to use in the circuit.
        max_gates: The maximum number of gates to use in the circuit.
        min_timeout: The minimum time to spend on the synthesis.
        max_timeout: The maximum time to spend on the synthesis.

    Returns:
        The synthesized encoding circuit and the qubits that are used to encode the logical qubits.
    """
    logger.info("Starting optimal encoding circuit synthesis.")
    checks, logicals, use_x_checks = _get_matrix_with_fewest_checks(code)
    assert checks is not None
    n_checks = checks.shape[0]
    checks_and_logicals = np.vstack((checks, logicals))
    rank = mod2.rank(checks_and_logicals)
    termination_criteria = functools.partial(
        _final_matrix_constraint_partially_full_reduction,
        full_reduction_rows=list(range(checks.shape[0], checks.shape[0] + logicals.shape[0])),
        rank=rank,
    )

    res = optimal_elimination(
        checks_and_logicals,
        termination_criteria,
        "column_ops",
        min_param=min_gates,
        max_param=max_gates,
        min_timeout=min_timeout,
        max_timeout=max_timeout,
    )
    if res is None:
        return None
    reduced_checks_and_logicals, cnots = res
    cnots = cnots[::-1]

    return _build_css_encoder_from_cnot_list(n_checks, reduced_checks_and_logicals, cnots, use_x_checks)


def depth_optimal_encoding_circuit(
    code: CSSCode,
    min_depth: int = 1,
    max_depth: int = 10,
    min_timeout: int = 1,
    max_timeout: int = 3600,
) -> QuantumCircuit:
    """Synthesize an encoding circuit for the given CSS code using minimal depth.

    Args:
        code: The CSS code to synthesize the encoding circuit for.
        min_depth: The minimum number of gates to use in the circuit.
        max_depth: The maximum number of gates to use in the circuit.
        min_timeout: The minimum time to spend on the synthesis.
        max_timeout: The maximum time to spend on the synthesis.

    Returns:
        The synthesized encoding circuit and the qubits that are used to encode the logical qubits.
    """
    logger.info("Starting optimal encoding circuit synthesis.")
    checks, logicals, use_x_checks = _get_matrix_with_fewest_checks(code)
    assert checks is not None
    n_checks = checks.shape[0]
    checks_and_logicals = np.vstack((checks, logicals))
    rank = mod2.rank(checks_and_logicals)
    termination_criteria = functools.partial(
        _final_matrix_constraint_partially_full_reduction,
        full_reduction_rows=list(range(checks.shape[0], checks.shape[0] + logicals.shape[0])),
        rank=rank,
    )
    res = optimal_elimination(
        checks_and_logicals,
        termination_criteria,
        "parallel_ops",
        min_param=min_depth,
        max_param=max_depth,
        min_timeout=min_timeout,
        max_timeout=max_timeout,
    )
    if res is None:
        return None
    reduced_checks_and_logicals, cnots = res
    cnots = cnots[::-1]

    return _build_css_encoder_from_cnot_list(n_checks, reduced_checks_and_logicals, cnots, use_x_checks)


# def gottesmann_encoding_circuit(code: StabilizerCode) -> QuantumCircuit:
#     """Construct an encoding circuit for the given stabilizer code using the scheme in Gottesmann's book.

#     Args:
#             code: The stabilizer code to construct the encoding circuit for.

#     Returns:
#             The encoding circuit.
#     """
#     qc = QuantumCircuit(code.n)
#     matrix = code.symplectic.copy()

#     def move_to_diagonal(matrix:npt.NDArray[np.int], row: int, col: int) -> None:
#         if col < code.n:
#             qc.swap(row, col)
#             matrix[[row, col]] = matrix[[col, row]]
#         else:
#             qc.h(row)
#             qc.swap(row, col - code.n)
#             matrix[[row, col]] = matrix[[col, row]]
#             matrix[row, :] = (matrix[row, :] + matrix[col, :]) % 2
#     for row in range(code.n-code.k):
#         # find first non-zero entry in row
#         col = np.where(matrix[row, :])[0][0]
#         if row != col:        # move to diagonal
#             move_to_diagonal(matrix, row, col)
#         # reduce row
#         non_zero_x = np.where(matrix[row, :code.n])[0][1:]
#         non_zero_z = np.where(matrix[row, code.n:])[0]
#         if row in non_zero_z:
#             non_zero_z.remove(row)
#             qc.s(row)
#         qc.cx(row, non_zero_x)
#         qc.cz(row, non_zero_z)
#         matrix[row, :] = 0
#         matrix[:, row] = 0 # reduce columns (change stabilizer generators)
#         matrix[row, row] = 1 # reset the 1 entry

#     # perform final hadamards
#     qc.h(range(code.n-code.k))

#     # correct sign
#     tableau = StabilizerTableau.identity(code.n)
#     updated_tableau = apply_clifford_circuit(tableau, qc)
#     corrections = np.where(updated_tableau.phase == 1)[0]

#     qc = qc.inverse()
#     qc.x(corrections[:code.n-code.k])
#     tableau.phase = 0

#     if code.z_logicals is None:
#         return qc

#     # logicals are given, so compute the difference between the desired logicals and the actual logicals
#     x_tableau = StabilizerTableau.identity(code.n)
#     for i in range(code.n-code.k, code.n):
#         x_tableau.apply_x(i)
#     x_tableau = apply_clifford_circuit(x_tableau, qc)
#     z_matrix = tableau.tableau.matrix
#     x_matrix = x_tableau.tableau.matrix
#     z_diff = (code.z_logicals.tableau.matrix - z_matrix)%2
#     x_diff = (code.x_logicals.tableau.matrix - x_matrix)%2
#     matrix = np.hstack((z_diff, x_diff))

#     right_qc = QuantumCircuit(code.n)
#     left_qc = QuantumCircuit(code.n)
#     for row in range(code.n - code.k, code.n):
#         # find first non-zero entry in row
#         col = np.where(matrix[row, :])[0][0]
#         if row != col:        # move to diagonal
#             move_to_diagonal(matrix, row, col)
#         # reduce row
#         non_zero_x = np.where(matrix[row, :code.n])[0][1:]
#         non_zero_z = np.where(matrix[row, code.n:])[0]
#         if row in non_zero_z:
#             non_zero_z.remove(row)
#             right_qc.s(row)
#         right_qc.cx(row, non_zero_x)
#         right_qc.cz(row, non_zero_z)
#         matrix[row, :] = 0
#         matrix[row, row] = 1 # reset the 1 entry

#         # reduce column
#         non_zero_stab = np.where(matrix[row, :code.n])[0][1:]
#         left_qc.cx(row, non_zero_stab)
#         # qc.cz(row, non_zero_z)
#         matrix[:, row] = 0
#         matrix[row, row] = 1 # reset the 1 entry

#     # perform final hadamards
#     right_qc.h(range(code.n-code.k, code.n))

# def standard_encoding_circuit(code:StabilizerCode) -> QuantumCircuit:
#     """Construct an encoding circuit for the given stabilizer code using the standard method.

#     Args:
#             code: The stabilizer code to construct the encoding circuit for.

#     Returns:
#             The encoding circuit.
#     """
#     stabs = [str(stab) for stab in code.generators.to_pauli_list()]
#     if code.z_logicals is not None:
#         logicals = [str(logical) for logical in code.z_logicals.to_pauli_list()]
#         stabs += logicals
#     circ = StabilizerState.from_stabilizer_list(stabs).clifford.to_circuit()
#     message_qubits = list(range(code.n - code.k, code.n))
#     return


def _get_matrix_with_fewest_checks(code: CSSCode) -> tuple[npt.NDArray[np.int8], npt.NDArray[np.int8], bool]:
    """Return the stabilizer matrix with the fewest checks, the corresponding logicals and a bool indicating whether X- or Z-checks have been returned."""
    if code.Hx is None or code.Hz is None:
        msg = "The code must have both X and Z stabilizers defined."
        raise InvalidCSSCodeError(msg)

    use_x_checks = code.Hx.shape[0] < code.Hz.shape[0]
    checks = code.Hx if use_x_checks else code.Hz
    logicals = code.Lx if use_x_checks else code.Lz
    return checks, logicals, use_x_checks


def _final_matrix_constraint_partially_full_reduction(
    columns: npt.NDArray[z3.BoolRef | bool], full_reduction_rows: list[int], rank: int
) -> z3.BoolRef:
    assert len(columns.shape) == 3

    partial_reduction_rows = list(set(range(columns.shape[1])) - set(full_reduction_rows))

    # assert that the partial_reduction_rows are partially reduced, i.e. there are at least columns.shape[2] - (columns.shape[1] - len(full_reduction_rows)) non-zero columns
    partially_reduced = z3.PbEq(
        [(z3.Not(z3.Or(list(columns[-1, partial_reduction_rows, col]))), 1) for col in range(columns.shape[2])],
        columns.shape[2] - (rank - len(full_reduction_rows)),
    )

    # assert that there is no overlap between the full_reduction_rows and the partial_reduction_rows
    overlap_constraints = [True]
    for col in range(columns.shape[2]):
        has_entry_partial = z3.Or(list(columns[-1, partial_reduction_rows, col]))
        has_entry_full = z3.Or(list(columns[-1, full_reduction_rows, col]))
        overlap_constraints.append(z3.Not(z3.And(has_entry_partial, has_entry_full)))

    # assert that the full_reduction_rows are fully reduced
    fully_reduced = z3.PbEq(
        [
            (z3.PbEq([(columns[-1, row, col], 1) for col in range(columns.shape[2])], 1), 1)
            for row in full_reduction_rows
        ],
        len(full_reduction_rows),
    )

    return z3.And(fully_reduced, partially_reduced, z3.And(overlap_constraints))


def _build_css_encoder_from_cnot_list(
    n_checks: int, checks_and_logicals: npt.NDArray[np.int8], cnots: list[tuple[int, int]], use_x_checks: bool
) -> tuple[QuantumCircuit, list[int]]:
    encoding_qubits = np.where(checks_and_logicals[n_checks:, :].sum(axis=0) != 0)[0]
    if use_x_checks:
        hadamards = np.where(checks_and_logicals[:n_checks, :].sum(axis=0) != 0)[0]
    else:
        hadamards = np.where(checks_and_logicals[:n_checks, :].sum(axis=0) == 0)[0]
        cnots = [(j, i) for i, j in cnots]

    hadamards = np.setdiff1d(hadamards, encoding_qubits)
    circ = build_css_circuit_from_cnot_list(checks_and_logicals.shape[1], cnots, list(hadamards))
    return circ, encoding_qubits


def _balance_matrix(m: npt.NDArray[np.int8]) -> None:
    """Balance the columns of the matrix.

    Try to balance the number of 1's in each column via row operations without increasing the row-weight.
    """
    variance = np.var(m.sum(axis=0))
    reduced = False

    while not reduced:
        reduced = True
        # compute row operations that do not increase the row-weights
        row_ops = []
        for i, row_1 in enumerate(m):
            for j, row_2 in enumerate(m):
                if i == j:
                    continue
                s = (row_1 + row_2) % 2
                if s.sum() > row_1.sum() or s.sum() > row_2.sum():
                    continue
                # compute associated column weights
                m[j] = s  # noqa: B909

                new_variance = np.var(m.sum(axis=0))
                if new_variance < variance:
                    row_ops.append((i, j, new_variance))

                m[j] = row_2  # noqa: B909
        # sort by lowest variance
        row_ops.sort(key=operator.itemgetter(2))
        # apply best row operation
        if row_ops:
            i, j = row_ops[0][:2]
            m[i] = (m[i] + m[j]) % 2
            reduced = False
            variance = row_ops[0][2]
