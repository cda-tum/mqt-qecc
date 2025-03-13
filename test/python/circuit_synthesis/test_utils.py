"""Test utility functions for the circuit synthesis module."""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pytest
from ldpc import mod2
from qiskit import AncillaRegister, ClassicalRegister, QuantumCircuit, QuantumRegister
from stim import Flow, PauliString

from mqt.qecc.circuit_synthesis.state_prep import final_matrix_constraint
from mqt.qecc.circuit_synthesis.synthesis_utils import (
    gaussian_elimination_min_column_ops,
    gaussian_elimination_min_parallel_eliminations,
    measure_flagged,
    measure_one_flagged,
    measure_stab_unflagged,
    measure_two_flagged,
    qiskit_to_stim_circuit,
)

if TYPE_CHECKING:
    import numpy.typing as npt
    from qiskit import AncillaQubit, ClBit, Qubit


class MatrixTest(NamedTuple):
    """Test matrix and expected results."""

    matrix: npt.NDArray[np.int8]
    min_col_ops: int
    max_parallel_steps: int


def check_correct_elimination(
    matrix: npt.NDArray[np.int8], final_matrix: npt.NDArray[np.int8], col_ops: list[tuple[int, int]]
) -> bool:
    """Check if the matrix is correctly eliminated."""
    matrix = matrix.copy()
    for op in col_ops:
        matrix[:, op[1]] = (matrix[:, op[0]] + matrix[:, op[1]]) % 2
    return np.array_equal(matrix, final_matrix)


def get_n_parallel_layers(ops: list[tuple[int, int]]) -> int:
    """Get the number of parallel layers in the elimination."""
    used_cols: set[int] = set()
    layer = 0
    for op in ops:
        if op[0] in used_cols or op[1] in used_cols:
            layer += 1
            used_cols = set()
        used_cols.add(op[0])
        used_cols.add(op[1])
    return layer


@pytest.fixture
def identity_matrix() -> MatrixTest:
    """Return a 4x4 identity matrix."""
    return MatrixTest(np.eye(4, dtype=np.int8), 0, 0)


@pytest.fixture
def full_matrix() -> MatrixTest:
    """Return a 4x4 matrix with all ones."""
    return MatrixTest(np.array([[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]], dtype=np.int8), 3, 2)


class MeasurementTest(NamedTuple):
    """Class containing all information for a measurement test."""

    qc: QuantumCircuit
    stab: list[Qubit]
    ancilla: AncillaQubit
    measurement_bit: ClBit
    t: int
    z_measurement: bool


def _make_measurement_test(n: int, stab: list[int], t: int, z_measurements: bool) -> MeasurementTest:
    q = QuantumRegister(n, "q")
    c = ClassicalRegister(1, "c")
    anc = AncillaRegister(1, "anc")
    qc = QuantumCircuit(q, c, anc)
    stab_qubits = [q[i] for i in stab]
    ancilla = anc[0]
    measurement_bit = c[0]
    return MeasurementTest(qc, stab_qubits, ancilla, measurement_bit, t, z_measurements)


@pytest.fixture
def weight_4_z() -> MeasurementTest:
    """Return a measurement test for a weight 4 stabilizer."""
    return _make_measurement_test(4, [0, 1, 2, 3], 1, True)


@pytest.fixture
def weight_4_x() -> MeasurementTest:
    """Return a measurement test for a weight 4 stabilizer."""
    return _make_measurement_test(4, [3, 2, 1, 0], 1, False)


@pytest.fixture
def weight_4_z_6q() -> MeasurementTest:
    """Return a measurement test for a weight 4 stabilizer."""
    return _make_measurement_test(6, [0, 1, 2, 3], 1, True)


@pytest.fixture
def weight_5_z_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 6 stabilizer."""
    return _make_measurement_test(5, [0, 1, 2, 3, 4], 1, True)


@pytest.fixture
def weight_5_x_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 6 stabilizer."""
    return _make_measurement_test(5, [0, 1, 2, 3, 4], 1, False)


@pytest.fixture
def weight_6_z_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 6 stabilizer."""
    return _make_measurement_test(6, [0, 1, 2, 3, 4, 5], 1, True)


@pytest.fixture
def weight_6_x_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 6 stabilizer."""
    return _make_measurement_test(6, [0, 1, 2, 3, 4, 5], 1, False)


@pytest.fixture
def weight_7_z_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 8 stabilizer."""
    return _make_measurement_test(7, [6, 5, 4, 3, 2, 1, 0], 1, True)


@pytest.fixture
def weight_7_x_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 8 stabilizer."""
    return _make_measurement_test(7, [6, 5, 4, 3, 2, 1, 0], 1, False)


@pytest.fixture
def weight_8_x_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 8 stabilizer."""
    return _make_measurement_test(8, [7, 6, 5, 4, 3, 2, 1, 0], 1, False)


@pytest.fixture
def weight_8_z_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 8 stabilizer."""
    return _make_measurement_test(8, [7, 6, 5, 4, 3, 2, 1, 0], 1, True)


@pytest.fixture
def weight_9_z_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 8 stabilizer."""
    return _make_measurement_test(100, [8, 7, 6, 5, 4, 3, 2, 1, 0], 2, True)


@pytest.fixture
def weight_11_z_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 8 stabilizer."""
    return _make_measurement_test(11, list(range(11)), 2, True)


@pytest.fixture
def weight_11_x_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 8 stabilizer."""
    return _make_measurement_test(11, list(range(11)), 2, False)


@pytest.fixture
def weight_12_z_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 8 stabilizer."""
    return _make_measurement_test(12, list(range(12)), 2, True)


@pytest.fixture
def weight_12_x_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 8 stabilizer."""
    return _make_measurement_test(12, list(range(12)), 2, False)


@pytest.fixture
def weight_16_x_measurement() -> MeasurementTest:
    """Return a measurement test for a weight 16 stabilizer."""
    return _make_measurement_test(16, [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0], 2, False)


@pytest.mark.parametrize("test_vals", ["identity_matrix", "full_matrix"])
def test_min_column_ops(test_vals: MatrixTest, request) -> None:  # type: ignore[no-untyped-def]
    """Check correct number of column operations are returned."""
    fixture = request.getfixturevalue(test_vals)
    matrix = fixture.matrix
    min_col_ops = fixture.min_col_ops
    rank = mod2.rank(matrix)
    res = gaussian_elimination_min_column_ops(
        matrix, lambda checks: final_matrix_constraint(checks, rank), max_eliminations=fixture.min_col_ops
    )
    assert res is not None
    reduced, ops = res
    assert len(ops) == min_col_ops
    assert check_correct_elimination(matrix, reduced, ops)


@pytest.mark.parametrize("test_vals", ["identity_matrix", "full_matrix"])
def test_min_parallel_eliminations(test_vals: MatrixTest, request) -> None:  # type: ignore[no-untyped-def]
    """Check correct number of parallel eliminations are returned."""
    fixture = request.getfixturevalue(test_vals)
    matrix = fixture.matrix
    rank = mod2.rank(matrix)
    max_parallel_steps = fixture.max_parallel_steps
    res = gaussian_elimination_min_parallel_eliminations(
        matrix, lambda checks: final_matrix_constraint(checks, rank), max_parallel_steps=fixture.max_parallel_steps
    )
    assert res is not None
    reduced, ops = res

    assert check_correct_elimination(matrix, reduced, ops)

    n_parallel_layers = get_n_parallel_layers(ops)
    assert n_parallel_layers <= max_parallel_steps


def correct_stabilizer_propagation(
    qc: QuantumCircuit, stab: list[Qubit], ancilla: AncillaQubit, z_measurement: bool
) -> bool:
    """Check that the stabilizer is propagated correctly."""
    circ = qiskit_to_stim_circuit(qc)
    pauli = "Z" if z_measurement else "X"
    anc_idx = qc.find_bit(ancilla).index
    initial_pauli = PauliString("_" * (anc_idx) + "Z" + "_" * (len(qc.qubits) - anc_idx - 1))
    final_pauli = PauliString(
        "".join([pauli if q in stab else "_" for q in qc.qubits]) + "_" * (len(qc.qubits) - anc_idx - 1)
    )
    f = Flow(input=initial_pauli, output=final_pauli, measurements=[qc.num_ancillas - 1])
    return bool(circ.has_flow(f))


@pytest.mark.parametrize(
    "test", ["weight_4_z", "weight_4_z_6q", "weight_4_x", "weight_6_z_measurement", "weight_8_x_measurement"]
)
def test_one_flag_measurements(test: MeasurementTest, request):  # type: ignore[no-untyped-def]
    """Test measurement circuits."""
    fixture = request.getfixturevalue(test)
    qc = fixture.qc
    stab = fixture.stab
    ancilla = fixture.ancilla
    measurement_bit = fixture.measurement_bit
    z_measurement = fixture.z_measurement

    measure_one_flagged(qc, stab, ancilla, measurement_bit, z_measurement)
    assert qc.depth() == len(stab) + 3 + 2 * int(not z_measurement)  # 6 CNOTs + Measurement + 2 possible hadamards
    assert qc.count_ops().get("cx", 0) == len(stab) + 2  # CNOTs from measurement + 2 flagging CNOTs
    assert correct_stabilizer_propagation(qc, stab, ancilla, z_measurement)


@pytest.mark.parametrize(
    "test",
    [
        "weight_4_z",
        "weight_4_z_6q",
        "weight_4_x",
        "weight_5_z_measurement",
        "weight_5_x_measurement",
        "weight_6_z_measurement",
        "weight_6_x_measurement",
        "weight_7_z_measurement",
        "weight_7_x_measurement",
        "weight_8_x_measurement",
        "weight_9_z_measurement",
        "weight_11_z_measurement",
        "weight_11_x_measurement",
        "weight_12_z_measurement",
        "weight_12_x_measurement",
        "weight_16_x_measurement",
    ],
)
def test_two_flag(test: MeasurementTest, request):  # type: ignore[no-untyped-def]
    """Test two flag measurement circuits."""
    fixture = request.getfixturevalue(test)
    qc = fixture.qc
    stab = fixture.stab
    ancilla = fixture.ancilla
    measurement_bit = fixture.measurement_bit
    z_measurement = fixture.z_measurement

    measure_two_flagged(qc, stab, ancilla, measurement_bit, z_measurement)
    assert correct_stabilizer_propagation(qc, stab, ancilla, z_measurement)


@pytest.mark.parametrize(
    "test",
    [
        "weight_4_z",
        "weight_4_z_6q",
        "weight_4_x",
        "weight_6_z_measurement",
        "weight_8_x_measurement",
        "weight_9_z_measurement",
    ],
)
def test_flag_measurements(test: MeasurementTest, request) -> None:  # type: ignore[no-untyped-def]
    """Test flag measurement circuits."""
    fixture = request.getfixturevalue(test)
    qc = fixture.qc
    stab = fixture.stab
    ancilla = fixture.ancilla
    t = fixture.t
    measurement_bit = fixture.measurement_bit
    z_measurement = fixture.z_measurement

    measure_flagged(qc, stab, ancilla, measurement_bit, t, z_measurement)
    assert correct_stabilizer_propagation(qc, stab, ancilla, z_measurement)


@pytest.mark.parametrize(
    "test",
    [
        "weight_4_z",
        "weight_4_z_6q",
        "weight_4_x",
        "weight_6_z_measurement",
        "weight_8_x_measurement",
        "weight_9_z_measurement",
    ],
)
def test_unflagged_measurements(test: MeasurementTest, request) -> None:  # type: ignore[no-untyped-def]
    """Test flag measurement circuits."""
    fixture = request.getfixturevalue(test)
    qc = fixture.qc
    stab = fixture.stab
    ancilla = fixture.ancilla
    measurement_bit = fixture.measurement_bit
    z_measurement = fixture.z_measurement

    measure_stab_unflagged(qc, stab, ancilla, measurement_bit, z_measurement)
    assert correct_stabilizer_propagation(qc, stab, ancilla, z_measurement)
