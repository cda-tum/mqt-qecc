"""Test synthesis of state preparation and verification circuits for FT state preparation."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from ldpc import mod2
from qiskit.quantum_info import Clifford

from mqt.qecc import CSSCode
from mqt.qecc.codes import SquareOctagonColorCode
from mqt.qecc.ft_stateprep import (
    depth_optimal_prep_circuit,
    gate_optimal_prep_circuit,
    gate_optimal_verification_circuit,
    gate_optimal_verification_stabilizers,
    heuristic_prep_circuit,
    heuristic_verification_circuit,
    heuristic_verification_stabilizers,
)

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt
    from qiskit import QuantumCircuit

    from mqt.qecc.ft_stateprep import StatePrepCircuit


@pytest.fixture
def steane_code() -> CSSCode:
    """Return the Steane code."""
    return CSSCode.from_code_name("Steane")


@pytest.fixture
def css_4_2_2_code() -> CSSCode:
    """Return the 4,2,2  code."""
    return CSSCode(2, np.array([[1] * 4]), np.array([[1] * 4]))


@pytest.fixture
def css_6_2_2_code() -> CSSCode:
    """Return the 4,2,2  code."""
    return CSSCode(2, np.array([[1] * 6]), np.array([[1] * 6]))


@pytest.fixture
def surface_code() -> CSSCode:
    """Return the distance 3 rotated Surface Code."""
    return CSSCode.from_code_name("surface", 3)


@pytest.fixture
def tetrahedral_code() -> CSSCode:
    """Return the tetrahedral code."""
    return CSSCode.from_code_name("tetrahedral")


@pytest.fixture
def cc_4_8_8_code() -> CSSCode:
    """Return the d=5 4,8,8 color code."""
    return SquareOctagonColorCode(5)


@pytest.fixture
def steane_code_sp(steane_code: CSSCode) -> StatePrepCircuit:
    """Return a non-ft state preparation circuit for the Steane code."""
    sp_circ = heuristic_prep_circuit(steane_code)
    sp_circ.compute_fault_sets()
    return sp_circ


@pytest.fixture
def tetrahedral_code_sp(tetrahedral_code: CSSCode) -> StatePrepCircuit:
    """Return a non-ft state preparation circuit for the tetrahedral code."""
    sp_circ = heuristic_prep_circuit(tetrahedral_code)
    sp_circ.compute_fault_sets()
    return sp_circ


@pytest.fixture
def color_code_d5_sp(cc_4_8_8_code: CSSCode) -> StatePrepCircuit:
    """Return a non-ft state preparation circuit for the d=5 4,8,8 color code."""
    sp_circ = heuristic_prep_circuit(cc_4_8_8_code)
    sp_circ.compute_fault_sets()
    return sp_circ


def eq_span(a: npt.NDArray[np.int_], b: npt.NDArray[np.int_]) -> bool:
    """Check if two matrices have the same row space."""
    return a.shape == b.shape and mod2.rank(np.vstack((a, b))) == mod2.rank(a) == mod2.rank(b)


def in_span(m: npt.NDArray[np.int_], v: npt.NDArray[np.int_]) -> bool:
    """Check if a vector is in the row space of a matrix."""
    return bool(mod2.rank(np.vstack((m, v))) == mod2.rank(m))


def get_stabs(qc: QuantumCircuit) -> tuple[npt.NDArray[np.int_], npt.NDArray[np.int_]]:
    """Return the stabilizers of a quantum circuit."""
    cliff = Clifford(qc)
    x = cliff.stab_x.astype(int)
    x = x[np.where(np.logical_not(np.all(x == 0, axis=1)))[0]]
    z = cliff.stab_z.astype(int)
    z = z[np.where(np.logical_not(np.all(z == 0, axis=1)))[0]]
    return x, z


@pytest.mark.parametrize(
    "code", ["steane_code", "css_4_2_2_code", "css_6_2_2_code", "tetrahedral_code", "surface_code"]
)
def test_heuristic_prep_consistent(code: CSSCode, request) -> None:  # type: ignore[no-untyped-def]
    """Check that heuristic_prep_circuit returns a valid circuit with the correct stabilizers."""
    code = request.getfixturevalue(code)

    sp_circ = heuristic_prep_circuit(code)
    circ = sp_circ.circ
    max_cnots = np.sum(code.Hx) + np.sum(code.Hz)  # type: ignore[arg-type]

    assert circ.num_qubits == code.n
    assert circ.num_nonlocal_gates() <= max_cnots

    x, z = get_stabs(circ)
    assert eq_span(code.Hx, x)  # type: ignore[arg-type]
    assert eq_span(np.vstack((code.Hz, code.Lz)), z)  # type: ignore[arg-type]


@pytest.mark.parametrize("code", ["steane_code", "css_4_2_2_code", "css_6_2_2_code"])
def test_gate_optimal_prep_consistent(code: CSSCode, request) -> None:  # type: ignore[no-untyped-def]
    """Check that gate_optimal_prep_circuit returns a valid circuit with the correct stabilizers."""
    code = request.getfixturevalue(code)
    sp_circ = gate_optimal_prep_circuit(code, max_timeout=3)
    assert sp_circ is not None
    assert sp_circ.zero_state

    circ = sp_circ.circ
    max_cnots = np.sum(code.Hx) + np.sum(code.Hz)  # type: ignore[arg-type]

    assert circ.num_qubits == code.n
    assert circ.num_nonlocal_gates() <= max_cnots

    x, z = get_stabs(circ)
    assert eq_span(code.Hx, x)  # type: ignore[arg-type]
    assert eq_span(np.vstack((code.Hz, code.Lz)), z)  # type: ignore[arg-type]


@pytest.mark.parametrize("code", ["steane_code", "css_4_2_2_code", "css_6_2_2_code"])
def test_depth_optimal_prep_consistent(code: CSSCode, request) -> None:  # type: ignore[no-untyped-def]
    """Check that depth_optimal_prep_circuit returns a valid circuit with the correct stabilizers."""
    code = request.getfixturevalue(code)

    sp_circ = depth_optimal_prep_circuit(code, max_timeout=3)
    assert sp_circ is not None
    circ = sp_circ.circ
    max_cnots = np.sum(code.Hx) + np.sum(code.Hz)  # type: ignore[arg-type]

    assert circ.num_qubits == code.n
    assert circ.num_nonlocal_gates() <= max_cnots

    x, z = get_stabs(circ)
    assert eq_span(code.Hx, x)  # type: ignore[arg-type]
    assert eq_span(np.vstack((code.Hz, code.Lz)), z)  # type: ignore[arg-type]


@pytest.mark.parametrize("code", ["steane_code", "css_4_2_2_code", "css_6_2_2_code"])
def test_plus_state_gate_optimal(code: CSSCode, request) -> None:  # type: ignore[no-untyped-def]
    """Test synthesis of the plus state."""
    code = request.getfixturevalue(code)
    sp_circ_plus = gate_optimal_prep_circuit(code, max_timeout=3, zero_state=False)

    assert sp_circ_plus is not None
    assert not sp_circ_plus.zero_state

    circ_plus = sp_circ_plus.circ
    max_cnots = np.sum(code.Hx) + np.sum(code.Hz)  # type: ignore[arg-type]

    assert circ_plus.num_qubits == code.n
    assert circ_plus.num_nonlocal_gates() <= max_cnots

    x, z = get_stabs(circ_plus)
    assert eq_span(code.Hz, z)  # type: ignore[arg-type]
    assert eq_span(np.vstack((code.Hx, code.Lx)), x)  # type: ignore[arg-type]

    sp_circ_zero = gate_optimal_prep_circuit(code, max_timeout=5, zero_state=True)

    assert sp_circ_zero is not None

    circ_zero = sp_circ_zero.circ
    x_zero, z_zero = get_stabs(circ_zero)

    if code.is_self_dual():
        assert np.array_equal(x, z_zero)
        assert np.array_equal(z, x_zero)
    else:
        assert not np.array_equal(x, z_zero)
        assert not np.array_equal(z, x_zero)


@pytest.mark.parametrize(
    "code", ["steane_code", "css_4_2_2_code", "css_6_2_2_code", "surface_code", "tetrahedral_code"]
)
def test_plus_state_heuristic(code: CSSCode, request) -> None:  # type: ignore[no-untyped-def]
    """Test synthesis of the plus state."""
    code = request.getfixturevalue(code)
    sp_circ_plus = heuristic_prep_circuit(code, zero_state=False)

    assert sp_circ_plus is not None
    assert not sp_circ_plus.zero_state

    circ_plus = sp_circ_plus.circ
    max_cnots = np.sum(code.Hx) + np.sum(code.Hz)  # type: ignore[arg-type]

    assert circ_plus.num_qubits == code.n
    assert circ_plus.num_nonlocal_gates() <= max_cnots

    x, z = get_stabs(circ_plus)
    assert eq_span(code.Hz, z)  # type: ignore[arg-type]
    assert eq_span(np.vstack((code.Hx, code.Lx)), x)  # type: ignore[arg-type]

    sp_circ_zero = heuristic_prep_circuit(code, zero_state=True)
    circ_zero = sp_circ_zero.circ
    x_zero, z_zero = get_stabs(circ_zero)

    if code.is_self_dual():
        assert np.array_equal(x, z_zero)
        assert np.array_equal(z, x_zero)
    else:
        assert not np.array_equal(x, z_zero)
        assert not np.array_equal(z, x_zero)


def test_optimal_steane_verification_circuit(steane_code_sp: StatePrepCircuit) -> None:
    """Test that the optimal verification circuit for the Steane code is correct."""
    circ = steane_code_sp
    ver_stabs_layers = gate_optimal_verification_stabilizers(circ, x_errors=True, max_timeout=5)

    assert len(ver_stabs_layers) == 1  # 1 Ancilla measurement

    ver_stabs = ver_stabs_layers[0]

    assert np.sum(ver_stabs) == 3  # 3 CNOTs
    z_gens = circ.z_checks

    for stab in ver_stabs:
        assert in_span(z_gens, stab)

    errors = circ.compute_fault_set(1)
    non_detected = np.where(np.all(ver_stabs @ errors.T % 2 == 0, axis=1))[0]
    assert len(non_detected) == 0

    # Check that circuit is correct
    circ_ver = gate_optimal_verification_circuit(circ)
    assert circ_ver.num_qubits == circ.num_qubits + 1
    assert circ_ver.num_nonlocal_gates() == np.sum(ver_stabs) + circ.circ.num_nonlocal_gates()
    assert circ_ver.depth() == np.sum(ver_stabs) + circ.circ.depth() + 1  # 1 for the measurement


def test_heuristic_steane_verification_circuit(steane_code_sp: StatePrepCircuit) -> None:
    """Test that the optimal verification circuit for the Steane code is correct."""
    circ = steane_code_sp

    ver_stabs_layers = heuristic_verification_stabilizers(circ, x_errors=True)

    assert len(ver_stabs_layers) == 1  # 1 layer of verification measurements

    ver_stabs = ver_stabs_layers[0]
    assert len(ver_stabs) == 1  # 1 Ancilla measurement
    assert np.sum(ver_stabs[0]) == 3  # 3 CNOTs
    z_gens = circ.z_checks

    for stab in ver_stabs:
        assert in_span(z_gens, stab)

    errors = circ.compute_fault_set(1)
    non_detected = np.where(np.all(ver_stabs @ errors.T % 2 == 0, axis=1))[0]
    assert len(non_detected) == 0

    # Check that circuit is correct
    circ_ver = heuristic_verification_circuit(circ)
    assert circ_ver.num_qubits == circ.num_qubits + 1
    assert circ_ver.num_nonlocal_gates() == np.sum(ver_stabs) + circ.circ.num_nonlocal_gates()
    assert circ_ver.depth() == np.sum(ver_stabs) + circ.circ.depth() + 1  # 1 for the measurement


def test_optimal_tetrahedral_verification_circuit(tetrahedral_code_sp: StatePrepCircuit) -> None:
    """Test the optimal verification circuit for the tetrahedral code is correct.

    The tetrahedral code has an x-distance of 7. We expect that the verification only checks for a single propagated error since the tetrahedral code has a distance of 3.
    """
    circ = tetrahedral_code_sp

    ver_stabs_layers = gate_optimal_verification_stabilizers(circ, x_errors=True, max_ancillas=1, max_timeout=5)

    assert len(ver_stabs_layers) == 1  # 1 layer of verification measurements

    ver_stabs = ver_stabs_layers[0]
    assert len(ver_stabs) == 1  # 1 Ancilla measurement
    assert np.sum(ver_stabs[0]) == 3  # 3 CNOTs
    z_gens = circ.z_checks

    for stab in ver_stabs:
        assert in_span(z_gens, stab)

    errors = circ.compute_fault_set(1)
    non_detected = np.where(np.all(ver_stabs @ errors.T % 2 == 0, axis=1))[0]
    assert len(non_detected) == 0

    # Check that circuit is correct
    circ_ver = gate_optimal_verification_circuit(circ, max_ancillas=1, max_timeout=5)
    assert circ_ver.num_qubits == circ.num_qubits + 1
    assert circ_ver.num_nonlocal_gates() == np.sum(ver_stabs) + circ.circ.num_nonlocal_gates()
    assert circ_ver.depth() == np.sum(ver_stabs) + circ.circ.depth() + 1  # 1 for the measurement


def test_heuristic_tetrahedral_verification_circuit(tetrahedral_code_sp: StatePrepCircuit) -> None:
    """Test the optimal verification circuit for the tetrahedral code is correct.

    The tetrahedral code has an x-distance of 7. We expect that the verification only checks for a single propagated error since the tetrahedral code has a distance of 3.
    """
    circ = tetrahedral_code_sp

    ver_stabs_layers = heuristic_verification_stabilizers(circ, x_errors=True)

    assert len(ver_stabs_layers) == 1  # 1 layer of verification measurements

    ver_stabs = ver_stabs_layers[0]
    assert len(ver_stabs) == 1  # 1 Ancilla measurement
    assert np.sum(ver_stabs[0]) == 3  # 3 CNOTs
    z_gens = circ.z_checks

    for stab in ver_stabs:
        assert in_span(z_gens, stab)

    errors = circ.compute_fault_set(1)
    non_detected = np.where(np.all(ver_stabs @ errors.T % 2 == 0, axis=1))[0]
    assert len(non_detected) == 0

    # Check that circuit is correct
    circ_ver = heuristic_verification_circuit(circ)
    assert circ_ver.num_qubits == circ.num_qubits + 1
    assert circ_ver.num_nonlocal_gates() == np.sum(ver_stabs) + circ.circ.num_nonlocal_gates()
    assert circ_ver.depth() == np.sum(ver_stabs) + circ.circ.depth() + 1  # 1 for the measurement


def test_not_full_ft_opt_cc5(color_code_d5_sp: StatePrepCircuit) -> None:
    """Test that the optimal verification is also correct for higher distance.

    Ignore Z errors.
    Due to time constraints, we set the timeout for each search to 2 seconds.
    """
    circ = color_code_d5_sp

    ver_stabs_layers = gate_optimal_verification_stabilizers(circ, x_errors=True, max_ancillas=3, max_timeout=5)

    assert len(ver_stabs_layers) == 2  # 2 layers of verification measurements

    ver_stabs_1 = ver_stabs_layers[0]
    assert len(ver_stabs_1) == 2  # 2 Ancilla measurements
    assert np.sum(ver_stabs_1) == 9  # 9 CNOTs

    ver_stabs_2 = ver_stabs_layers[1]
    assert len(ver_stabs_2) == 3  # 2 Ancilla measurements
    assert np.sum(ver_stabs_2) <= 14  # less than 14 CNOTs (sometimes 13, sometimes 14 depending on how fast the CPU is)

    z_gens = circ.z_checks

    for stab in np.vstack((ver_stabs_1, ver_stabs_2)):
        assert in_span(z_gens, stab)

    errors_1 = circ.compute_fault_set(1)
    non_detected = np.where(np.all(ver_stabs_1 @ errors_1.T % 2 == 0, axis=1))[0]
    assert len(non_detected) == 0

    errors_2 = circ.compute_fault_set(2)
    non_detected = np.where(np.all(ver_stabs_2 @ errors_2.T % 2 == 0, axis=1))[0]
    assert len(non_detected) == 0

    # Check that circuit is correct
    n_cnots = np.sum(ver_stabs_1) + np.sum(ver_stabs_2)
    circ_ver = gate_optimal_verification_circuit(circ, max_ancillas=3, max_timeout=5, full_fault_tolerance=True)
    assert circ_ver.num_qubits > circ.num_qubits + 5  # overhead from the flags
    assert circ_ver.num_nonlocal_gates() > n_cnots + circ.circ.num_nonlocal_gates()  # Overhead from Flag CNOTS


def test_not_full_ft_heuristic_cc5(color_code_d5_sp: StatePrepCircuit) -> None:
    """Test that the optimal verification circuit for the Steane code is correct.

    Ignore Z errors.
    """
    circ = color_code_d5_sp
    ver_stabs_layers = heuristic_verification_stabilizers(circ, x_errors=True)

    assert len(ver_stabs_layers) == 2  # 2 layers of verification measurements

    ver_stabs_1 = ver_stabs_layers[0]
    ver_stabs_2 = ver_stabs_layers[1]

    z_gens = circ.z_checks

    for stab in np.vstack((ver_stabs_1, ver_stabs_2)):
        assert in_span(z_gens, stab)

    errors_1 = circ.compute_fault_set(1)
    non_detected = np.where(np.all(ver_stabs_1 @ errors_1.T % 2 == 0, axis=1))[0]
    assert len(non_detected) == 0

    errors_2 = circ.compute_fault_set(2)
    non_detected = np.where(np.all(ver_stabs_2 @ errors_2.T % 2 == 0, axis=1))[0]
    assert len(non_detected) == 0

    # Check that circuit is correct
    circ_ver = heuristic_verification_circuit(circ, full_fault_tolerance=False)
    n_cnots = np.sum(ver_stabs_1) + np.sum(ver_stabs_2)
    assert circ_ver.num_qubits == circ.num_qubits + len(ver_stabs_1) + len(ver_stabs_2)
    assert circ_ver.num_nonlocal_gates() == n_cnots + circ.circ.num_nonlocal_gates()


def test_full_ft_opt_cc5(color_code_d5_sp: StatePrepCircuit) -> None:
    """Test that the optimal verification is also correct for higher distance.

    Include Z errors.
    Due to time constraints, we set the timeout for each search to 2 seconds.
    """
    circ = color_code_d5_sp

    circ_ver_full_ft = gate_optimal_verification_circuit(circ, max_ancillas=3, max_timeout=5, full_fault_tolerance=True)
    circ_ver_x_ft = gate_optimal_verification_circuit(circ, max_ancillas=3, max_timeout=5, full_fault_tolerance=False)
    assert circ_ver_full_ft.num_nonlocal_gates() > circ_ver_x_ft.num_nonlocal_gates()
    assert circ_ver_full_ft.depth() > circ_ver_x_ft.depth()


def test_full_ft_heuristic_cc5(color_code_d5_sp: StatePrepCircuit) -> None:
    """Test that the optimal verification is also correct for higher distance.

    Include Z errors.
    """
    circ = color_code_d5_sp

    circ_ver_full_ft = heuristic_verification_circuit(circ, full_fault_tolerance=True)
    circ_ver_x_ft = heuristic_verification_circuit(circ, full_fault_tolerance=False)
    assert circ_ver_full_ft.num_nonlocal_gates() > circ_ver_x_ft.num_nonlocal_gates()
    assert circ_ver_full_ft.depth() > circ_ver_x_ft.depth()


def test_error_detection_code() -> None:
    """Test that different circuits are obtained when using an error detection code."""
    code = CSSCode.from_code_name("carbon")
    circ = heuristic_prep_circuit(code)

    circ_ver_correction = gate_optimal_verification_circuit(
        circ, max_ancillas=3, max_timeout=5, full_fault_tolerance=False
    )

    circ.set_error_detection(True)
    circ_ver_detection = gate_optimal_verification_circuit(
        circ, max_ancillas=3, max_timeout=5, full_fault_tolerance=False
    )

    assert circ_ver_detection.num_qubits > circ_ver_correction.num_qubits
    assert circ_ver_detection.num_nonlocal_gates() > circ_ver_correction.num_nonlocal_gates()
