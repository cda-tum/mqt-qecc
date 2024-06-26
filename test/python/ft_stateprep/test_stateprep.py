"""Test synthesis of state preparation and verification circuits for FT state preparation."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from ldpc import mod2
from qiskit.quantum_info import Clifford

from mqt.qecc import CSSCode
from mqt.qecc.ft_stateprep import (
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


@pytest.fixture()
def steane_code() -> CSSCode:
    """Return the Steane code."""
    return CSSCode.from_code_name("steane")


@pytest.fixture()
def surface_code() -> CSSCode:
    """Return the surface code."""
    return CSSCode.from_code_name("surface")


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


@pytest.mark.parametrize("code_name", ["steane", "tetrahedral", "surface", "cc_4_8_8"])
def test_heuristic_prep_consistent(code_name: str) -> None:
    """Check that heuristic_prep_circuit returns a valid circuit with the correct stabilizers."""
    code = CSSCode.from_code_name(code_name)
    sp_circ = heuristic_prep_circuit(code)
    circ = sp_circ.circ
    max_cnots = np.sum(code.Hx) + np.sum(code.Hz)

    assert circ.num_qubits == code.n
    assert circ.num_nonlocal_gates() <= max_cnots

    x, z = get_stabs(circ)
    assert eq_span(code.Hx, x)
    assert eq_span(np.vstack((code.Hz, code.Lz)), z)


@pytest.mark.parametrize("code_name", ["steane", "surface"])
def test_gate_optimal_prep_consistent(code_name: str) -> None:
    """Check that gate_optimal_prep_circuit returns a valid circuit with the correct stabilizers."""
    code = CSSCode.from_code_name(code_name)
    sp_circ = gate_optimal_prep_circuit(code, max_timeout=2)
    assert sp_circ is not None
    circ = sp_circ.circ
    max_cnots = np.sum(code.Hx) + np.sum(code.Hz)

    assert circ.num_qubits == code.n
    assert circ.num_nonlocal_gates() <= max_cnots

    x, z = get_stabs(circ)
    assert eq_span(code.Hx, x)
    assert eq_span(np.vstack((code.Hz, code.Lz)), z)


@pytest.mark.parametrize("code_name", ["steane", "surface"])
def test_depth_optimal_prep_consistent(code_name: str) -> None:
    """Check that depth_optimal_prep_circuit returns a valid circuit with the correct stabilizers."""
    code = CSSCode.from_code_name(code_name)
    sp_circ = gate_optimal_prep_circuit(code, max_timeout=2)
    assert sp_circ is not None
    circ = sp_circ.circ
    max_cnots = np.sum(code.Hx) + np.sum(code.Hz)

    assert circ.num_qubits == code.n
    assert circ.num_nonlocal_gates() <= max_cnots

    x, z = get_stabs(circ)
    assert eq_span(code.Hx, x)
    assert eq_span(np.vstack((code.Hz, code.Lz)), z)


def test_optimal_steane_verification_circuit(steane_code: CSSCode) -> None:
    """Test that the optimal verification circuit for the Steane code is correct."""
    circ = heuristic_prep_circuit(steane_code)
    ver_stabs_layers = gate_optimal_verification_stabilizers(circ, x_errors=True, max_timeout=2)

    assert len(ver_stabs_layers) == 1  # 1 Ancilla measurement

    ver_stabs = ver_stabs_layers[0]

    assert np.sum(ver_stabs) == 3  # 3 CNOTs
    z_gens = np.vstack((steane_code.Hz, steane_code.Lz))

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


def test_heuristic_steane_verification_circuit(steane_code: CSSCode) -> None:
    """Test that the optimal verification circuit for the Steane code is correct."""
    circ = heuristic_prep_circuit(steane_code)
    ver_stabs_layers = heuristic_verification_stabilizers(circ, x_errors=True)

    assert len(ver_stabs_layers) == 1  # 1 Ancilla measurement

    ver_stabs = ver_stabs_layers[0]

    assert np.sum(ver_stabs) == 3  # 3 CNOTs
    z_gens = np.vstack((steane_code.Hz, steane_code.Lz))

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
