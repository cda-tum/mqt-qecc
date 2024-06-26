from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from mqt.qecc import CSSCode
from mqt.qecc.ft_stateprep import (
    NoisyNDFTStatePrepSimulator,
    heuristic_prep_circuit,
    heuristic_verification_circuit,
)

if TYPE_CHECKING:  # pragma: no cover
    from qiskit import QuantumCircuit


@pytest.fixture()
def steane_code() -> CSSCode:
    return CSSCode.from_code_name("steane")


@pytest.fixture()
def non_ft_steane_zero(steane_code: CSSCode) -> QuantumCircuit:
    return heuristic_prep_circuit(steane_code).circ


@pytest.fixture()
def ft_steane_circ(steane_code: CSSCode) -> QuantumCircuit:
    circ = heuristic_prep_circuit(steane_code)
    return heuristic_verification_circuit(circ)


def test_non_ft_sim(steane_code: CSSCode, non_ft_steane_zero: QuantumCircuit):
    tol = 5e-4
    p = 1e-3
    lower = 1e-4
    simulator = NoisyNDFTStatePrepSimulator(non_ft_steane_zero, steane_code, p=p)
    p_l, _, _, _ = simulator.logical_error_rate()

    assert p_l - tol > lower


def test_ft_sim(steane_code: CSSCode, ft_steane_circ: QuantumCircuit):
    tol = 5e-4
    p = 1e-3
    lower = 1e-4
    simulator = NoisyNDFTStatePrepSimulator(ft_steane_circ, steane_code, p=p)
    p_l, _, _, _ = simulator.logical_error_rate()

    assert p_l - tol < lower
