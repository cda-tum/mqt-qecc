from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
import numpy as np

from mqt.qecc import CSSCode
from mqt.qecc.ft_stateprep import (
    NoisyNDFTStatePrepSimulator,
    heuristic_prep_circuit,
    heuristic_verification_circuit,
    LutDecoder,
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


def test_lut(steane_code: CSSCode):
    """Test the LutDecoder class."""
    lut = LutDecoder(steane_code, init_luts=False)

    assert len(lut.x_lut) == 0
    assert len(lut.z_lut) == 0

    lut.generate_x_lut()
    lut.generate_z_lut()
    
    assert len(lut.x_lut) != 0
    assert lut.x_lut is lut.z_lut  # Code is self dual so luts should be the same

    error_1 = np.zeros(steane_code.n, dtype=np.int8)
    error_1[0] = 1
    
    error_w1 = (steane_code.Hx[0] + error_1) % 2
    syndrome_1 = steane_code.get_x_syndrome(error_w1)
    estimate_1 = lut.decode_x(syndrome_1.astype(np.int8))
    assert steane_code.stabilizer_eq_x_error(estimate_1, error_1)
    assert steane_code.stabilizer_eq_z_error(estimate_1, error_1)
    
    error_2 = np.zeros(steane_code.n, dtype=np.int8)
    error_2[0] = 1
    error_2[1] = 1
    error_w2 = (steane_code.Hx[0] + error_2) % 2
    syndrome_2 = steane_code.get_x_syndrome(error_w2)
    estimate_2 = lut.decode_x(syndrome_2.astype(np.int8))

    # Weight 2 error should have be estimated to be weight 1
    assert not steane_code.stabilizer_eq_x_error(estimate_2, error_2)
    assert np.sum(estimate_2) == 1
    
    error_3 = np.ones((steane_code.n), dtype=np.int8)
    error_w3 = (steane_code.Hx[0] + error_3) % 2
    syndrome_3 = steane_code.get_x_syndrome(error_w3)
    estimate_3 = lut.decode_x(syndrome_3.astype(np.int8))
    # Weight 3 error should have be estimated to be weight 0
    assert not steane_code.stabilizer_eq_x_error(estimate_3, error_3)
    assert steane_code.stabilizer_eq_x_error(estimate_3, np.zeros(steane_code.n))
    assert np.sum(estimate_3) == 0
    

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

