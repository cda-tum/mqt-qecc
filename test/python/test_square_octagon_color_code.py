"""Test the square octagon color code. Created by Peter-Jan Derks."""

from __future__ import annotations

import numpy as np
import pytest
from mqt.qecc.cc_decoder.square_octagon_color_code import SquareOctagonColorCode


@pytest.fixture()
def code(request) -> SquareOctagonColorCode:
    """Color code on 4.8.8 lattice."""
    return SquareOctagonColorCode(distance=request.param)


@pytest.mark.parametrize("code", list(range(3, 23, 2)), indirect=True)
def test_number_of_qubits(code: SquareOctagonColorCode) -> None:
    """Test the number of qubits for larger distances."""
    assert len(code.data_qubits) == 1 / 2 * code.distance**2 + code.distance - 1 / 2
    assert len(code.ancilla_qubits) == (1 / 2 * code.distance**2 + code.distance - 1 / 2) // 2


@pytest.mark.parametrize("code", [3], indirect=True)
def test_d3(code: SquareOctagonColorCode) -> None:
    """Test coordinates of qubits for distance 3."""
    assert code.data_qubits == {(0, 1), (2, 1), (6, 1), (3, 2), (5, 2), (3, 4), (5, 4)}
    assert code.ancilla_qubits == {(4, 0), (1, 3), (4, 3)}


@pytest.mark.parametrize("code", [3], indirect=True)
def test_H(code: SquareOctagonColorCode) -> None:
    """Test the parity check matrix for distance 3."""
    assert np.array_equal(code.H, np.array([[0, 0, 1, 0, 1, 1, 1], [0, 1, 0, 1, 0, 1, 1], [1, 1, 1, 0, 0, 1, 0]]))
