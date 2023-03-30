import numpy as np
from square_octagon_color_code import SquareOctagonColorCode


def test_d3():
    d3_layout = SquareOctagonColorCode(3)

    assert d3_layout.data_qubits == {(0, 1), (2, 1), (6, 1), (3, 2), (5, 2), (3, 4), (5, 4)}
    assert d3_layout.ancilla_qubits == {(4, 0), (1, 3), (4, 3)}


def test_number_of_qubits():
    for d in range(3, 16, 2):
        layout = SquareOctagonColorCode(d)
        assert len(layout.data_qubits) == 1 / 2 * d**2 + d - 1 / 2
        assert len(layout.ancilla_qubits) == (1 / 2 * d**2 + d - 1 / 2) // 2


def test_H():
    d3_layout = SquareOctagonColorCode(3)
    assert (d3_layout.H == np.array([[0, 0, 1, 0, 1, 1, 1], [0, 1, 0, 1, 0, 1, 1], [1, 1, 1, 0, 0, 1, 0]])).all()
