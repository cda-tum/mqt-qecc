from mqt.qecc.cc_decoder import hexagonal_color_code
from mqt.qecc.cc_decoder.decoder import LightsOut


def test_hex_layout() -> None:
    code = hexagonal_color_code.HexagonalColorCode(distance=3)
    assert len(code.data_qubits) == 7
    code.compute_logical()
    assert len(code.L) == 1
    assert code.H.shape == (3, 7)
    assert code.distance == 3


def test_lo() -> None:
    code = hexagonal_color_code.HexagonalColorCode(distance=3)
    lo = LightsOut(code.faces_to_qubits, code.qubits_to_faces)

    assert len(lo.lights_to_switches) == 3
    assert len(lo.lights_to_switches[0]) == 4
    assert lo.optimizer is not None
