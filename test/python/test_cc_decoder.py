from unittest.mock import MagicMock, patch

from mqt.qecc.cc_decoder import hexagonal_color_code
from mqt.qecc.cc_decoder.decoder import LightsOut, simulate_error_rate


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


def test_maxsat_construction() -> None:
    code = hexagonal_color_code.HexagonalColorCode(distance=3)
    lo = LightsOut(code.faces_to_qubits, code.qubits_to_faces)
    lo.preconstruct_z3_instance()
    lights = [True, False, False]
    est, _, _ = lo.solve(lights=lights, solver_path="z3")
    assert len(est) == 7
    assert est[1] == 1


def test_simulate() -> None:
    distance = 3
    p = 0.1
    nr_sims = 1
    solver_path = "z3"

    res = simulate_error_rate(distance, p, nr_sims, solver_path)
    assert res is not None
    assert res["distance"] == 3
    assert res["p"] == 0.1


@patch("mqt.qecc.cc_decoder.decoder.subprocess.run")
def test_solve_non_z3(mock_run):
    mock_stdout = MagicMock()
    mock_stdout.configure_mock(**{"stdout.decode.return_value": "solved"})

    mock_run.return_value = mock_stdout
    code = hexagonal_color_code.HexagonalColorCode(distance=3)
    lo = LightsOut(code.faces_to_qubits, code.qubits_to_faces)
    lo.preconstruct_z3_instance()
    lights = [True, False, False]
    lo.solve(lights=lights, solver_path="exact")
    assert True
