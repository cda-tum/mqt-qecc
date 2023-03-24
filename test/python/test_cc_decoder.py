import pytest_mock
from mqt.qecc import cc_decoder
from mqt.qecc.cc_decoder import cli, hexagonal_color_code
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


def test_solve_non_z3(mocker: pytest_mock.MockerFixture) -> None:
    distance: int = 3
    error_rate: float = 0.1
    nr_sims: int = 1
    results_dir: str = "./results"
    solver: str = "z3"
    decoder: str = "maxsat"

    class args:  # noqa: B903
        def __init__(self) -> None:
            self.distance = distance
            self.error_rate = error_rate
            self.nr_sims = nr_sims
            self.results_dir = results_dir
            self.solver = solver
            self.decoder = decoder

    mocker.patch("argparse.ArgumentParser.parse_args", return_value=args())
    mocker.patch("mqt.qecc.cc_decoder.decoder.run", return_value=None)

    cli()
    cc_decoder.decoder.run.assert_called_once_with(distance, error_rate, nr_sims, results_dir, solver)  # type: ignore[attr-defined]
