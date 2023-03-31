from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import pytest_mock
from pytest_console_scripts import ScriptRunner

from mqt.qecc import cc_decoder
from mqt.qecc.cc_decoder.decoder import LightsOut, simulate_error_rate


@pytest.fixture
def distance() -> int:
    """Distance of the hexagonal color code."""
    return 3


@pytest.fixture
def p() -> float:
    """Error rate."""
    return 0.1


@pytest.fixture
def nr_sims() -> int:
    """Number of simulations."""
    return 1


@pytest.fixture
def results_dir() -> str:
    """Directory to store results."""
    return "./results"


@pytest.fixture
def code(distance: int) -> cc_decoder.HexagonalColorCode:
    """Hexagonal color code."""
    return cc_decoder.HexagonalColorCode(distance=distance)


def test_hex_layout(code: cc_decoder.HexagonalColorCode, distance: int) -> None:
    """Test the construction of the hexagonal color code."""
    assert len(code.data_qubits) == 7
    code.compute_logical()
    assert len(code.L) == 1
    assert code.H.shape == (distance, 7)
    assert code.distance == distance


def test_lo(code: cc_decoder.HexagonalColorCode, distance: int) -> None:
    """Test the construction of the lights out problem."""
    lo = LightsOut(code.faces_to_qubits, code.qubits_to_faces)

    assert len(lo.lights_to_switches) == distance
    assert len(lo.lights_to_switches[0]) == 4
    assert lo.optimizer is not None


def test_maxsat_decoder(code: cc_decoder.HexagonalColorCode) -> None:
    """Test the MaxSAT decoder on a small example."""
    lo = LightsOut(code.faces_to_qubits, code.qubits_to_faces)
    lo.preconstruct_z3_instance()
    lights = [True, False, False]
    est, _, _ = lo.solve(lights=lights)
    assert len(est) == 7
    assert est[1] == 1


def test_simulate(distance: int, p: float, nr_sims: int) -> None:
    """Test the general simulation function."""
    res = simulate_error_rate(distance, p, nr_sims)
    assert res is not None
    assert res["distance"] == distance
    assert res["p"] == p


def check_and_load_json(file_name: str, results_dir: str) -> dict[str, Any]:
    """Check that the results directory contains exactly one file with the given name and load it as JSON."""
    results_path = Path(results_dir)
    assert results_path.exists()
    assert results_path.is_dir()
    assert len(list(results_path.iterdir())) == 1
    result_file = results_path / file_name
    assert result_file.exists()
    assert result_file.is_file()
    with result_file.open("r") as f:
        result = json.load(f)

    for file in results_path.iterdir():
        file.unlink()
    results_path.rmdir()

    return result


def test_z3_solver(script_runner: ScriptRunner, distance: int, p: float, nr_sims: int, results_dir: str) -> None:
    """Test the Z3 solver."""
    ret = script_runner.run(
        "mqt.qecc.cc-decoder",
        str(distance),
        str(p),
        "--nr_sims",
        str(nr_sims),
        "--results_dir",
        results_dir,
    )
    assert ret.success
    assert ret.stderr == ""

    result = check_and_load_json(f"distance={distance},p={round(p, 4)},solver=z3.json", results_dir)
    assert result is not None
    assert result["distance"] == distance
    assert result["p"] == p
    assert result["logical_error_rate"] is not None
    assert result["logical_error_rate_eb"] is not None
    assert result["min_wt_logical_err"] is not None
    assert result["preconstr_time"] > 0.0
    assert result["avg_constr_time"] > 0.0
    assert result["avg_solve_time"] > 0.0
    assert result["avg_total_time"] > 0.0


def test_external_solver(
    script_runner: ScriptRunner,
    mocker: pytest_mock.MockerFixture,
    distance: int,
    p: float,
    nr_sims: int,
    results_dir: str,
) -> None:
    """Mock the call to an external solver."""
    solver = "mock_solver"
    decoder = "maxsat"

    # mock the subprocess.run call in the `solve` method of the `LightsOut` class
    mocker.patch("subprocess.run", return_value=None)
    ret = script_runner.run(
        "mqt.qecc.cc-decoder",
        str(distance),
        str(p),
        "--nr_sims",
        str(nr_sims),
        "--results_dir",
        results_dir,
        "--solver",
        solver,
        "--decoder",
        decoder,
    )
    assert ret.success
    assert ret.stderr == ""

    result = check_and_load_json(f"distance={distance},p={round(p, 4)},solver={solver}.json", results_dir)

    solver_output_file = Path(f"solver-out_{solver}.txt")
    if solver_output_file.exists():
        solver_output_file.unlink()

    assert result is not None
    assert result["distance"] == distance
    assert result["p"] == p
    assert result["logical_error_rate"] == 0.0
    assert result["logical_error_rate_eb"] == 0.0
    assert result["min_wt_logical_err"] == 0
    assert result["preconstr_time"] > 0.0
    assert result["avg_constr_time"] > 0.0
    assert result["avg_solve_time"] > 0.0
    assert result["avg_total_time"] > 0.0


def test_tn_decoder(script_runner: ScriptRunner, distance: int, p: float, nr_sims: int, results_dir: str) -> None:
    """Test the TN decoder."""
    decoder = "tn"

    ret = script_runner.run(
        "mqt.qecc.cc-decoder",
        str(distance),
        str(p),
        "--nr_sims",
        str(nr_sims),
        "--results_dir",
        results_dir,
        "--decoder",
        decoder,
    )
    assert ret.success
    assert ret.stderr == ""

    result = check_and_load_json(f"distance={distance},p={round(p, 4)}.json", results_dir)
    assert result is not None
    assert result["n_k_d"][-1] == distance
    assert result["error_probability"] == p
    assert result["n_run"] == nr_sims
    assert result["wall_time"] > 0.0
