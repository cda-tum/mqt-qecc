from typing import Any, cast

import pytest
from pytest_console_scripts import ScriptRunner


@pytest.fixture(
    params=[
        "none",
        "Id",
        "Q3Shor",
        "Q7Steane",
    ]
)
def simulator(request: Any) -> str:
    return cast(str, request.param)


@pytest.fixture(params=["B", "P", "D"])
def noise_models(request: Any) -> str:
    return cast(str, request.param)


def test_with_stab_simulator(simulator: str, noise_models: str, script_runner: ScriptRunner) -> None:
    """Testing the script with different parameters"""
    print("Type is: " + str(type(script_runner)))
    ret = script_runner.run(
        "ecc_framework_qiskit_wrapper",
        "-m",
        noise_models,
        "-p",
        "0.001",
        "-n",
        "1000",
        "-ecc",
        simulator,
        "-f",
        "test/python/ExampleCircuit.qasm",
    )
    assert ret.success


def test_failing_simulators(script_runner: ScriptRunner) -> None:
    """Testing the script with unsupported ecc."""
    ret = script_runner.run(
        "ecc_framework_qiskit_wrapper",
        "-m",
        "D",
        "-p",
        "0.001",
        "-n",
        "1000",
        "-s",
        "1",
        "-fs",
        "aer_simulator_extended_stabilizer",
        "-ecc",
        "UnsupportedEcc",
        "-f",
        "test/python/ExampleCircuit.qasm",
    )

    assert not ret.success
    assert "No ECC found for" in ret.stderr

    def test_unavailable_backend() -> None:
        """Testing the script with unsupported ecc."""

    ret = script_runner.run(
        "ecc_framework_qiskit_wrapper",
        "-fs",
        "dummyBackedn",
        "-f",
        "test/python/ExampleCircuit.qasm",
    )

    assert not ret.success
    assert "Unknown backend specified" in ret.stderr


def test_statevector_simulators(script_runner: ScriptRunner) -> None:
    """Testing the script with unsupported ecc."""
    ret = script_runner.run(
        "ecc_framework_qiskit_wrapper",
        "-m",
        "BAPD",
        "-p",
        "0.001",
        "-n",
        "0",
        "-s",
        "1",
        "-fs",
        "aer_simulator_statevector",
        "-ecc",
        "Id",
        "-f",
        "test/python/ExampleCircuit.qasm",
    )

    assert ret.success
