import pathlib
from typing import Any, cast

import pytest
from pytest_console_scripts import ScriptRunner
from qiskit import QuantumCircuit

qasm_circuit = (
        "OPENQASM 2.0;\n"
        + 'include "qelib1.inc";\n'
        + "qreg q[1];\n"
        + "creg c[1];\n"
        + "x q[0];\n"
        + "measure q[0] -> c[0];\n"
)


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
    circ = QuantumCircuit().from_qasm_str(qasm_circuit)
    circ.qasm(filename="dummyCircuit.qasm")
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
        "dummyCircuit.qasm",
    )
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert ret.success


def test_failing_simulators(script_runner: ScriptRunner) -> None:
    """Testing the script with unsupported ecc."""
    circ = QuantumCircuit().from_qasm_str(qasm_circuit)
    circ.qasm(filename="dummyCircuit.qasm")
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
        "dummyCircuit.qasm",
    )
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert not ret.success
    assert "No ECC found for" in ret.stderr


def test_unavailable_backend(script_runner: ScriptRunner) -> None:
    """Testing the script with unsupported ecc."""
    circ = QuantumCircuit().from_qasm_str(qasm_circuit)
    circ.qasm(filename="dummyCircuit.qasm")
    ret = script_runner.run(
        "ecc_framework_qiskit_wrapper",
        "-fs",
        "dummyBackedn",
        "-f",
        "dummyCircuit.qasm",
    )
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert not ret.success
    assert "Unknown backend specified" in ret.stderr


def test_unavailable_error(script_runner: ScriptRunner) -> None:
    """Testing the script with unsupported ecc."""
    circ = QuantumCircuit().from_qasm_str(qasm_circuit)
    circ.qasm(filename="dummyCircuit.qasm")
    ret = script_runner.run(
        "ecc_framework_qiskit_wrapper",
        "-m",
        "K",
        "-f",
        "dummyCircuit.qasm",
    )
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert not ret.success
    assert "Unknown error typ provided" in ret.stderr


def test_statevector_simulators(script_runner: ScriptRunner) -> None:
    """Testing the script with another backend."""
    circ = QuantumCircuit().from_qasm_str(qasm_circuit)
    circ.qasm(filename="dummyCircuit.qasm")
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
        "dummyCircuit.qasm",
    )
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert ret.success


def test_save_circuit(script_runner: ScriptRunner) -> None:
    """Testing the script with another backend."""
    circ = QuantumCircuit().from_qasm_str(qasm_circuit)
    circ.qasm(filename="dummyCircuit.qasm")
    ret = script_runner.run(
        "ecc_framework_qiskit_wrapper",
        "-e",
        "dummyCircuitWithEcc.qasm",
        "-f",
        "dummyCircuit.qasm",
    )
    for circuit_to_delete in ["dummyCircuit.qasm", "dummyCircuitWithEcc.qasm"]:
        file_to_remove = pathlib.Path(circuit_to_delete)
        file_to_remove.unlink()
    assert ret.success
