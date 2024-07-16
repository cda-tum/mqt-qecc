"""Test the ecc framework."""

from __future__ import annotations

import pathlib
from typing import TYPE_CHECKING

import pytest
from qiskit.qasm2 import dump

if TYPE_CHECKING:
    from pytest_console_scripts import ScriptRunner

import locale

from qiskit import QuantumCircuit


@pytest.fixture
def circ() -> QuantumCircuit:
    """Fixture for a quantum circuit."""
    qasm_circuit = 'OPENQASM 2.0;\n include "qelib1.inc";\n qreg q[1];\n creg c[1];\n x q[0];\n measure q[0] -> c[0];\n'
    return QuantumCircuit().from_qasm_str(qasm_circuit)


@pytest.fixture
def circ_no_measure() -> QuantumCircuit:
    """Fixture for a quantum circuit without measure."""
    qasm_circuit_no_measure = 'OPENQASM 2.0;\n include "qelib1.inc";\n qreg q[1];\n creg c[1];\n x q[0];\n'
    return QuantumCircuit().from_qasm_str(qasm_circuit_no_measure)


@pytest.mark.parametrize("simulator", ["none", "Id", "Q7Steane"])
@pytest.mark.parametrize("noise_models", ["B", "P", "D"])
def test_with_stab_simulator(
    circ: QuantumCircuit, simulator: str, noise_models: str, script_runner: ScriptRunner
) -> None:
    """Testing the script with different parameters."""
    with pathlib.Path("dummyCircuit.qasm").open("w", encoding=locale.getpreferredencoding(False)) as f:
        dump(circ, f)
    ret = script_runner.run([
        "ecc_qiskit_wrapper",
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
    ])
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert ret.success


def test_failing_simulators(circ: QuantumCircuit, script_runner: ScriptRunner) -> None:
    """Testing the script with unsupported ecc."""
    with pathlib.Path("dummyCircuit.qasm").open("w", encoding=locale.getpreferredencoding(False)) as f:
        dump(circ, f)
    ret = script_runner.run([
        "ecc_qiskit_wrapper",
        "-m",
        "D",
        "-p",
        "0.001",
        "-n",
        "1000",
        "-s",
        "1",
        "-fs",
        "extended_stabilizer",
        "-ecc",
        "UnsupportedEcc",
        "-f",
        "dummyCircuit.qasm",
    ])
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert not ret.success
    assert "No ECC found for" in ret.stderr


def test_unavailable_backend(circ: QuantumCircuit, script_runner: ScriptRunner) -> None:
    """Testing the script with unsupported backend."""
    with pathlib.Path("dummyCircuit.qasm").open("w", encoding=locale.getpreferredencoding(False)) as f:
        dump(circ, f)
    ret = script_runner.run(["ecc_qiskit_wrapper", "-fs", "dummyBackend", "-f", "dummyCircuit.qasm"])
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert not ret.success
    assert "Available methods are" in ret.stderr


def test_unavailable_error_type(circ: QuantumCircuit, script_runner: ScriptRunner) -> None:
    """Testing the script with unsupported ecc."""
    with pathlib.Path("dummyCircuit.qasm").open("w", encoding=locale.getpreferredencoding(False)) as f:
        dump(circ, f)
    ret = script_runner.run(["ecc_qiskit_wrapper", "-m", "K", "-f", "dummyCircuit.qasm"])
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert not ret.success
    assert "Unknown error type in noise model: " in ret.stderr


def test_statevector_simulators(circ: QuantumCircuit, script_runner: ScriptRunner) -> None:
    """Testing the simulator with a different simulator."""
    with pathlib.Path("dummyCircuit.qasm").open("w", encoding=locale.getpreferredencoding(False)) as f:
        dump(circ, f)
    ret = script_runner.run([
        "ecc_qiskit_wrapper",
        "-m",
        "BAPD",
        "-p",
        "0.001",
        "-n",
        "50",
        "-s",
        "1",
        "-fs",
        "statevector",
        "-ecc",
        "Id",
        "-f",
        "dummyCircuit.qasm",
    ])
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert ret.success


def test_save_circuit(circ: QuantumCircuit, script_runner: ScriptRunner) -> None:
    """Saving a circuit after applying an ECC."""
    with pathlib.Path("dummyCircuit.qasm").open("w", encoding=locale.getpreferredencoding(False)) as f:
        dump(circ, f)
    ret = script_runner.run([
        "ecc_qiskit_wrapper",
        "-e",
        "dummyCircuitWithEcc.qasm",
        "-f",
        "dummyCircuit.qasm",
    ])
    for circuit_to_delete in ["dummyCircuit.qasm", "dummyCircuitWithEcc.qasm"]:
        file_to_remove = pathlib.Path(circuit_to_delete)
        file_to_remove.unlink()
    assert ret.success


def test_circuit_without_measurements(circ_no_measure: QuantumCircuit, script_runner: ScriptRunner) -> None:
    """Testing circuit without ecc."""
    with pathlib.Path("dummyCircuit.qasm").open("w", encoding=locale.getpreferredencoding(False)) as f:
        dump(circ_no_measure, f)
    ret = script_runner.run(["ecc_qiskit_wrapper", "-f", "dummyCircuit.qasm"])
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert ret.success


def test_trying_to_use_stabilizer_simulator(circ: QuantumCircuit, script_runner: ScriptRunner) -> None:
    """Testing circuit without ecc."""
    with pathlib.Path("dummyCircuit.qasm").open("w", encoding=locale.getpreferredencoding(False)) as f:
        dump(circ, f)
    ret = script_runner.run(["ecc_qiskit_wrapper", "-f", "dummyCircuit.qasm", "-m", "A"])
    file_to_remove = pathlib.Path("dummyCircuit.qasm")
    file_to_remove.unlink()
    assert not ret.success
    assert "Simulation failed" in ret.stderr
