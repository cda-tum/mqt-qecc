from __future__ import annotations

import argparse

from mqt import qecc
from qiskit import Aer, QuantumCircuit, execute, providers
from qiskit.result import counts
from qiskit_aer.noise import (
    NoiseModel,
    QuantumError,
    amplitude_damping_error,
    depolarizing_error,
)
from qiskit_aer.noise.errors import pauli_error


def compose_error(error: QuantumError, new_error: QuantumError) -> QuantumError:
    if error is None:
        error = new_error
    else:
        error = error.compose(new_error)
    return error


def create_noise_model(n_model: str, p_error: float) -> NoiseModel:
    # Create an empty noise model
    noise_model = NoiseModel()
    error = None
    for char in n_model:
        if char == "B":
            # Add a bit flip error channel
            new_error = pauli_error([("X", p_error), ("I", 1 - p_error)])
            error = compose_error(error, new_error)
        elif char == "D":
            # Add a depolarization error channel
            new_error = depolarizing_error(p_error, 1)
            error = compose_error(error, new_error)
        elif char == "A":
            # simple amplitude damping error channel which mimics energy loss to the environment (T1-error)
            new_error = amplitude_damping_error(p_error)
            error = compose_error(error, new_error)
        elif char == "P":
            # Add a phase flip error channel
            new_error = pauli_error([("Z", p_error), ("I", 1 - p_error)])
            error = compose_error(error, new_error)

        else:
            raise ValueError("Unknown error type in noise model: " + char)

    assert error is not None

    # single qubit and multi qubit operations are treated the same at the moment
    noise_model.add_all_qubit_quantum_error(
        error, ["u1", "u2", "u3", "h", "id", "t", "tdg", "sdg", "rx", "ry", "rz", "s"]
    )
    noise_model.add_all_qubit_quantum_error(error.tensor(error), ["cx", "swap"])
    noise_model.add_all_qubit_quantum_error(error.tensor(error).tensor(error), ["cswap"])

    return noise_model


def print_simulation_results(result_counts: counts, n_shots: int, threshold_probability: float = 0) -> None:
    printed_results = 0
    summarized_counts: dict[str, int] = {}
    for result_id in result_counts:
        sub_result = result_id.split(" ")[-1]
        if sub_result not in summarized_counts.keys():
            summarized_counts[sub_result] = 0
        summarized_counts[sub_result] += result_counts[result_id]

    for result_id in sorted(summarized_counts.keys()):
        # Print all results > threshold_probability
        if summarized_counts[result_id] / n_shots > threshold_probability or printed_results == 0:
            result_string = str(result_id)
            print("State |" + result_string + "> probability " + str(summarized_counts[result_id] / n_shots))
            printed_results += 1
            if printed_results == 1000:
                break


def main() -> None:
    parser = argparse.ArgumentParser(description="Error correcting codes Qiskit wrapper")
    parser.add_argument(
        "-m",
        type=str,
        default="D",
        help="Define the error_channels (e.g., -m APD), available errors channels are amplitude "
        'damping (A), phase flip (P), bit flip (B), and depolarization (D) (Default="D")',
    )
    parser.add_argument("-p", type=float, default=0.001, help="Set the noise probability (Default=0.001)")
    parser.add_argument(
        "-n", type=int, default=2000, help="Set the number of shots. 0 for deterministic simulation (" "Default=2000)"
    )
    parser.add_argument("-s", type=int, default=0, help="Set a seed (Default=0)")
    parser.add_argument("-f", type=str, required=True, help="Path to a OpenQASM file")
    parser.add_argument(
        "-e",
        type=str,
        required=False,
        default=None,
        help="Export circuit, with error correcting code applied, as OpenQASM circuit instead of "
        'simulation it (e.g., -e "/path/to/new/openqasm_file") (Default=None)',
    )
    parser.add_argument(
        "-fs",
        type=str,
        default="none",
        help='Specify a simulator (Default: "statevector_simulator" for simulation without noise, '
        '"aer_simulator_density_matrix", for deterministic noise-aware simulation'
        '"aer_simulator_statevector", for stochastic noise-aware simulation). Available: ' + str(Aer.backends()),
    )
    parser.add_argument(
        "-ecc",
        type=str,
        default="Q7Steane",
        help="Specify a ecc to be applied to the circuit. Currently available are none, Q3Shor, Q5Laflamme, "
        "Q7Steane, Q9Shor, Q9Surface, and Q18Surface (Default=Q7Steane)",
    )
    parser.add_argument(
        "-fq",
        type=int,
        default=100,
        help="Specify after how many qubit usages error correction is " "applied to it (Default=100)",
    )

    args = parser.parse_args()

    error_channels = args.m
    error_probability = args.p
    number_of_shots = args.n
    seed = args.s
    open_qasm_file = args.f

    if args.fs.lower() == "none":
        forced_simulator = None
    else:
        forced_simulator = args.fs

    if args.ecc.lower() == "none":
        ecc = None
    else:
        ecc = args.ecc

    ecc_frequency = args.fq
    ecc_export_filename = args.e

    if number_of_shots > 0:
        n_shots = number_of_shots
    else:
        n_shots = 2000

    # Creating the noise model
    if error_probability > 0:
        noise_model = create_noise_model(n_model=error_channels, p_error=error_probability)
    else:
        noise_model = NoiseModel()

    circ = QuantumCircuit.from_qasm_file(open_qasm_file)

    if not any(gate[0].name == "measure" for gate in circ.data):
        print("No measurement gates found in circuit. Adding measurement gates to all qubits.")
        circ.measure_all()

    # Initializing the quantum circuit
    if ecc is not None:
        # Applying error correction to the circuit
        result = qecc.apply_ecc(circ, ecc, ecc_frequency)
        circ = QuantumCircuit().from_qasm_str(result["circ"])

    if ecc_export_filename is not None:
        print("Exporting circuit to: " + str(ecc_export_filename))
        circ.qasm(filename=ecc_export_filename)
        return

    size = circ.num_qubits
    simulator_backend = None
    print(
        "_____Trying to simulate with "
        + str(error_channels)
        + " (prob="
        + str(error_probability)
        + ", shots="
        + str(n_shots)
        + ", n_qubits="
        + str(size)
        + ", error correction="
        + str(ecc)
        + ") Error______",
        flush=True,
    )

    if forced_simulator is not None:
        # Setting the simulator backend to the requested one
        try:
            simulator_backend = Aer.get_backend(forced_simulator)
        except providers.exceptions.QiskitBackendNotFoundError:
            raise ValueError(
                "Simulator " + str(forced_simulator) + " not found! Available simulators are: " + str(Aer.backends())
            ) from None
    else:
        print(
            "No backend specified. Setting backend to `aer_simulator_statevector`, which is fast but does not support non-Clifford gates."
        )
        simulator_backend = Aer.get_backend("aer_simulator_statevector")

    result = execute(
        circ,
        backend=simulator_backend,
        shots=n_shots,
        seed_simulator=seed,
        noise_model=noise_model,
    )

    if result.result().status != "COMPLETED":
        raise RuntimeError("Simulation exited with status: " + str(result.result().status))

    result_counts = result.result().get_counts()
    print_simulation_results(result_counts, n_shots)
