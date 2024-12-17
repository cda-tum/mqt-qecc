"""General circuit constructions."""

from __future__ import annotations

from qiskit import QuantumCircuit, QuantumRegister


def reorder_qubits(circ: QuantumCircuit, qubit_mapping: dict[int, int]) -> QuantumCircuit:
    """Reorders the qubits in a QuantumCircuit based on the given mapping.

    Parameters:
        circuit (QuantumCircuit): The original quantum circuit.
        qubit_mapping (dict[int, int]): A dictionary mapping original qubit indices to new qubit indices.

    Returns:
        QuantumCircuit: A new quantum circuit with qubits reordered.
    """
    # Validate the qubit_mapping
    if sorted(qubit_mapping.keys()) != list(range(len(circ.qubits))) or sorted(qubit_mapping.values()) != list(
        range(len(circ.qubits))
    ):
        msg = "Invalid qubit_mapping: It must be a permutation of the original qubit indices."
        raise ValueError(msg)

    # Create a new quantum register
    num_qubits = len(circ.qubits)
    new_register = QuantumRegister(num_qubits, "q")
    new_circuit = QuantumCircuit(new_register)

    # Remap instructions based on the qubit_mapping
    for instruction, qubits, clbits in circ.data:
        new_qubits = [new_register[qubit_mapping[circ.find_bit(q)[0]]] for q in qubits]
        new_circuit.append(instruction, new_qubits, clbits)

    return new_circuit
