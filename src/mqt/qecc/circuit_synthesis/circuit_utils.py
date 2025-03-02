"""General circuit constructions."""

from __future__ import annotations

import stim
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


def relabel_qubits(circ: stim.Circuit, qubit_mapping: dict[int, int] | int) -> stim.Circuit:
    """Relabels the qubits in a stim circuit based on the given mapping.

    Parameters:
        circ (stim.Circuit): The original stim circuit.
        qubit_mapping (dict[int, int] | int): Either a dictionary mapping original qubit indices to new qubit indices or a constant offset to add to all qubit indices.

    Returns:
        stim.Circuit: A new stim circuit with qubits relabeled.
    """
    new_circ = stim.Circuit()
    for op in circ:
        if isinstance(qubit_mapping, dict):
            relabelled_qubits = [qubit_mapping[q.value] for q in op.targets_copy()]
        else:
            relabelled_qubits = [q.value + qubit_mapping for q in op.targets_copy()]
        new_circ.append(op.name, relabelled_qubits)
    return new_circ
