"""General circuit constructions."""

from __future__ import annotations

from qiskit import QuantumCircuit, QuantumRegister
import stim

from typing import TYPE_CHECKING

from ..codes.pauli import StabilizerTableau

    
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


def apply_clifford_circuit(stabs: StabilizerTableau, circ: QuantumCircuit | stim.Circuit) -> SymplecticMatrix:
    """Apply a Clifford circuit to a stabilizer tableau.

    Args:
        stabs (StabilizerTableau): The stabilizer tableau.
        circ (QuantumCircuit | stim.Circuit): The Clifford circuit.

    Returns:
        The tableau after applying the Clifford circuit.
    """
    if isinstance(circ, stim.Circuit):
        circ = QuantumCircuit.from_qasm_str(circ.to_qasm(open_qasm_version="2.0"))

    n = QuantumCircuit.num_qubits
    assert n == stabs.n, "The number of qubits in the circuit must match the number of qubits in the tableau."
        
    # Initialize the new tableau
    new_stabs = stabs.copy()
    
    for gate in circ:
        name = gate.name
        qubit = gate.qubits[0]

        if name == "H":
            new_stabs.apply_h(qubit)
        elif name == "S":
            new_stabs.apply_s(qubit)
        elif name == "X":
            new_stabs.apply_x(qubit)
        elif name == "Y":
            new_stabs.apply_y(qubit)
        elif name == "Z":
            new_stabs.apply_z(qubit)
        elif name == "CX":
            ctrl, tgt = qubit, gate.qubits[1]
            new_stabs.apply_cx(ctrl, tgt)
        elif name == "CZ":
            ctrl, tgt = qubit, gate.qubits[1]
            new_stabs.apply_cz(ctrl, tgt)

    return new_stabs
    
