# Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

"""Circuit representations."""

from __future__ import annotations

from typing import TYPE_CHECKING

import stim
from qiskit import QuantumCircuit

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Iterable


class CNOTCircuit:
    """Represents a restricted quantum circuit composed of CNOT gates with optional qubit initialization."""

    def __init__(self) -> None:
        """Initialize an empty CNOT circuit."""
        self.cnots: list[tuple[int, int]] = []
        self.initializations: dict[int, str] = {}  # Dictionary mapping qubit index to initialization type ('Z' or 'X')

    def add_cnot(self, control: int, target: int) -> None:
        """Add a single CNOT gate to the circuit.

        Args:
            control: The control qubit index.
            target: The target qubit index.
        """
        if control < 0 or target < 0:
            msg = "Control and target qubits must have non-negative indices."
            raise ValueError(msg)
        if control == target:
            msg = "Control and target qubits cannot be the same."
            raise ValueError(msg)
        self.cnots.append((control, target))

    def add_cnots(self, cnot_pairs: Iterable[tuple[int, int]]) -> None:
        """Add multiple CNOT gates to the circuit.

        Args:
            cnot_pairs: An iterable of (control, target) pairs.
        """
        for control, target in cnot_pairs:
            self.add_cnot(control, target)

    def initialize_qubit(self, qubit: int, basis: str) -> None:
        """Initialize a qubit in the specified basis.

        Args:
            qubit: The qubit index to initialize.
            basis: The basis for initialization ('Z' or 'X').
        """
        if qubit < 0:
            msg = "Qubit index must be non-negative."
            raise ValueError(msg)
        if basis.capitalize() not in {"Z", "X"}:
            msg = "Initialization basis must be 'Z' or 'X'."
            raise ValueError(msg)
        self.initializations[qubit] = basis

    def to_stim_circuit(self) -> stim.Circuit:
        """Convert the CNOT circuit to a stim.Circuit.

        Returns:
            A stim.Circuit representation of the CNOT circuit.
        """
        stim_circuit = stim.Circuit()

        # Add initializations
        for qubit, basis in self.initializations.items():
            stim_circuit.append("R" + basis, [qubit])

        # Add CNOT gates
        stim_circuit.append_operation("CX", [qubit for pair in self.cnots for qubit in pair])

        return stim_circuit

    def to_qiskit_circuit(self) -> QuantumCircuit:
        """Convert the CNOT circuit to a qiskit.QuantumCircuit.

        Returns:
            A qiskit.QuantumCircuit representation of the CNOT circuit.
        """
        return QuantumCircuit.from_qasm_str(self.to_stim_circuit().to_qasm(open_qasm_version=2))

    def is_state(self) -> bool:
        """Check if all qubits used in the circuit are initialized.

        Returns:
            True if all qubits involved in CNOT operations are initialized, False otherwise.
        """
        used_qubits = {qubit for control, target in self.cnots for qubit in (control, target)}
        return used_qubits.issubset(self.initializations.keys())
