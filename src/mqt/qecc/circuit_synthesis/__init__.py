"""Methods and utilities for synthesizing fault-tolerant circuits and gadgets."""

from __future__ import annotations

from .encoding import heuristic_encoding_circuit
from .simulation import LutDecoder, NoisyNDFTStatePrepSimulator
from .state_prep import (
    StatePrepCircuit,
    depth_optimal_prep_circuit,
    gate_optimal_prep_circuit,
    gate_optimal_verification_circuit,
    gate_optimal_verification_stabilizers,
    heuristic_prep_circuit,
    heuristic_verification_circuit,
    heuristic_verification_stabilizers,
    naive_verification_circuit,
)

__all__ = [
    "LutDecoder",
    "NoisyNDFTStatePrepSimulator",
    "StatePrepCircuit",
    "depth_optimal_prep_circuit",
    "gate_optimal_prep_circuit",
    "gate_optimal_verification_circuit",
    "gate_optimal_verification_stabilizers",
    "heuristic_encoding_circuit",
    "heuristic_prep_circuit",
    "heuristic_verification_circuit",
    "heuristic_verification_stabilizers",
    "naive_verification_circuit",
]
