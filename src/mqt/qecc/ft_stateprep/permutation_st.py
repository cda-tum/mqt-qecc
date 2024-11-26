"""This module implements the features of state preparation with qubit permutation."""
from __future__ import annotations

import logging

from qiskit import QuantumCircuit

logger = logging.getLogger(__name__)

def get_stean_prep_circ() -> QuantumCircuit:
    """Function that builds up a state preparation circuit for the Steane code."""
    circ = QuantumCircuit(7)
    circ.h(0)
    circ.h(1)
    circ.h(2)
    circ.cx(0,6)
    circ.cx(1,2)
    circ.cx(3,5)
    circ.cx(0,4)
    circ.cx(1,5)
    circ.cx(0,2)
    circ.cx(3,4)
    circ.cx(5,6)
    return circ
