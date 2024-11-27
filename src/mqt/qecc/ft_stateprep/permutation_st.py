"""This module implements the features of state preparation with qubit permutation."""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from qiskit import QuantumCircuit

from mqt.qecc.codes.css_code import CSSCode

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt

logger = logging.getLogger(__name__)

def get_stean_prep_circ() -> QuantumCircuit:
    """Function that builds up a state preparation circuit for the Steane code."""
    circ = QuantumCircuit(7)
    circ.h(0)
    circ.h(1)
    circ.h(3)
    circ.cx(0,6)
    circ.cx(1,2)
    circ.cx(3,5)
    circ.cx(0,4)
    circ.cx(1,5)
    circ.cx(0,2)
    circ.cx(3,4)
    circ.cx(5,6)
    return circ

def get_steane_css_object() -> CSSCode:
    """Return the Steane Code as CSS object."""
    distance: int = 3
    hx : npt.NDArray[np.int8] = np.array([[1,1,1,1,0,0,0],
                                          [0,1,1,0,1,1,0],
                                         [0,0,1,1,0,6,7]])
    hz : npt.NDArray[np.int8] = np.array([[1,1,1,1,0,0,0],
                                          [0,1,1,0,1,1,0],
                                         [0,0,1,1,0,6,7]])
    return CSSCode(distance, hx, hz, distance, distance)
