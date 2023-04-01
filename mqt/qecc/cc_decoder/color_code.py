""" A color code class"""

from __future__ import annotations

from enum import Enum

import numpy as np
from ldpc import mod2


class LatticeType(Enum):
    HEXAGON = "hexagon"
    SQUARE_OCTAGON = "square_octagon"


class ColorCode:
    def __init__(self, distance: int, type: LatticeType):
        self.distance = distance
        self.ancilla_qubits: set[tuple[int, int]] = set()
        self.data_qubits: set[tuple[int, int]] = set()
        self.qubits_to_faces: dict[int, list[int]] = {}
        self.faces_to_qubits: dict[int, list[int]] = {}
        self.lattice = type
        self.add_qubits()
        self.H = np.zeros((len(self.ancilla_qubits), len(self.data_qubits)), dtype=int)
        self.construct_layout()
        self.compute_logical()
        self.n = len(self.qubits_to_faces)
        self.k = self.L.shape[1]

    def add_qubits(self):
        print("base add qubits")
        pass

    def construct_layout(self):
        print("base construct layout")
        pass

    def compute_logical(self) -> None:
        """Compute the logical matrix L"""
        ker_hx = mod2.nullspace(self.H)  # compute the kernel basis of hx
        im_hz_transp = mod2.row_basis(self.H)  # compute the image basis of hz.T
        log_stack = np.vstack([im_hz_transp, ker_hx])
        pivots = mod2.row_echelon(log_stack.T)[3]
        log_op_indices = [i for i in range(im_hz_transp.shape[0], log_stack.shape[0]) if i in pivots]
        self.L = log_stack[log_op_indices]

    def get_syndrome(self, error: np.ndarray) -> np.ndarray:
        return self.H @ error % 2

    def check_if_logical_error(self, residual: np.ndarray) -> bool:
        return (self.L @ residual % 2).any()
