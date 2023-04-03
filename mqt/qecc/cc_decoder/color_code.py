"""A color code class."""

from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING

import numpy as np
from ldpc import mod2

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt


class LatticeType(Enum):
    """The type of lattice used for the color code."""

    HEXAGON = "hexagon"
    SQUARE_OCTAGON = "square_octagon"


class ColorCode:
    """A base class for color codes on a three-valent, three-colourable lattice."""

    def __init__(self, distance: int, lattice: LatticeType) -> None:
        """Initialize the color code."""
        self.distance = distance
        self.ancilla_qubits: set[tuple[int, int]] = set()
        self.data_qubits: set[tuple[int, int]] = set()
        self.qubits_to_faces: dict[int, list[int]] = {}
        self.faces_to_qubits: dict[int, list[int]] = {}
        self.lattice = lattice
        self.add_qubits()
        self.H = np.zeros((len(self.ancilla_qubits), len(self.data_qubits)), dtype=int)
        self.construct_layout()
        self.compute_logical()
        self.n = len(self.qubits_to_faces)
        self.k = self.L.shape[1]

    def add_qubits(self) -> None:
        """Compute the ancilla and data qubits lists from the lattice type."""
        pass

    def construct_layout(self) -> None:
        """Construct the adjacency lists of the code from the qubits lists. Assumes add_qubits was called."""
        pass

    def compute_logical(self) -> None:
        """Compute the logical matrix L."""
        ker_hx = mod2.nullspace(self.H)  # compute the kernel basis of hx
        im_hz_transp = mod2.row_basis(self.H)  # compute the image basis of hz.T
        log_stack = np.vstack([im_hz_transp, ker_hx])
        pivots = mod2.row_echelon(log_stack.T)[3]
        log_op_indices = [i for i in range(im_hz_transp.shape[0], log_stack.shape[0]) if i in pivots]
        self.L = log_stack[log_op_indices]

    def get_syndrome(self, error: npt.NDArray[np.int_]) -> npt.NDArray[np.int_]:
        """Compute the syndrome of the error."""
        return self.H @ error % 2

    def check_if_logical_error(self, residual: npt.NDArray[np.int_]) -> bool:
        """Check if the residual is a logical error."""
        return (self.L @ residual % 2).any() is True
