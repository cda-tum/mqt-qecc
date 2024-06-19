"""Class for representing quantum error correction codes."""

from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING

import numpy as np
from ldpc import mod2

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt

    
class CSSCode:
    """A class for representing CSS codes."""
    def __init__(self, distance: int, Hx: npt.NDArray[np.int_], Hz: npt.NDArray[np.int_]):
        """Initialize the code."""
        self.distance = distance
        
        assert Hx.shape[1] == Hz.shape[1], "Hx and Hz must have the same number of columns"
        
        self.Hx = Hx
        self.Hz = Hz
        self.n = Hx.shape[1]
        self.k = self.n - Hx.shape[0] - Hz.shape[0]
        
        self.Lx = CSSCode._compute_logical(self.Hx, self.Hz)
        self.Lz = CSSCode._compute_logical(self.Hz, self.Hx)
        
    def __hash__(self) -> int:
        """Compute a hash for the CSS code."""
        return hash(self.Hx.tobytes() ^ self.Hz.tobytes())

    def __eq__(self, other: object) -> bool:
        """Check if two CSS codes are equal."""
        if not isinstance(other, CSSCode):
            return NotImplemented
        return mod2.rank(self.Hx) == mod2.rank(np.vstack([self.Hx, other.Hx])) and \
               mod2.rank(self.Hz) == mod2.rank(np.vstack([self.Hz, other.Hz]))
    
    def _compute_logical(m1: npt.NDArray[np.int_], m2: npt.NDArray[np.int_]) -> npt.NDArray[np.int_]:
        """Compute the logical matrix L."""
        ker_m1 = mod2.nullspace(m1)  # compute the kernel basis of m1
        im_m2_transp = mod2.row_basis(m2)  # compute the image basis of m2
        log_stack = np.vstack([im_m2_transp, ker_m1])
        pivots = mod2.row_echelon(log_stack.T)[3]
        log_op_indices = [i for i in range(im_m2_transp.shape[0], log_stack.shape[0]) if i in pivots]
        return log_stack[log_op_indices]

    def get_x_syndrome(self, error: npt.NDArray[np.int_]) -> npt.NDArray[np.int_]:
        """Compute the x syndrome of the error."""
        return self.Hx @ error % 2

    def get_z_syndrome(self, error: npt.NDArray[np.int_]) -> npt.NDArray[np.int_]:
        """Compute the z syndrome of the error."""
        return self.Hz @ error % 2

    def is_self_dual(self) -> bool:
        """Check if the code is self-dual."""
        return mod2.rank(self.Hx) == mod2.rank(np.vstack([self.Hx, self.Hz]))
