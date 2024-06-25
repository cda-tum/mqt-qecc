"""Class for representing quantum error correction codes."""

from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np

from ldpc import mod2

try:
    from importlib import resources as impresources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as impresources

from . import sample_codes  # relative-import the *package* containing the templates


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
        self.Lx = CSSCode._compute_logical(self.Hz, self.Hx)
        self.Lz = CSSCode._compute_logical(self.Hx, self.Hz)

    def __hash__(self) -> int:
        """Compute a hash for the CSS code."""
        return hash(self.Hx.tobytes() ^ self.Hz.tobytes())

    def __eq__(self, other: object) -> bool:
        """Check if two CSS codes are equal."""
        if not isinstance(other, CSSCode):
            return NotImplemented
        return mod2.rank(self.Hx) == mod2.rank(np.vstack([self.Hx, other.Hx])) and \
               mod2.rank(self.Hz) == mod2.rank(np.vstack([self.Hz, other.Hz]))

    @staticmethod
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

    def check_if_logical_x_error(self, residual: npt.NDArray[np.int_]) -> bool:
        """Check if the residual is a logical error."""
        return (self.Lz @ residual % 2).any() is True

    def check_if_logical_z_error(self, residual: npt.NDArray[np.int_]) -> bool:
        """Check if the residual is a logical error."""
        return (self.Lx @ residual % 2).any() is True
    
    def is_self_dual(self) -> bool:
        """Check if the code is self-dual."""
        return self.Hx.shape[0] == self.Hz.shape[0] and mod2.rank(self.Hx) == mod2.rank(np.vstack([self.Hx, self.Hz]))

    def stabs_as_pauli_strings(self) -> tuple[list[str], list[str]]:
        """Return the stabilizers as Pauli strings."""
        return ["".join(["I" if x == 0 else "X" for x in row]) for row in self.Hx], \
               ["".join(["I" if z == 0 else "Z" for z in row]) for row in self.Hz]
    
    @staticmethod
    def from_code_name(code_name: str, distance: int=None) -> CSSCode:
        """Return CSSCode object for a known code.

        The following codes are supported:
        - [[7, 1, 3]] Steane (\"Steane\")
        - [[15, 1, 3]] tetrahedral code (\"Tetrahedral\")
        - [[15, 7, 3]] Hamming code (\"Hamming\")
        - [[9, 1, 3]] Shore code (\"Shor\")
        - [[9, 1, 3]] rotated surface code (\"Surface, 3\")
        - [[25, 1, 5]] rotated surface code (\"Surface, 5\")
        - [[17, 1, 5]] 4,8,8 color code (\"CC_4_8_8, 5\")
        - [[23, 1, 7]] golay code (\"Golay\")
        - 6,6,6 color code for arbitrary distances (\"CC_6_6_6, d\")
        - [[225, 9, 4]] hypergraph product code (\"HPG, 4\")


        Args:
            code_name: The name of the code.
            distance: The distance of the code.
        """
        prefix = impresources.files(sample_codes)
        paths = {
            "steane": prefix / "steane/",
            "tetrahedral": prefix / "tetrahedral/",
            "hamming": prefix / "hamming/",
            "shor": prefix / "shor/",
            "surface_3": prefix / "rotated_surface_d3/",
            "surface_5": prefix / "rotated_surface_d5/",
            "cc_4_8_8_5": prefix / "cc_4_8_8_d5/",
            "golay": prefix / "golay/",
            "hpg": prefix / "hpg_225_9_4/"
        }

        distances = {
            "steane": 3,
            "tetrahedral": 3,
            "hamming": 3,
            "shor": 3,
            "cc_4_8_8 5": 5,
            "golay": 7,
            "hpg": 4
        }

        code_name = code_name.lower()
        if code_name == "surface" or code_name == "cc_4_8_8":
            code_name = code_name + "_%d" % distance


        if code_name in paths:
            hx = np.load(paths[code_name] / "hx.npy")
            hz = np.load(paths[code_name] / "hz.npy")
                
            if code_name in distances:
                distance = distances[code_name]
            elif distance is None:
                raise ValueError(f"Distance is not specified for {code_name}")
            return CSSCode(distance, hx, hz)
        else:
            raise ValueError(f"Unknown code name: {code_name}")


class ClassicalCode:
    """A class for representing classical codes."""
    
    def __init__(self, distance: int, H: npt.NDArray[np.int_]):
        """Initialize the code."""
        self.distance = distance
        self.H = H
        self.n = H.shape[1]
        self.k = self.n - H.shape[0]

        
class HyperGraphProductCode(CSSCode):
    """A class for representing hypergraph product codes."""
    def __init__(self, c1: ClassicalCode, c2: ClassicalCode):
        """Initialize the code."""

        Hx = np.hstack((np.kron(c1.H.T, np.eye(c2.H.shape[0])),
                        np.kron(np.eye(c1.n), c2.H)))
        Hz = np.hstack((np.kron(np.eye(c1.H.shape[0]), c2.H.T),
                        np.kron(c1.H, np.eye(c2.n))))
        super().__init__(np.min(c1.distance, c2.distance), Hx, Hz)
        
