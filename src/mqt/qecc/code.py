"""Class for representing quantum error correction codes."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ldpc import mod2

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt


class CSSCode:
    """A class for representing CSS codes."""

    def __init__(
        self,
        distance: int,
        Hx: npt.NDArray[np.int8],
        Hz: npt.NDArray[np.int8],
        x_distance: int | None = None,
        z_distance: int | None = None,
    ) -> None:  # noqa: N803
        """Initialize the code."""
        self.distance = distance
        self.x_distance = x_distance if x_distance is not None else distance
        self.z_distance = z_distance if z_distance is not None else distance

        assert self.distance <= min(
            self.x_distance, self.z_distance
        ), "The distance must be less than or equal to the x and z distances"
        assert Hx.shape[1] == Hz.shape[1], "Hx and Hz must have the same number of columns"

        self.Hx = Hx
        self.Hz = Hz
        self.n = Hx.shape[1]
        self.k = self.n - Hx.shape[0] - Hz.shape[0]
        self.Lx = CSSCode._compute_logical(self.Hz, self.Hx)
        self.Lz = CSSCode._compute_logical(self.Hx, self.Hz)

    def __hash__(self) -> int:
        """Compute a hash for the CSS code."""
        return hash(int.from_bytes(self.Hx.tobytes(), sys.byteorder) ^ int.from_bytes(self.Hz.tobytes(), sys.byteorder))

    def __eq__(self, other: object) -> bool:
        """Check if two CSS codes are equal."""
        if not isinstance(other, CSSCode):
            return NotImplemented
        return bool(
            mod2.rank(self.Hx) == mod2.rank(np.vstack([self.Hx, other.Hx]))
            and mod2.rank(self.Hz) == mod2.rank(np.vstack([self.Hz, other.Hz]))
        )

    @staticmethod
    def _compute_logical(m1: npt.NDArray[np.int8], m2: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
        """Compute the logical matrix L."""
        ker_m1 = mod2.nullspace(m1)  # compute the kernel basis of m1
        im_m2_transp = mod2.row_basis(m2)  # compute the image basis of m2
        log_stack = np.vstack([im_m2_transp, ker_m1])
        pivots = mod2.row_echelon(log_stack.T)[3]
        log_op_indices = [i for i in range(im_m2_transp.shape[0], log_stack.shape[0]) if i in pivots]
        return log_stack[log_op_indices]

    def get_x_syndrome(self, error: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
        """Compute the x syndrome of the error."""
        return self.Hx @ error % 2

    def get_z_syndrome(self, error: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
        """Compute the z syndrome of the error."""
        return self.Hz @ error % 2

    def check_if_logical_x_error(self, residual: npt.NDArray[np.int8]) -> bool:
        """Check if the residual is a logical error."""
        return (self.Lz @ residual % 2).any() is True

    def check_if_logical_z_error(self, residual: npt.NDArray[np.int8]) -> bool:
        """Check if the residual is a logical error."""
        return (self.Lx @ residual % 2).any() is True

    def stabilizer_eq_x_error(self, error_1: npt.NDArray[np.int8], error_2: npt.NDArray[np.int8]) -> bool:
        """Check if two X errors are in the same coset."""
        m1 = np.vstack([self.Hx, error_1])
        m2 = np.vstack([self.Hx, error_2])
        m3 = np.vstack([self.Hx, error_1, error_2])
        return mod2.rank(m1) == mod2.rank(m2) == mod2.rank(m3)

    def stabilizer_eq_z_error(self, error_1: npt.NDArray[np.int8], error_2: npt.NDArray[np.int8]) -> bool:
        """Check if two Z errors are in the same coset."""
        m1 = np.vstack([self.Hz, error_1])
        m2 = np.vstack([self.Hz, error_2])
        m3 = np.vstack([self.Hz, error_1, error_2])
        return mod2.rank(m1) == mod2.rank(m2) == mod2.rank(m3)

    def is_self_dual(self) -> bool:
        """Check if the code is self-dual."""
        return self.Hx.shape[0] == self.Hz.shape[0] and mod2.rank(self.Hx) == mod2.rank(np.vstack([self.Hx, self.Hz]))

    def stabs_as_pauli_strings(self) -> tuple[list[str], list[str]]:
        """Return the stabilizers as Pauli strings."""
        return ["".join(["I" if x == 0 else "X" for x in row]) for row in self.Hx], [
            "".join(["I" if z == 0 else "Z" for z in row]) for row in self.Hz
        ]

    def z_logicals_as_pauli_string(self) -> str:
        """Return the logical Z operator as a Pauli string."""
        return "".join(["I" if z == 0 else "Z" for z in self.Lx[0]])

    def x_logicals_as_pauli_string(self) -> str:
        """Return the logical X operator as a Pauli string."""
        return "".join(["I" if x == 0 else "X" for x in self.Lz[0]])

    @staticmethod
    def from_code_name(code_name: str, distance: int | None = None) -> CSSCode:
        r"""Return CSSCode object for a known code.

        The following codes are supported:
        - [[7, 1, 3]] Steane (\"Steane\")
        - [[15, 1, 3]] tetrahedral code (\"Tetrahedral\")
        - [[15, 7, 3]] Hamming code (\"Hamming\")
        - [[9, 1, 3]] Shore code (\"Shor\")
        - [[9, 1, 3]] rotated surface code (\"Surface, 3\"), also default when no distance is given
        - [[25, 1, 5]] rotated surface code (\"Surface, 5\")
        - [[17, 1, 5]] 4,8,8 color code (\"CC_4_8_8\")
        - [[23, 1, 7]] golay code (\"Golay\")
        - 6,6,6 color code for arbitrary distances (\"CC_6_6_6, d\")


        Args:
            code_name: The name of the code.
            distance: The distance of the code.
        """
        prefix = (Path(__file__) / "../sample_codes/").resolve()
        paths = {
            "steane": prefix / "steane/",
            "tetrahedral": prefix / "tetrahedral/",
            "hamming": prefix / "hamming/",
            "shor": prefix / "shor/",
            "surface_3": prefix / "rotated_surface_d3/",
            "surface_5": prefix / "rotated_surface_d5/",
            "cc_4_8_8": prefix / "cc_4_8_8_d5/",
            "golay": prefix / "golay/",
        }

        distances = {
            "steane": (3, 3),
            "tetrahedral": (7, 3),
            "hamming": (3, 3),
            "shor": (3, 3),
            "cc_4_8_8": (5, 5),
            "golay": (7, 7),
            "surface_3": (3, 3),
            "surface_5": (5, 5),
        }  # X, Z distances

        code_name = code_name.lower()
        if code_name == "surface":
            if distance is None:
                distance = 3
            code_name += "_%d" % distance

        if code_name in paths:
            hx = np.load(paths[code_name] / "hx.npy")
            hz = np.load(paths[code_name] / "hz.npy")

            if code_name in distances:
                x_distance, z_distance = distances[code_name]
                distance = min(x_distance, z_distance)
            elif distance is None:
                msg = f"Distance is not specified for {code_name}"
                raise ValueError(msg)
            return CSSCode(distance, hx, hz, x_distance=x_distance, z_distance=z_distance)
        msg = f"Unknown code name: {code_name}"
        raise ValueError(msg)


class ClassicalCode:
    """A class for representing classical codes."""

    def __init__(self, distance: int, H: npt.NDArray[np.int8]) -> None:  # noqa: N803
        """Initialize the code."""
        self.distance = distance
        self.H = H
        self.n = H.shape[1]
        self.k = self.n - H.shape[0]


class HyperGraphProductCode(CSSCode):
    """A class for representing hypergraph product codes."""

    def __init__(self, c1: ClassicalCode, c2: ClassicalCode) -> None:
        """Initialize the code."""
        Hx = np.hstack((np.kron(c1.H.T, np.eye(c2.H.shape[0])), np.kron(np.eye(c1.n), c2.H)))  # noqa: N806
        Hz = np.hstack((np.kron(np.eye(c1.H.shape[0]), c2.H.T), np.kron(c1.H, np.eye(c2.n))))  # noqa: N806
        super().__init__(np.min(c1.distance, c2.distance), Hx, Hz)
