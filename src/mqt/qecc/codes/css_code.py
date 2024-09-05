"""Class for representing quantum error correction codes."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ldpc import mod2

from .stabilizer_code import StabilizerCode

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt


class CSSCode(StabilizerCode):
    """A class for representing CSS codes."""

    def __init__(
        self,
        distance: int,
        Hx: npt.NDArray[np.int8] | None = None,  # noqa: N803
        Hz: npt.NDArray[np.int8] | None = None,  # noqa: N803
        x_distance: int | None = None,
        z_distance: int | None = None,
    ) -> None:
        """Initialize the code."""
        self._check_valid_check_matrices(Hx, Hz)

        if Hx is None:
            assert Hz is not None
            self.n = Hz.shape[1]
            self.Hx = np.zeros((0, self.n), dtype=np.int8)
        else:
            self.Hx = Hx
        if Hz is None:
            assert Hx is not None
            self.n = Hx.shape[1]
            self.Hz = np.zeros((0, self.n), dtype=np.int8)
        else:
            self.Hz = Hz

        z_padding = np.zeros(self.Hx.shape, dtype=np.int8)
        x_padding = np.zeros(self.Hz.shape, dtype=np.int8)

        x_padded = np.hstack([self.Hx, z_padding])
        z_padded = np.hstack([x_padding, self.Hz])
        phases = np.zeros((x_padded.shape[0] + z_padded.shape[0], 1), dtype=np.int8)
        super().__init__(np.hstack((np.vstack((x_padded, z_padded)), phases)), distance)

        self.distance = distance
        self.x_distance = x_distance if x_distance is not None else distance
        self.z_distance = z_distance if z_distance is not None else distance

        if self.x_distance < self.distance or self.z_distance < self.distance:
            msg = "The x and z distances must be greater than or equal to the distance"
            raise InvalidCSSCodeError(msg)

        self.Lx = CSSCode._compute_logical(self.Hz, self.Hx)
        self.Lz = CSSCode._compute_logical(self.Hx, self.Hz)

    def x_checks_as_pauli_strings(self) -> list[str]:
        """Return the x checks as Pauli strings."""
        return ["".join("X" if bit == 1 else "I" for bit in row) for row in self.Hx]

    def z_checks_as_pauli_strings(self) -> list[str]:
        """Return the z checks as Pauli strings."""
        return ["".join("Z" if bit == 1 else "I" for bit in row) for row in self.Hz]

    def x_logicals_as_pauli_strings(self) -> list[str]:
        """Return the x logicals as a Pauli strings."""
        return ["".join("X" if bit == 1 else "I" for bit in row) for row in self.Lx]

    def z_logicals_as_pauli_strings(self) -> list[str]:
        """Return the z logicals as Pauli strings."""
        return ["".join("Z" if bit == 1 else "I" for bit in row) for row in self.Lz]

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
        if self.Hx is None:
            return np.empty((0, error.shape[0]), dtype=np.int8)
        return self.Hx @ error % 2

    def get_z_syndrome(self, error: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
        """Compute the z syndrome of the error."""
        if self.Hz is None:
            return np.empty((0, error.shape[0]), dtype=np.int8)
        return self.Hz @ error % 2

    def check_if_logical_x_error(self, residual: npt.NDArray[np.int8]) -> bool:
        """Check if the residual is a logical error."""
        return bool((self.Lz @ residual % 2 == 1).any())

    def check_if_x_stabilizer(self, pauli: npt.NDArray[np.int8]) -> bool:
        """Check if the Pauli is a stabilizer."""
        assert self.Hx is not None
        return bool(mod2.rank(np.vstack((self.Hx, pauli))) == mod2.rank(self.Hx))

    def check_if_logical_z_error(self, residual: npt.NDArray[np.int8]) -> bool:
        """Check if the residual is a logical error."""
        return bool((self.Lx @ residual % 2 == 1).any())

    def check_if_z_stabilizer(self, pauli: npt.NDArray[np.int8]) -> bool:
        """Check if the Pauli is a stabilizer."""
        assert self.Hz is not None
        return bool(mod2.rank(np.vstack((self.Hz, pauli))) == mod2.rank(self.Hz))

    def stabilizer_eq_x_error(self, error_1: npt.NDArray[np.int8], error_2: npt.NDArray[np.int8]) -> bool:
        """Check if two X errors are in the same coset."""
        if self.Hx is None:
            return bool(np.array_equal(error_1, error_2))
        m1 = np.vstack([self.Hx, error_1])
        m2 = np.vstack([self.Hx, error_2])
        m3 = np.vstack([self.Hx, error_1, error_2])
        return bool(mod2.rank(m1) == mod2.rank(m2) == mod2.rank(m3))

    def stabilizer_eq_z_error(self, error_1: npt.NDArray[np.int8], error_2: npt.NDArray[np.int8]) -> bool:
        """Check if two Z errors are in the same coset."""
        if self.Hz is None:
            return bool(np.array_equal(error_1, error_2))
        m1 = np.vstack([self.Hz, error_1])
        m2 = np.vstack([self.Hz, error_2])
        m3 = np.vstack([self.Hz, error_1, error_2])
        return bool(mod2.rank(m1) == mod2.rank(m2) == mod2.rank(m3))

    def is_self_dual(self) -> bool:
        """Check if the code is self-dual."""
        if self.Hx is None or self.Hz is None:
            return False
        return bool(
            self.Hx.shape[0] == self.Hz.shape[0] and mod2.rank(self.Hx) == mod2.rank(np.vstack([self.Hx, self.Hz]))
        )

    @staticmethod
    def _check_valid_check_matrices(Hx: npt.NDArray[np.int8] | None, Hz: npt.NDArray[np.int8] | None) -> None:  # noqa: N803
        """Check if the code is a valid CSS code."""
        if Hx is None and Hz is None:
            msg = "At least one of the check matrices must be provided"
            raise InvalidCSSCodeError(msg)

        if Hx is not None and Hz is not None:
            if Hx.shape[1] != Hz.shape[1]:
                msg = "Check matrices must have the same number of columns"
                raise InvalidCSSCodeError(msg)
            if np.any(Hx @ Hz.T % 2 != 0):
                msg = "The check matrices must be orthogonal"
                raise InvalidCSSCodeError(msg)

    @staticmethod
    def from_code_name(code_name: str, distance: int | None = None) -> CSSCode:
        r"""Return CSSCode object for a known code.

        The following codes are supported:
        - [[7, 1, 3]] Steane (\"Steane\")
        - [[15, 1, 3]] tetrahedral code (\"Tetrahedral\")
        - [[15, 7, 3]] Hamming code (\"Hamming\")
        - [[9, 1, 3]] Shore code (\"Shor\")
        - [[12, 2, 4]] Carbon Code (\"Carbon\")
        - [[9, 1, 3]] rotated surface code (\"Surface, 3\"), also default when no distance is given
        - [[25, 1, 5]] rotated surface code (\"Surface, 5\")
        - [[23, 1, 7]] golay code (\"Golay\")

        Args:
            code_name: The name of the code.
            distance: The distance of the code.
        """
        prefix = (Path(__file__) / "../").resolve()
        paths = {
            "steane": prefix / "steane/",
            "tetrahedral": prefix / "tetrahedral/",
            "hamming": prefix / "hamming/",
            "shor": prefix / "shor/",
            "surface_3": prefix / "rotated_surface_d3/",
            "surface_5": prefix / "rotated_surface_d5/",
            "golay": prefix / "golay/",
            "carbon": prefix / "carbon/",
        }

        distances = {
            "steane": (3, 3),
            "tetrahedral": (7, 3),
            "hamming": (3, 3),
            "shor": (3, 3),
            "golay": (7, 7),
            "surface_3": (3, 3),
            "surface_5": (5, 5),
            "carbon": (4, 4),
        }  # X, Z distances

        code_name = code_name.lower()
        if code_name == "surface":
            if distance is None:
                distance = 3
            code_name += f"_{distance}"

        if code_name in paths:
            hx = np.load(paths[code_name] / "hx.npy")
            hz = np.load(paths[code_name] / "hz.npy")

            if code_name in distances:
                x_distance, z_distance = distances[code_name]
                distance = min(x_distance, z_distance)
                return CSSCode(distance, hx, hz, x_distance=x_distance, z_distance=z_distance)

            if distance is None:
                msg = f"Distance is not specified for {code_name}"
                raise InvalidCSSCodeError(msg)
        msg = f"Unknown code name: {code_name}"
        raise InvalidCSSCodeError(msg)


class InvalidCSSCodeError(ValueError):
    """Raised when the CSS code is invalid."""
