"""Class for representing general stabilizer codes."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ldpc import mod2

if TYPE_CHECKING:
    import numpy.typing as npt

    from .pauli import Pauli, StabilizerTableau


class StabilizerCode:
    """A class for representing stabilizer codes."""

    def __init__(
        self,
        generators: StabilizerTableau,
        distance: int | None = None,
        Lz: StabilizerTableau | None = None,  # noqa: N803
        Lx: StabilizerTableau | None = None,  # noqa: N803
    ) -> None:
        """Initialize the code.

        Args:
            generators: The stabilizer generators of the code. We assume that stabilizers are ordered from left to right in ascending order of qubits.
            distance: The distance of the code.
            Lz: The logical Z-operators.
            Lx: The logical X-operators.
        """
        self.generators = generators
        self.n = generators.n
        self.k = self.n - mod2.rank(self.generators.tableau)

        if distance is not None and distance <= 0:
            msg = "Distance must be a positive integer."
            raise InvalidStabilizerCodeError(msg)

        self.distance = 1 if distance is None else distance  # default distance is 1

        self.Lz = None
        self.Lx = None

        if Lz is not None:
            self.Lz = Lz
        if Lx is not None:
            self.Lx = Lx

        self._check_code_correct()

    def __hash__(self) -> int:
        """Compute a hash for the stabilizer code."""
        return hash(self.generators)

    def __eq__(self, other: object) -> bool:
        """Check if two stabilizer codes are equal."""
        if not isinstance(other, StabilizerCode):
            return NotImplemented
        rnk = mod2.rank(self.generators)
        return bool(
            rnk == mod2.rank(other.generators) and rnk == mod2.rank(np.vstack((self.generators, other.generators)))
        )

    def get_syndrome(self, error: Pauli) -> npt.NDArray:
        """Compute the syndrome of the error.

        Args:
            error: The error as a pauli string or binary vector.
        """
        return self.generators.tableau @ error.symplectic

    def stabs_as_pauli_strings(self) -> list[str]:
        """Return the stabilizers as Pauli strings."""
        return [str(p) for p in self.generators]

    def stabilizer_equivalent(self, p1: Pauli, p2: Pauli) -> bool:
        """Check if two Pauli strings are equivalent up to stabilizers of the code."""
        return bool(
            mod2.rank(np.vstack((self.generators.as_matrix(), p1.as_vector(), p2.as_vector())))
            == mod2.rank(np.vstack((self.generators, p1.as_vector())))
        )

    def _check_code_correct(self) -> None:
        """Check if the code is correct. Throws an exception if not."""
        if self.Lz is not None or self.Lx is not None:
            if self.Lz is None:
                msg = "If logical X-operators are given, logical Z-operators must also be given."
                raise InvalidStabilizerCodeError(msg)
            if self.Lx is None:
                msg = "If logical Z-operators are given, logical X-operators must also be given."
                raise InvalidStabilizerCodeError(msg)

        if self.Lz is None:
            return

        if self.Lx is None:
            return

        if self.Lz.n != self.n:
            msg = "Logical operators must have the same number of qubits as the stabilizer generators."
            raise InvalidStabilizerCodeError(msg)

        if self.Lz.shape[0] > self.k:
            msg = "Number of logical Z-operators must be at most the number of logical qubits."
            raise InvalidStabilizerCodeError(msg)

        if self.Lx.n != self.n:
            msg = "Logical operators must have the same number of qubits as the stabilizer generators."
            raise InvalidStabilizerCodeError(msg)

        if self.Lx.shape[0] > self.k:
            msg = "Number of logical X-operators must be at most the number of logical qubits."
            raise InvalidStabilizerCodeError(msg)

        if not self.Lz.all_commute(self.generators):
            msg = "Logical Z-operators must anti-commute with the stabilizer generators."
            raise InvalidStabilizerCodeError(msg)
        if not self.Lx.all_commute(self.generators):
            msg = "Logical X-operators must commute with the stabilizer generators."
            raise InvalidStabilizerCodeError(msg)

        commutations = self.Lz.tableau @ self.Lx.tableau
        if not np.all(np.sum(commutations, axis=1) == 1):
            msg = "Every logical X-operator must anti-commute with exactly one logical Z-operator."
            raise InvalidStabilizerCodeError(msg)


# def commute(p1: npt.NDArray[np.int8], p2: npt.NDArray[np.int8]) -> bool:
#     """Check if two Paulistrings in binary representation commute."""
#     return bool(symplectic_inner_product(p1, p2) == 0)


# def anti_commute(p1: npt.NDArray[np.int8], p2: npt.NDArray[np.int8]) -> bool:
#     """Check if two Paulistrings in binary representation anti-commute."""
#     return not commute(p1, p2)


# def all_commute(ps1: npt.NDArray[np.int8], ps2: npt.NDArray[np.int8]) -> bool:
#     """Check if all Paulistrings in binary representation commute."""
#     return bool((symplectic_matrix_product(ps1, ps2) == 0).all())


# def symplectic_inner_product(p1: npt.NDArray[np.int8], p2: npt.NDArray[np.int8]) -> int:
#     """Compute the symplectic inner product of two symplectic vectors."""
#     n = p1.shape[0] // 2
#     return int((p1[:n] @ p2[n:] + p1[n:] @ p2[:n]) % 2)


# def symplectic_matrix_product(m1: npt.NDArray[np.int8], m2: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
#     """Compute the symplectic matrix product of two symplectic matrices."""
#     n = m1.shape[1] // 2
#     return ((m1[:, :n] @ m2[:, n:].T) + (m1[:, n:] @ m2[:, :n].T)) % 2


# def symplectic_matrix_mul(m: npt.NDArray[np.int8], v: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
#     """Compute the symplectic matrix product of symplectic matrix with symplectic vector."""
#     n = m.shape[1] // 2
#     return (m[:, :n] @ v[n:] + m[:, n:] @ v[:n]) % 2


# def pauli_to_symplectic_vec(p: Pauli) -> npt.NDArray:
#     """Convert a Pauli string to a symplectic vector."""
#     return pauli_to_binary(p)[:-1]


class InvalidStabilizerCodeError(ValueError):
    """Raised when the stabilizer code is invalid."""
