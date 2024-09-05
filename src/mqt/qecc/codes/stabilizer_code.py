"""Class for representing general stabilizer codes."""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

if sys.version_info >= (3, 10):
    from typing import TypeAlias

    Pauli: TypeAlias = npt.NDArray[np.int8] | list[str]
else:
    from typing import Union

    from typing_extensions import TypeAlias

    Pauli: TypeAlias = Union[npt.NDArray[np.int8], list[str]]

from ldpc import mod2

if TYPE_CHECKING:
    from collections.abc import Iterable


class StabilizerCode:
    """A class for representing stabilizer codes."""

    def __init__(
        self,
        generators: npt.NDArray | list[str],
        distance: int | None = None,
        Lz: Pauli | None = None,  # noqa: N803
        Lx: Pauli | None = None,  # noqa: N803
    ) -> None:
        """Initialize the code.

        Args:
            generators: The stabilizer generators of the code. Qiskit has a reverse order of qubits in PauliList. We assume that stabilizers are ordered from left to right in ascending order of qubits.
            distance: The distance of the code.
            Lz: The logical Z-operators.
            Lx: The logical X-operators.
        """
        self._check_stabilizer_generators(generators)
        self.n = get_n_qubits_from_pauli(generators[0])
        self.generators = paulis_to_binary(generators)
        self.symplectic_matrix = self.generators[:, :-1]  # discard the phase
        self.phases = self.generators[:, -1]
        self.k = self.n - mod2.rank(self.generators)

        if distance is not None and distance <= 0:
            msg = "Distance must be a positive integer."
            raise InvalidStabilizerCodeError(msg)

        self.distance = 1 if distance is None else distance  # default distance is 1

        if Lz is not None:
            self.Lz = paulis_to_binary(Lz)
            self.Lz_symplectic = self.Lz[:, :-1]
        else:
            self.Lz = None
            self.Lz_symplectic = None

        if Lx is not None:
            self.Lx = paulis_to_binary(Lx)
            self.Lx_symplectic = self.Lx[:, :-1]
        else:
            self.Lx = None
            self.Lx_symplectic = None

        self._check_code_correct()

    def __hash__(self) -> int:
        """Compute a hash for the stabilizer code."""
        return hash(int.from_bytes(self.generators.tobytes(), sys.byteorder))

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
        return symplectic_matrix_mul(self.symplectic_matrix, pauli_to_symplectic_vec(error))

    def stabs_as_pauli_strings(self) -> list[str]:
        """Return the stabilizers as Pauli strings."""
        return [binary_to_pauli_string(s) for s in self.generators]

    def stabilizer_equivalent(self, p1: Pauli, p2: Pauli) -> bool:
        """Check if two Pauli strings are equivalent up to stabilizers of the code."""
        v1 = pauli_to_binary(p1)
        v2 = pauli_to_binary(p2)
        return bool(mod2.rank(np.vstack((self.generators, v1, v2))) == mod2.rank(np.vstack((self.generators, v1))))

    @staticmethod
    def _check_stabilizer_generators(generators: npt.NDArray[np.int8] | list[str]) -> None:
        """Check if the stabilizer generators are valid. Throws an exception if not."""
        if len(generators) == 0:
            msg = "Stabilizer code must have at least one generator."
            raise InvalidStabilizerCodeError(msg)
        if not all(len(generators[0]) == len(g) for g in generators):
            msg = "All stabilizer generators must have the same length."
            raise InvalidStabilizerCodeError(msg)

        if not isinstance(generators[0], str):
            return

        if not all(is_pauli_string(g) for g in generators):
            msg = "When providing stabilizer generators as strings, they must be valid Pauli strings."
            raise InvalidStabilizerCodeError(msg)

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

        if get_n_qubits_from_pauli(self.Lz[0]) != self.n:
            msg = "Logical operators must have the same number of qubits as the stabilizer generators."
            raise InvalidStabilizerCodeError(msg)

        if self.Lz.shape[0] > self.k:
            msg = "Number of logical Z-operators must be at most the number of logical qubits."
            raise InvalidStabilizerCodeError(msg)

        if get_n_qubits_from_pauli(self.Lx[0]) != self.n:
            msg = "Logical operators must have the same number of qubits as the stabilizer generators."
            raise InvalidStabilizerCodeError(msg)

        if self.Lx.shape[0] > self.k:
            msg = "Number of logical X-operators must be at most the number of logical qubits."
            raise InvalidStabilizerCodeError(msg)

        if not all_commute(self.Lz_symplectic, self.symplectic_matrix):
            msg = "Logical Z-operators must anti-commute with the stabilizer generators."
            raise InvalidStabilizerCodeError(msg)
        if not all_commute(self.Lx_symplectic, self.symplectic_matrix):
            msg = "Logical X-operators must commute with the stabilizer generators."
            raise InvalidStabilizerCodeError(msg)

        commutations = symplectic_matrix_product(self.Lz_symplectic, self.Lx_symplectic)
        if not np.all(np.sum(commutations, axis=1) == 1):
            msg = "Every logical X-operator must anti-commute with exactly one logical Z-operator."
            raise InvalidStabilizerCodeError(msg)


def pauli_to_binary(p: Pauli) -> npt.NDArray:
    """Convert a Pauli string to a binary array."""
    if isinstance(p, np.ndarray):
        return p

    # check if there is a sign
    phase = 0
    if p[0] in {"+", "-"}:
        phase = 0 if p[0] == "+" else 1
        p = p[1:]
    x_part = np.array([int(p == "X") for p in p])
    z_part = np.array([int(p == "Z") for p in p])
    y_part = np.array([int(p == "Y") for p in p])
    x_part += y_part
    z_part += y_part
    return np.hstack((x_part, z_part, np.array([phase])))


def paulis_to_binary(ps: Iterable[Pauli]) -> npt.NDArray:
    """Convert a list of Pauli strings to a 2d binary array."""
    return np.array([pauli_to_binary(p) for p in ps])


def binary_to_pauli_string(b: npt.NDArray) -> str:
    """Convert a binary array to a Pauli string."""
    x_part = b[: len(b) // 2]
    z_part = b[len(b) // 2 : -1]
    phase = b[-1]

    pauli = ["X" if x and not z else "Z" if z and not x else "Y" if x and z else "I" for x, z in zip(x_part, z_part)]
    return f"{'' if phase == 0 else '-'}" + "".join(pauli)


def is_pauli_string(p: str) -> bool:
    """Check if a string is a valid Pauli string."""
    return len(p) > 0 and all(c in {"I", "X", "Y", "Z"} for c in p[1:]) and p[0] in {"+", "-", "I", "X", "Y", "Z"}


def get_n_qubits_from_pauli(p: Pauli) -> int:
    """Get the number of qubits from a Pauli string."""
    if isinstance(p, np.ndarray):
        return int(p.shape[0] // 2)
    if p[0] in {"+", "-"}:
        return len(p) - 1
    return len(p)


def commute(p1: npt.NDArray[np.int8], p2: npt.NDArray[np.int8]) -> bool:
    """Check if two Paulistrings in binary representation commute."""
    return bool(symplectic_inner_product(p1, p2) == 0)


def anti_commute(p1: npt.NDArray[np.int8], p2: npt.NDArray[np.int8]) -> bool:
    """Check if two Paulistrings in binary representation anti-commute."""
    return not commute(p1, p2)


def all_commute(ps1: npt.NDArray[np.int8], ps2: npt.NDArray[np.int8]) -> bool:
    """Check if all Paulistrings in binary representation commute."""
    return bool((symplectic_matrix_product(ps1, ps2) == 0).all())


def symplectic_inner_product(p1: npt.NDArray[np.int8], p2: npt.NDArray[np.int8]) -> int:
    """Compute the symplectic inner product of two symplectic vectors."""
    n = p1.shape[0] // 2
    return int((p1[:n] @ p2[n:] + p1[n:] @ p2[:n]) % 2)


def symplectic_matrix_product(m1: npt.NDArray[np.int8], m2: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
    """Compute the symplectic matrix product of two symplectic matrices."""
    n = m1.shape[1] // 2
    return ((m1[:, :n] @ m2[:, n:].T) + (m1[:, n:] @ m2[:, :n].T)) % 2


def symplectic_matrix_mul(m: npt.NDArray[np.int8], v: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
    """Compute the symplectic matrix product of symplectic matrix with symplectic vector."""
    n = m.shape[1] // 2
    return (m[:, :n] @ v[n:] + m[:, n:] @ v[:n]) % 2


def pauli_to_symplectic_vec(p: Pauli) -> npt.NDArray:
    """Convert a Pauli string to a symplectic vector."""
    return pauli_to_binary(p)[:-1]


class InvalidStabilizerCodeError(ValueError):
    """Raised when the stabilizer code is invalid."""
