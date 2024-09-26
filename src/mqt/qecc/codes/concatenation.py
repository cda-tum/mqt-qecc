"""Concatenated quantum codes."""

from __future__ import annotations

from typing import cast

from .pauli import Pauli
from .stabilizer_code import InvalidStabilizerCodeError, StabilizerCode
from .symplectic import SymplecticVector


class ConcatenatedCode(StabilizerCode):
    """A concatenated quantum code."""

    def __init__(self, outer_code: StabilizerCode, inner_code: StabilizerCode | list[StabilizerCode]) -> None:
        """Initialize a concatenated quantum code.

        Args:
            outer_code: The outer code.
            inner_code: The inner code. If a list of codes is provided, the qubits of the outer code are encoded by the different inner codes in the list.
        """
        self.outer_code = outer_code
        if isinstance(inner_code, list):
            self.inner_codes = inner_code
        else:
            self.inner_codes = [inner_code] * outer_code.n
        if not all(code.k == 1 for code in self.inner_codes):
            msg = "The inner codes must be stabilizer codes with a single logical qubit."
            raise InvalidStabilizerCodeError(msg)

        x_s = [code.x_logicals for code in self.inner_codes]
        z_s = [code.z_logicals for code in self.inner_codes]
        if None in x_s or None in z_s:
            msg = "The inner codes must have valid logical operators."
            raise InvalidStabilizerCodeError(msg)
        self.n = sum(code.n for code in self.inner_codes)
        generators = [self._outer_pauli_to_physical(p) for p in outer_code.generators]
        x_logicals = [self._outer_pauli_to_physical(p) for p in cast(list[Pauli], x_s)]
        z_logicals = [self._outer_pauli_to_physical(p) for p in cast(list[Pauli], z_s)]
        d = min(code.distance * outer_code.distance for code in self.inner_codes)
        super().__init__(generators, d, x_logicals, z_logicals)

    def _outer_pauli_to_physical(self, p: Pauli) -> Pauli:
        """Convert a Pauli operator on the outer code to the operator on the concatenated code.

        Args:
            p: The Pauli operator.

        Returns:
            The Pauli operator on the physical qubits.
        """
        if len(p) != self.outer_code.n:
            msg = "The Pauli operator must have the same number of qubits as the outer code."
            raise InvalidStabilizerCodeError(msg)
        concatenated = SymplecticVector.zeros(self.n)
        phase = 0
        offset = 0
        for i in range(self.outer_code.n):
            c = self.inner_codes[i]
            new_offset = offset + c.n
            assert c.x_logicals is not None
            assert c.z_logicals is not None
            if p[i] == "X":
                concatenated[offset:new_offset] = c.x_logicals[0].x_part()
                concatenated[offset + self.n : new_offset + self.n] = c.x_logicals[0].z_part()
                phase += c.x_logicals[0].phase
            elif p[i] == "Z":
                concatenated[offset:new_offset] = c.z_logicals[0].x_part()
                concatenated[offset + self.n : new_offset + self.n] = c.z_logicals[0].z_part()
                phase += c.z_logicals[0].phase

            elif p[i] == "Y":
                concatenated[offset:new_offset] = c.x_logicals[0].x_part ^ c.z_logicals[0].x_part()
                concatenated[offset + self.n : new_offset + self.n] = c.x_logicals[0].z_part ^ c.z_logicals[0].z_part()
                phase += c.x_logicals[0].phase + c.z_logicals[0].phase

            offset = new_offset
        return Pauli(concatenated, phase)


# def _valid_logicals(lst: list[StabilizerTableau | None]) -> TypeGuard[list[StabilizerTableau]]:
#     return None not in lst
