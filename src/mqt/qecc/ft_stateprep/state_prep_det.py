"""Synthesizeing optimal deterministic state preparation circuits for d=3 CSS codes."""

from __future__ import annotations

import logging

import numpy as np
import z3
from typing import TYPE_CHECKING, Tuple

from .state_prep import StatePrepCircuit



logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    import numpy.typing as npt
    DeterministicVerification = dict[int, tuple[list[npt.NDArray[np.int8]], dict[int, npt.NDArray[np.int8]]]]

def array_to_int(array: npt.NDArray[np.int8]) -> int:
    return int("".join(map(str, array.tolist())), 2)

# class DeterministicVerification:
#     """Class to store the additional verification circuits and the corresponding corrections for deterministic state preparation"""

#     def __init__(self) -> None:
#         """
#         Initialize the DeterministicVerification object.
#         Data is stored as a dictionary with the key being the ND-verification outcome (as int) and the value being, tuple of the verification stabilizers and the corrections.
#         The corrections themselves are stored as a dictionary with the key being the D-verifcation outcome (as int) and the value being the x/z corrections on the data qubits.
#         """
#         self.data = {0 : None}

#     def add_determinisitc_verification(self, measurement_outcome: npt.NDArray[np.int8], verification_stabilizers: npt.NDArray[np.int8], corrections: npt.NDArray[np.int8]) -> None:
#         self.data.update({array_to_int(measurement_outcome): (verification_stabilizers, corrections)})

#     def get_deterministic_verification(self, measurement_outcome: npt.NDArray[np.int8]) -> Tuple[npt.NDArray[np.int8], npt.NDArray[np.int8]]:
#         return self.data[array_to_int(measurement_outcome)]




def gate_optimal_deterministic_verification(
        sp_circ: StatePrepCircuit,
        nd_verification_stabilizers: list[npt.NDArray[np.int8]],
        min_timeout: int = 1,
        max_timeout: int = 3600,
        max_ancillas: int | None = None) -> DeterministicVerification:
        return {1: ([np.array([0,0,1,0,1,1,0])], dict({1 : np.array([0,0,0,0,0,1,0]), 0 : np.array([0,0,0,0,0,0,0])}))}
        # return {1: ([np.array([1,1,0,1,0,0,1])], dict({1 : np.array([0,0,0,0,0,0,1]), 0 : np.array([0,0,0,0,0,0,0])}))}



def determinisitc_verification(
        sp_circ: StatePrepCircuit,
        nd_verification_stabilizers: list[npt.NDArray[np.int8]],
        fault_set: npt.NDArray[np.int8],
        max_num_anc: int,
        max_num_cnot: int,
        x_errors: bool = True,
) -> DeterministicVerification:
    pass 
    



