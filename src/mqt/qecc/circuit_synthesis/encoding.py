"""Methods for synthesizing encoding circuits for CSS codes."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from ..codes import InvalidCSSCodeError
from .synthesis_utils import build_css_circuit_from_list_and_checks, heuristic_gaussian_elimination

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt
    from qiskit import QuantumCircuit

    from ..codes import CSSCode

logger = logging.getLogger(__name__)


def heuristic_encoding_circuit(code: CSSCode, optimize_depth: bool = True) -> QuantumCircuit:
    """Synthesize an encoding circuit for the given CSS code using a heuristic greedy search.

    Args:
        code: The CSS code to synthesize the encoding circuit for.
        optimize_depth: Whether to optimize the depth of the circuit.

    Returns:
        The synthesized encoding circuit and the qubits that are used to encode the logical qubits.
    """
    logging.info("Starting encoding circuit synthesis.")
    checks, logicals, use_x_checks = _get_matrix_with_fewest_checks(code)
    n_checks = checks.shape[0]
    checks, cnots = heuristic_gaussian_elimination(np.vstack((checks, logicals)), parallel_elimination=optimize_depth)
    cnots = cnots[::-1]

    encoding_qubits = np.where(checks[n_checks:, :].sum(axis=0) != 0)[0]
    if use_x_checks:
        hadamards = np.where(checks[:n_checks, :].sum(axis=0) != 0)[0]
    else:
        hadamards = np.where(checks[:n_checks, :].sum(axis=0) == 0)[0]
        cnots = [(j, i) for i, j in cnots]

    hadamards = np.setdiff1d(hadamards, encoding_qubits)
    circ = build_css_circuit_from_list_and_checks(checks.shape[1], cnots, list(hadamards))
    return circ, encoding_qubits


def gate_optimal_encoding_circuit(code: CSSCode, optimize_depth: bool = True) -> QuantumCircuit:
    """Synthesize an encoding circuit for the given CSS code using the minimal number of gates.

    Args:
        code: The CSS code to synthesize the encoding circuit for.
        optimize_depth: Whether to optimize the depth of the circuit.

    Returns:
        The synthesized encoding circuit and the qubits that are used to encode the logical qubits.
    """
    logging.info("Starting optimal encoding circuit synthesis.")
    checks, logicals, _use_x_checks = _get_matrix_with_fewest_checks(code)
    checks, cnots = heuristic_gaussian_elimination(np.vstack((checks, logicals)), parallel_elimination=optimize_depth)
    cnots = cnots[::-1]


def _get_matrix_with_fewest_checks(code: CSSCode) -> tuple[npt.NDArray[np.int8], npt.NDArray[np.int8], bool]:
    """Return the stabilizer matrix with the fewest checks, the corresponding logicals and a bool indicating whether X- or Z-checks have been returned."""
    if code.Hx is None or code.Hz is None:
        msg = "The code must have both X and Z stabilizers defined."
        raise InvalidCSSCodeError(msg)

    use_x_checks = code.Hx.shape[0] < code.Hz.shape[0]
    checks = code.Hx if use_x_checks else code.Hz
    logicals = code.Lx if use_x_checks else code.Lz
    return checks, logicals, use_x_checks
