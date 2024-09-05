"""Synthesizeing optimal deterministic state preparation circuits for d=3 CSS codes."""

from __future__ import annotations

import logging

import numpy as np
import z3
from typing import TYPE_CHECKING

from qiskit import QuantumCircuit
from .state_prep import StatePrepCircuit, _vars_to_stab, _odd_overlap, _symbolic_vector_eq, _coset_leader, iterative_search_with_timeout, _run_with_timeout


logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    import numpy.typing as npt
    DeterministicVerification = dict[int, tuple[list[npt.NDArray[np.int8]], dict[int, npt.NDArray[np.int8]]]]

def optimal_deterministic_verification(
        sp_circ: StatePrepCircuit,
        nd_d3_verification_stabilizers: list[npt.NDArray[np.int8]],
        min_timeout: int = 1,
        max_timeout: int = 3600,
        max_ancillas: int | None = None,
        zero_state: bool = True) -> DeterministicVerification:
        """
        Returns an gate and ancilla optimal deterministic verification circuit for a given state preparation circuit and 
        non-deterministic verification stabilizers.
        It computes the corresponding fault set and then solves the problem if finding optimal deterministic verification
        stabilizers for each non-deterministic verification outcome separately.
        """
        num_nd_stabs = len(nd_d3_verification_stabilizers)
        num_qubits = sp_circ.code.n

        # get the fault set
        fault_set = sp_circ.compute_fault_set(1, x_errors=zero_state)
        # append single-qubit errors that could also trigger the nd verification
        for s in nd_d3_verification_stabilizers:
            for i in range(num_qubits):
                if np.any(s[i] == 1):
                    # if not already in the fault set
                    if not np.any(np.all(fault_set == np.eye(num_qubits, dtype=np.int8)[i], axis=1)):
                        fault_set = np.vstack((fault_set, np.eye(num_qubits, dtype=np.int8)[i]))
        
        det_verify = dict()
        for verify_outcome_int in range(1, 2**num_nd_stabs):
            verify_outcome = _int_to_int8_array(verify_outcome_int, num_nd_stabs)

            # only consider errors that triggered the verification pattern
            errors_filtered = np.array([error for error in fault_set if np.array_equal(verify_outcome, [np.sum(m * error) % 2 for m in nd_d3_verification_stabilizers])])
            # add the no-error case for the error beeing on one of the verification ancillae
            if np.sum(verify_outcome) == 1:
                errors_filtered = np.vstack((errors_filtered, np.zeros(num_qubits, dtype=np.int8)))
            if errors_filtered.shape[0] == 0:
                det_verify[verify_outcome_int] = np.zeros((num_qubits, 0), dtype=np.int8), {0: np.zeros(num_qubits, dtype=np.int8), 1: np.zeros(num_qubits, dtype=np.int8)}
            elif errors_filtered.shape[0] == 1:
                det_verify[verify_outcome_int] = ([np.zeros(num_qubits, dtype=np.int8)], {0: errors_filtered[0], 1: errors_filtered[0]})
            else:
                det_verify[verify_outcome_int] = optimal_deterministic_verification_single_outcome(sp_circ, errors_filtered, min_timeout, max_timeout, max_ancillas, zero_state)
        return det_verify

def optimal_deterministic_verification_single_outcome(
        sp_circ: StatePrepCircuit,
        fault_set: npt.NDArray[np.int8],
        min_timeout: int,
        max_timeout: int,
        max_ancillas: int | None = None,
        zero_state: bool = True) -> tuple[list[npt.NDArray[np.int8]], dict[int, npt.NDArray[np.int8]]]:
    """
    Returns an gate and ancilla optimal deterministic verification circuit for a given fault set.
    Geometrically increases the number of ancilla qubits until a solution is found.
    Then, first the number of ancillae is optimized and then the number of CNOTs.
    """

    num_anc = 1
    num_qubits = sp_circ.code.n

    def _func(num_anc: int):
        return deterministic_verification(sp_circ, fault_set, num_anc, num_anc*num_qubits, x_errors=zero_state)

    res = iterative_search_with_timeout(_func, num_anc, max_ancillas, min_timeout, max_timeout)

    if res is None:
        return None
    optimal_det_verify, num_anc = res
    logger.info(f"Found deterministic verification with {num_anc} ancillae.")

    while num_anc > 1:
        logger.info(f"Trying to reduce the number of ancillae to {num_anc-1}.")
        det_verify = _run_with_timeout(_func, num_anc-1, timeout=max_timeout)
        if det_verify is None or (isinstance(det_verify, str) and det_verify == "timeout"):
            break
        optimal_det_verify = det_verify
        num_anc -= 1
    logger.info(f"Optimal number of ancillae: {num_anc}.")

    # try to reduce the number of CNOTs
    def min_cnot_func(num_cnots: int):
        return deterministic_verification(sp_circ, fault_set, num_anc, num_cnots, x_errors=zero_state)

    num_cnots = 2
    while num_cnots > 1:
        # set the max number of CNOTs to the number returned by the previous step
        num_cnots = np.sum([np.sum(m) for m in optimal_det_verify[0]])

        logger.info(f"Trying to reduce the number of CNOTs to {num_cnots-1}.")
        det_verify = _run_with_timeout(min_cnot_func, num_cnots-1, timeout=max_timeout)
        if det_verify is None or (isinstance(det_verify, str) and det_verify == "timeout"):
            break
        optimal_det_verify = det_verify
        num_cnots -= 1
    logger.info(f"Optimal number of CNOTs: {num_cnots}.")
    return optimal_det_verify



def deterministic_verification(
        sp_circ: StatePrepCircuit,
        fault_set: npt.NDArray[np.int8],
        num_anc: int,
        num_cnot: int,
        x_errors: bool = True,
) -> dict[int, tuple[list[npt.NDArray[np.int8]], dict[int, npt.NDArray[np.int8]]]]:
    """
    Return deterministic verification stabilizers with corresponding corrections using z3.
    """

    gens = sp_circ.z_checks if x_errors else sp_circ.x_checks
    correction_gens = sp_circ.x_checks if x_errors else sp_circ.z_checks

    n_gens = gens.shape[0]
    n_corr_gens = correction_gens.shape[0]
    n_qubits = sp_circ.code.n
    n_errors = fault_set.shape[0]

    # Measurements are written as sums of generators
    # The variables indicate which generators are non-zero in the sum
    measurement_vars = [[z3.Bool(f"m_{anc}_{i}") for i in range(n_gens)] for anc in range(num_anc)]
    measurement_stabs = [_vars_to_stab(vars_, gens) for vars_ in measurement_vars]

    # create "stabilizer degree of freedom" variables
    free_var = [[z3.Bool(f"free_{e}_{g}") for g in range(n_corr_gens)] for e in range(n_errors)]
    free_stabs = [_vars_to_stab(vars_, correction_gens) for vars_ in free_var]

    # correction variables for each possible deterministic verification outcome
    corrections = [[z3.Bool("c_{}_{}".format(anc, i)) for i in range(n_qubits)] for anc in range(2**num_anc)]

    solver = z3.Solver()

    # for each error, the pattern is computed and the corresponding correction is applied
    error_patterns = []
    det_patterns = []
    triggereds = []
    final_errors = []
    for idx_error, error in enumerate(fault_set):
        error_pattern = [_odd_overlap(measurement, error) for measurement in measurement_stabs]
        error_patterns.append(error_pattern)
        det_patterns.append([])
        triggereds.append([])
        final_errors.append([])
        for det_pattern, correction in enumerate(corrections):
            det_pattern = _int_to_bool_array(det_pattern, num_anc)
            det_patterns[idx_error].append(det_pattern)
            # check if error triggers the pattern
            triggered = _symbolic_vector_eq(error_pattern, det_pattern)
            triggereds[idx_error].append(triggered)
            # constraint: weight(error + correction + arbitrary free stabilizer) <= 1
            final_error = [z3.Xor(correction[i] if error[i] == 0 else z3.Not(correction[i]), free_stabs[idx_error][i]) for i in range(n_qubits)]
            final_errors[idx_error].append(final_error)
            solver.add(z3.If(triggered, z3.Sum(final_error) <= 1, True))

    # assert that not too many CNOTs are used
    solver.add(z3.PbLe([(measurement[q], 1) for measurement in measurement_stabs for q in range(n_qubits)], num_cnot))
    
    if solver.check() == z3.sat:
         return _extract_measurement_and_correction(solver.model(), gens, correction_gens, n_qubits, num_anc, measurement_vars, corrections)
    return None


def _extract_measurement_and_correction(model : z3.Model,
                                        gens: list[npt.NDArray[np.int8]],
                                        correction_gens: list[npt.NDArray[np.int8]],
                                        n_qubits: int, 
                                        num_anc: int,
                                        measurement_vars: list[list[z3.BoolRef]],
                                        corrections: list[list[z3.BoolRef]]) -> tuple[list[npt.NDArray[np.int8]], dict[int, npt.NDArray[np.int8]]]:
    """ Extract deterministic verification stabilizers and corrections from sat z3 solver. """
    # get measurements
    actual_measurements = []
    for m in measurement_vars:
        v = np.zeros(len(gens[0]), dtype=np.int8)
        for g in range(len(gens)):
            if model[m[g]]:
                v += gens[g]
        actual_measurements.append(v%2)

    # get corrections for each pattern
    actual_corrections = dict()
    for outcome in range(2**num_anc):
        actual_correction = np.array([int(bool(model[corrections[outcome][i]])) for i in range(n_qubits)], dtype=np.int8)

        if np.sum(actual_correction) == 0:
            actual_corrections[outcome] = actual_correction
        else:
            actual_corrections[outcome] = _coset_leader(actual_correction, correction_gens)
    return actual_measurements, actual_corrections

def _int_to_bool_array(num: int, num_anc: int) -> npt.NDArray[np.bool_]:
    """ Convert an integer to a boolean array of length num_anc corresponding to the binary representation of the integer."""
    return np.array([bool(num & (1 << i)) for i in range(num_anc)])[::-1]

def _int_to_int8_array(num: int, n_qubits: int) -> npt.NDArray[np.int8]:
    """ Convert an integer to an int8 array of length n_qubits"""
    return np.array([int(bool(num & (1 << i))) for i in range(n_qubits)], dtype=np.int8)[::-1]
