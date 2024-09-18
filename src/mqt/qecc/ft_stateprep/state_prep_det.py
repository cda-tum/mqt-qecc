"""Synthesizeing optimal deterministic state preparation circuits for d=3 CSS codes."""

from __future__ import annotations

import logging

import numpy as np
import galois
import z3
from typing import TYPE_CHECKING

from .state_prep import StatePrepCircuit, _vars_to_stab, _odd_overlap, _symbolic_vector_eq, _coset_leader, iterative_search_with_timeout, _run_with_timeout, gate_optimal_verification_stabilizers, verification_stabilizers, _hook_errors


logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    import numpy.typing as npt
    Recovery = dict[int, tuple[list[npt.NDArray[np.int8]], dict[int, npt.NDArray[np.int8]]]]
    DeterministicCorrection = tuple[list[npt.NDArray[np.int8]], Recovery]
    Verification = list[npt.NDArray[np.int8]]

class DeterministicVerification:
    """
    """

    def __init__(self, nd_verification_stabs: Verification, det_correction: DeterministicCorrection | None = None, hook_corrections: list[DeterministicCorrection] | None = None):
        self.stabs = nd_verification_stabs
        self.det_correction = det_correction
        self.hook_corrections = [] if hook_corrections is None else hook_corrections

    @staticmethod
    def _num_cnots_correction(correction: DeterministicCorrection) -> int:
        if correction is None:
            return 0
        return np.sum([np.sum([np.sum(m) for m in v[0]]) for v in correction.values()])
    @staticmethod
    def _num_anc_correction(correction: DeterministicCorrection) -> int:
        if correction is None:
            return 0
        return np.sum([len(v[0]) for v in correction.values()])
    @staticmethod
    def _num_branches_correction(correction: DeterministicCorrection) -> int:
        if correction is None:
            return 0
        return sum([len(v[1]) for v in correction.values()]) 

    def is_flagged(self, stab_idx : int) -> bool:
        if self.hook_corrections[stab_idx] == None:
            return False
        return True
    def are_flagged(self) -> list[bool]:
        return [self.is_flagged(i) for i in range(len(self.stabs))]


    # Statistics methods
    def num_ancillae_verification(self) -> int:
        return len(self.stabs)
    def num_cnots_verification(self) -> int:
        return np.sum([np.sum(m) for m in self.stabs])

    def num_ancillae_correction(self) -> int:
        return self._num_anc_correction(self.det_correction)

    def num_cnots_correction(self) -> int:
        return self._num_cnots_correction(self.det_correction)

    def num_ancillae_hooks(self) -> int:
        return len([v for v in self.hook_corrections if v is not None])
    def num_cnots_hooks(self) -> int:
        return len(self.hook_corrections) * 2
    def num_ancillae_hook_corrections(self) -> int:
        return sum([self._num_anc_correction(c) for c in self.hook_corrections])
    def num_cnots_hook_corrections(self) -> int:
        return sum([self._num_cnots_correction(c) for c in self.hook_corrections])
    
    def num_ancillae_total(self) -> int:
        return self.num_ancillae_verification() + self.num_ancillae_correction() + self.num_ancillae_hooks() + self.num_ancillae_hook_corrections()
    def num_cnots_total(self) -> int:
        return self.num_cnots_verification() + self.num_cnots_correction() + self.num_cnots_hooks() + self.num_cnots_hook_corrections()
    
    def num_branches_det_correction(self) -> int:
        return self._num_branches_correction(self.det_correction)
    def num_branches_hook_corrections(self) -> int:
        return sum([self._num_branches_correction(c) for c in self.hook_corrections])
    def num_branches_total(self) -> int:
        return self.num_branches_det_correction() + self.num_branches_hook_corrections()


class DeterministicVerificationHelper:
    """
    Class for storing deterministic verification stabilizers and corrections. 
    The methods can be used to compute deterministic verifications for X or Z errors seperatly for d=3 CSS codes. For d=4 it also supports the computation of the deterministic hook error correction.
    """

    def __init__(self, state_prep: StatePrepCircuit, full_fault_tolerance: bool = True):
        self.state_prep = state_prep
        self.code = state_prep.code
        self.full_fault_tolerance = full_fault_tolerance
        assert self.code.distance < 5, "Only d=3 and d=4 codes are supported."
        self.num_qubits = self.code.n

        self._nd_layers = [[], []]

    def compute_nd_stabs(self, min_timeout: int = 1, max_timeout: int = 3600, max_ancilla: int | None = None, compute_all_solutions: bool = False):
        """
        Returns the gate non-deterministic verification stabilizers for X and Z errors. 
        """
        for idx, x_errors in enumerate([self.state_prep.zero_state, not self.state_prep.zero_state]):
            logger.info(f"Computing non-deterministic verification stabilizers for layer {idx + 1} / 2.")
            stabs = gate_optimal_verification_stabilizers(self.state_prep, x_errors=x_errors, min_timeout=min_timeout, max_timeout=max_timeout, max_ancillas=max_ancilla, return_all_solutions=compute_all_solutions)[0]
            if not compute_all_solutions:
                verify = DeterministicVerification(stabs)
                self._nd_layers[idx] = [verify]
            else:
                self._nd_layers[idx] = [DeterministicVerification(stabs) for stabs in stabs]
    
    def compute_det_corrections(self, min_timeout: int = 1, max_timeout: int = 3600, max_ancilla: int | None = None):
        """
        Returns the deterministic verification stabilizers for the first layer of non-deterministic verification stabilizers.
        """
        for layer_idx in range(2):
            logger.info(f"Computing deterministic verification for layer {layer_idx + 1} / 2.")
            for verify_idx, verify in enumerate(self._nd_layers[layer_idx]):
                logger.info(f"Computing deterministic verification for non-det verification {verify_idx + 1} / {len(self._nd_layers[layer_idx])}.")
                self._nd_layers[layer_idx][verify_idx].det_correction = deterministic_correction(self.state_prep, verify.stabs, min_timeout=min_timeout, max_timeout=max_timeout, max_ancillas=max_ancilla, zero_state=self.state_prep.zero_state)

    def compute_hook_corrections(self, min_timeout: int = 1, max_timeout: int = 3600, max_ancilla: int | None = None):
        """
        Computes the additional stabilizers to measure with corresponding corrections for the hook errors of each stabilizer measurement in layer 2.
        """
        for layer_idx in range(2):
            logger.info(f"Computing deterministic verification for hook errors of layer {layer_idx + 1} / 2.")
            for verify_idx, verify in enumerate(self._nd_layers[layer_idx]):
                if verify.stabs == []:
                    self._nd_layers[layer_idx][verify_idx].hook_corrections = []
                    continue
                for stab in verify.stabs:
                    hook_errors = _hook_errors([stab])

                    # check if the hook error is trivial
                    errors_trivial = []
                    code_stabs = np.vstack((self.code.Hz, self.code.Lz)) if self.state_prep.zero_state else np.vstack((self.code.Hx, self.code.Lx))
                    rank = np.linalg.matrix_rank(code_stabs)
                    GF = galois.GF(2)
                    for error in hook_errors:
                        single_qubit_deviation = [(error + np.eye(self.num_qubits, dtype=np.int8)[i]) % 2 for i in range(self.num_qubits)]
                        stabs_plus_single_qubit = [GF(np.vstack((code_stabs, single_qubit_deviation[i]))) for i in range(self.num_qubits)]
                        trivial = any([np.linalg.matrix_rank(m.row_reduce()) == rank for m in stabs_plus_single_qubit])
                        errors_trivial.append(trivial)
                    if all(errors_trivial):
                        self._nd_layers[layer_idx][verify_idx].hook_corrections.append(None)
                        continue
                    # TODO only keep non-trivial errors
                    # hook_errors = hook_errors[np.logical_not(errors_trivial)]

                    # hook errors are non-trivial
                    # add case of error on hook ancilla
                    hook_errors = np.vstack((hook_errors, np.zeros(self.num_qubits, dtype=np.int8)))
                    self._nd_layers[layer_idx][verify_idx].hook_corrections.append({1: deterministic_correction_single_outcome(self.state_prep, hook_errors, min_timeout=min_timeout, max_timeout=max_timeout, max_ancillas=max_ancilla, zero_state= not self.state_prep.zero_state)})
    
    def get_solution(self, min_timeout: int = 1, max_timeout: int = 3600, max_ancilla: int | None = None) -> tuple[tuple[Verification,DeterministicCorrection], tuple[Verification,DeterministicCorrection], list[DeterministicCorrection]]:
        """
        Returns a tuple representing the first layer, second layer and hook error corrections.
        """
        if max_ancilla is None:
            max_ancilla = self.code.Hx.shape[0] + self.code.Hz.shape[0]

        self.compute_nd_stabs(min_timeout=min_timeout, max_timeout=max_timeout, max_ancilla=max_ancilla)
        self.compute_hook_corrections(min_timeout=min_timeout, max_timeout=max_timeout, max_ancilla=max_ancilla)
        self.compute_det_corrections(min_timeout=min_timeout, max_timeout=max_timeout, max_ancilla=max_ancilla)
        return self._nd_layers[0][0], self._nd_layers[1][0]

    def _filter_nd_stabs(self):
        """
        Only keep the best non-deterministic verification stabilizers with minimal number of ancillae and CNOTs.
        """
        for layer_idx in range(2):
            # get best numbers
            best_num_anc = int(1e6)
            best_num_cnots = int(1e6)
            best_case_indices = []
            for idx_verify, verify in enumerate(self._nd_layers[layer_idx]):
                num_anc = verify.num_ancillae_verification() + verify.num_ancillae_hook_corrections()
                num_cnot = verify.num_cnots_verification() + verify.num_cnots_hook_corrections()
                if best_num_anc > num_anc or (best_num_anc == num_anc and best_num_cnots > num_cnot):
                    best_num_anc = num_anc
                    best_num_cnots = num_cnot
                    best_case_indices = [idx_verify]
                elif best_num_anc == num_anc and best_num_cnots == num_cnot:
                    best_case_indices.append(idx_verify)
            # filter out all but the best case
            self._nd_layers[layer_idx] = [self._nd_layers[layer_idx][idx] for idx in best_case_indices]

    def get_global_solution(self, min_timeout: int = 1, max_timeout: int = 3600, max_ancilla: int | None = None) -> tuple[Verification, Verification]:
        """
        Returns the optimal non-deterministic verification stabilizers for the first and second layer regarding the number of ancillae and CNOTs.
        """

        if max_ancilla is None:
            max_ancilla = self.code.Hx.shape[0] + self.code.Hz.shape[0]

        self.compute_nd_stabs(min_timeout=min_timeout, max_timeout=max_timeout, max_ancilla=max_ancilla, compute_all_solutions=True)
        self.compute_hook_corrections(min_timeout=min_timeout, max_timeout=max_timeout, max_ancilla=max_ancilla)
        self._filter_nd_stabs()
        self.compute_det_corrections(min_timeout=min_timeout, max_timeout=max_timeout, max_ancilla=max_ancilla)

        best_stab_indices = [0, 0]
        for layer_idx in range(2):
            best_num_anc = 3 * max_ancilla
            best_num_cnots = 3 * max_ancilla * self.num_qubits
            for idx_verify, verify in enumerate(self._nd_layers[layer_idx]):
                num_anc = verify.num_ancillae_total()
                num_cnots = verify.num_cnots_total()
                if best_num_anc > num_anc or (best_num_anc == num_anc and best_num_cnots > num_cnots):
                    best_num_anc = num_anc
                    best_num_cnots = num_cnots
                    best_stab_indices[layer_idx] = idx_verify
        if len(self._nd_layers[1]) == 0:
            return self._nd_layers[0][best_stab_indices[0]], DeterministicVerification([])
        return self._nd_layers[0][best_stab_indices[0]], self._nd_layers[1][best_stab_indices[1]]
                


def global_deterministic_verification(
        sp_circ: StatePrepCircuit,
        min_timeout: int = 1,
        max_timeout: int = 3600,
        max_ancillas: int | None = None,
        zero_state: bool = True) -> DeterministicCorrection:
        """
        Returns the gate and ancilla optimal deterministic verification circuit for a given state preparation circuit.

        First, the optimal non-deterministic verification stabilizers are computed.
        Then, the optimal deterministic verification stabilizers for each non-deterministic verification outcome are computed and the best overall solution is returned.

        Args:
            sp_circ: The state preparation circuit.
            min_timeout: Minimum timeout for the z3 solver.
            max_timeout: Maximum timeout for the z3 solver.
            max_ancillas: Maximum number of ancillae to use.
            optimize_cnots: Whether to optimize the number of CNOTs, otherwise the number of ancillae is optimized.
            zero_state: Whether to optimize for the zero state.
        """
        num_qubits = sp_circ.code.n


        # get optimal number of ancillae for non-deterministic verification
        sp_circ.max_errors = 1
        opt_nd_verif_stabs = gate_optimal_verification_stabilizers(sp_circ, x_errors=zero_state, min_timeout=min_timeout, max_timeout=max_timeout, max_ancillas=max_ancillas)[0]
        num_nd_stabs = len(opt_nd_verif_stabs)
        num_cnots = np.sum([np.sum(m) for m in opt_nd_verif_stabs])

        logger.info(f"Found optimal non-det verification stabilizers with {num_nd_stabs} stabilizers and {num_cnots} CNOTs.\n Computing all possible solutions...")

        # get the fault set
        fault_set = sp_circ.compute_fault_set(1, x_errors=zero_state)
        # get all solutions with the optimal number of ancillae and cnots
        all_nd_verify_stabs = verification_stabilizers(sp_circ,
                                                       fault_set,
                                                       num_anc=num_nd_stabs,
                                                       num_cnots=num_cnots,
                                                       x_errors=zero_state,
                                                       return_all_solutions=True)
        logger.info(f"Found {len(all_nd_verify_stabs)} optimal non-det solutions.")
        # compute det verification for each solution and return the best one
        best_nd_verify = None
        best_det_verify = None
        best_det_verify_num_anc = max_ancillas
        best_det_verify_num_cnots = max_ancillas * num_qubits

        for idx, nd_verify_stabs in enumerate(all_nd_verify_stabs):
            logger.info(f"Computing deterministic verification for non-det solution {idx+1}/{len(all_nd_verify_stabs)}")
            det_verify = deterministic_correction(sp_circ=sp_circ, nd_d3_verification_stabilizers=nd_verify_stabs,
                                                            min_timeout=min_timeout,max_timeout=max_timeout,max_ancillas=best_det_verify_num_anc,zero_state=zero_state)
            det_num_anc = np.sum([len(v[0]) for v in det_verify.values()])
            det_num_cnots = np.sum([np.sum([np.sum(m) for m in v[0]]) for v in det_verify.values()])
            logger.info(f"Found deterministic verification with {det_num_anc} ancillae and {det_num_cnots} CNOTs.")
            if best_det_verify_num_anc > det_num_anc or (best_det_verify_num_anc == det_num_anc and best_det_verify_num_cnots > det_num_cnots):
                best_nd_verify = nd_verify_stabs
                best_det_verify = det_verify
                best_det_verify_num_anc = det_num_anc
                best_det_verify_num_cnots = det_num_cnots
        logger.info(f"Optimal deterministic verification with {best_det_verify_num_anc} ancillae and {best_det_verify_num_cnots} CNOTs.")
        return best_nd_verify, best_det_verify 


def deterministic_correction(
        sp_circ: StatePrepCircuit,
        nd_d3_verification_stabilizers: list[npt.NDArray[np.int8]],
        min_timeout: int = 1,
        max_timeout: int = 3600,
        max_ancillas: int | None = None,
        zero_state: bool = True,
        additional_faults: npt.NDArray[np.int8] | None = None) -> DeterministicCorrection:
        """
        Returns an gate and ancilla optimal deterministic verification circuit for a given state preparation circuit and 
        non-deterministic verification stabilizers.
        It computes the corresponding fault set and then solves the problem if finding optimal deterministic verification
        stabilizers for each non-deterministic verification outcome separately.
        """
        num_nd_stabs = len(nd_d3_verification_stabilizers)
        num_qubits = sp_circ.code.n
        if max_ancillas is None:
            max_ancillas = sp_circ.code.Hx.shape[0] + sp_circ.code.Hz.shape[0]

        # get the fault set
        if additional_faults is not None:
            fault_set = sp_circ.combine_faults(additional_faults=additional_faults, x_errors=zero_state)
        else:
            fault_set = sp_circ.compute_fault_set(1, x_errors=zero_state)

        det_verify = dict()
        for verify_outcome_int in range(1, 2**num_nd_stabs):
            verify_outcome = _int_to_int8_array(verify_outcome_int, num_nd_stabs)
            logger.info(f"Computing deterministic verification for non-det outcome {verify_outcome}: {verify_outcome_int}/{2**num_nd_stabs-1}")

            # only consider errors that triggered the verification pattern
            errors_filtered = np.array([error for error in fault_set if np.array_equal(verify_outcome, [np.sum(m * error) % 2 for m in nd_d3_verification_stabilizers])])

            # append single-qubit errors that could have triggered the verification pattern
            for qubit in range(num_qubits):
                # compute error pattern of single-qubit error on qubit i
                error_pattern = [np.sum(m * np.eye(num_qubits, dtype=np.int8)[qubit]) % 2 for m in nd_d3_verification_stabilizers]
                for i in range(num_nd_stabs):
                    if np.array_equal(verify_outcome, error_pattern):
                        # if not already in the fault set
                        if not np.any(np.all(errors_filtered == np.eye(num_qubits, dtype=np.int8)[qubit], axis=1)):
                            errors_filtered = np.vstack((errors_filtered, np.eye(num_qubits, dtype=np.int8)[qubit]))
                    else:
                        error_pattern[i] = 0

            # add the no-error case for the error beeing on one of the verification ancillae
            if np.sum(verify_outcome) == 1:
                errors_filtered = np.vstack((errors_filtered, np.zeros(num_qubits, dtype=np.int8)))
            # case of no errors or only one error is trivial
            if errors_filtered.shape[0] == 0:
                det_verify[verify_outcome_int] = np.zeros((num_qubits, 0), dtype=np.int8), {0: np.zeros(num_qubits, dtype=np.int8), 1: np.zeros(num_qubits, dtype=np.int8)}
            elif errors_filtered.shape[0] == 1:
                det_verify[verify_outcome_int] = ([np.zeros(num_qubits, dtype=np.int8)], {0: errors_filtered[0], 1: errors_filtered[0]})
            else:
                det_verify[verify_outcome_int] = deterministic_correction_single_outcome(sp_circ, errors_filtered, min_timeout, max_timeout, max_ancillas, zero_state)
        return det_verify

def deterministic_correction_single_outcome(
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
        return correction_stabilizers(sp_circ, fault_set, num_anc, num_anc*num_qubits, x_errors=zero_state)

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
        return correction_stabilizers(sp_circ, fault_set, num_anc, num_cnots, x_errors=zero_state)

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



def correction_stabilizers(
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
