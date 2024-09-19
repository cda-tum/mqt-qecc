""" Simulation of deterministic state preparation circuits using qsample (https://github.com/dpwinter/qsample)"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from functools import partial
import qsample as qs
import matplotlib.pyplot as plt

from ..codes import InvalidCSSCodeError
from .simulation import LutDecoder, _support

if TYPE_CHECKING:
    import numpy.typing as npt
    from qiskit import QuantumCircuit
    from qsample import ErrorModel
    from .state_prep_det import DeterministicVerification, DeterministicCorrection

    from ..codes import CSSCode


class UnsupportedCode(ValueError):
    """Raised when the code is not supported by the simulator."""

def _support_int(array: npt.NDArray[np.int8]) -> set[int]:
    """
    Return the indices of the non-zero elements of the array.
    """
    return [int(i) for i in np.where(array)[0]]

def return_correction(outcome: int, corrections: dict, zero_state: bool) -> qs.Circuit:
    """
    Return the det correction circuit for the given outcome of the D-verification.
    """
    correction = _support_int(corrections[outcome])
    if len(correction) == 0:
        return qs.Circuit([], noisy=False)
    correction_circuit = []
    if zero_state:
        correction_circuit.append({"X": set(correction)})
    else:
        correction_circuit.append({"Z": set(correction)})
    return qs.Circuit(correction_circuit, noisy=False)

def ndv_outcome_check(outcome: int, correction_index: int, num_measurements: int, cutoff: int, flag: bool) -> bool:
    """
    Check if the outcome of the ND-verification is equal to the given outcome.
    """
    outcome_bitstring = format(outcome, f'0{num_measurements}b')
    if not flag:
        # return False if any of the flags trigged
        if '1' in outcome_bitstring[cutoff:]:
            return False 
        correction_outcome_bitstring = format(correction_index, f'0{cutoff}b')
        return outcome_bitstring[:cutoff] == correction_outcome_bitstring
    correction_outcome_bitstring = format(correction_index, f'0{num_measurements - cutoff}b')
    return outcome_bitstring[cutoff:] == correction_outcome_bitstring

# perform decoding using LUT
def decode_failure(measurements: int, code: CSSCode, decoder: LutDecoder, zero_state: bool, num_measurements: int) -> bool:
    """
    Check if the decoding failed.
    """
    # convert outcome to binary numpy int8 array
    measurements = np.array([int(x) for x in f"{measurements:0{num_measurements}b}"], dtype=np.int8) 
    num_syndromes = num_measurements - code.n
    syndrome = measurements[:num_syndromes]
    state = measurements[num_syndromes:]
    if zero_state:
        estimate = decoder.decode_x(syndrome)
        observables = code.Lz
    else:
        estimate = decoder.decode_z(syndrome)
        observables = code.Lx
    corrected = state + estimate

    # check if logical error
    if np.any(corrected @ observables.T % 2 != 0):
        return True
    return False 


class NoisyDFTStatePrepSimulator:
    """
    Class to simulate the state preparation circuits using qsample (https://github.com/dpwinter/qsample) using 
    with one of the predefined noise models.
    """

    def __init__(self,
        state_prep_circuit: QuantumCircuit,
        verifications: tuple[DeterministicVerification], 
        code: CSSCode,
        err_model: qs.ErrorModel = qs.noise.E1_1,
        zero_state: bool = True) -> None:

        if code.Hx is None or code.Hz is None:
            msg = "The code must have both X and Z checks."
            raise InvalidCSSCodeError(msg)
        
        if code.distance >= 5:
            msg = "Only distance <5 CSS codes are supported."
            raise UnsupportedCode(msg) 

        self.code = code
        self.err_model = err_model
        self.zero_state = zero_state
        self._ancilla_index = code.n

        # create LUT
        self.decoder = LutDecoder(code)
        if zero_state:
            self.decoder.generate_x_lut()
        else:
            self.decoder.generate_z_lut()

        # create protocol
        self.protocol = qs.Protocol()
        self.create_det_protocol(state_prep_circuit, verifications, code)

    def create_det_protocol(self, nd_state_prep_circuit: QuantumCircuit, verifications: tuple[DeterministicVerification], code: CSSCode) -> qs.Protocol:
        """
        Create the protocol for the noisy deterministic state preparation circuit.
        """
        circ = qiskit_to_qsample(nd_state_prep_circuit)
        self.protocol.add_node(name="PREP", circuit=circ)
        self.protocol.add_edge("START", "PREP", check='True')

        # create ND-verifications
        self.protocol.add_node(name=f"NDV_0", circuit=self._create_stab_measurement_circuit(verifications[0].stabs, are_flagged=verifications[0].are_flagged(), z_stabs=self.zero_state, noisy=True))
        self.protocol.add_node(name=f"NDV_1", circuit=self._create_stab_measurement_circuit(verifications[1].stabs, are_flagged=verifications[1].are_flagged(), z_stabs=not self.zero_state, noisy=True))

        self.protocol.add_edge("PREP", "NDV_0", check='True')


        # create decoding circuit
        matrix = code.Hz if self.zero_state else code.Hx
        decoding_circuit = self._create_stab_measurement_circuit(matrix, z_stabs=self.zero_state, noisy = False)
        if not self.zero_state:
            decoding_circuit.append({"H": set(range(self.code.n))})
        decoding_circuit.append({"measure": set(range(self.code.n))})

        self.protocol.add_node(name="DEC", circuit=decoding_circuit)

        # add layers
        self._append_verification(verifications[0], layer = "0", end_node="NDV_1", z_stabs=self.zero_state)
        self._append_verification(verifications[1], layer = "1", end_node="DEC", z_stabs=not self.zero_state)
        
        # add decoding and failure check
        self._append_decoding(code)

    def _append_verification(self, verification: DeterministicVerification, layer: str, end_node: str, z_stabs: bool) -> None:
        """
        Append the different deterministic verification circuits to the protocol and connect them to the ND-verification 
        depending on the outcome of the ND-verification.
        """
        # # case of no errors detected
        self.protocol.add_edge(f"NDV_{layer}", end_node, check=f"NDV_{layer}[-1] == 0")
            
        num_measurements = verification.num_ancillae_verification() + verification.num_ancillae_hooks()
        num_nd_measurements = verification.num_ancillae_verification()

        for outcome, (cor_stabs, rec) in verification.det_correction.items():
            # det corrections
            self.protocol.add_node(name=f"COR_{layer}_{outcome}", circuit=self._create_stab_measurement_circuit(cor_stabs,z_stabs=z_stabs, noisy=True))
            self.protocol.check_functions[f"ndv_outcome_check_{layer}_{outcome}"] = partial(ndv_outcome_check, correction_index=outcome, num_measurements=num_measurements, cutoff=num_nd_measurements, flag=False)
            self.protocol.add_edge(f"NDV_{layer}", f"COR_{layer}_{outcome}", check=f"ndv_outcome_check_{layer}_{outcome}(NDV_{layer}[-1])")
            # recovery
            self.protocol.check_functions[f"return_correction_circuit_{layer}_{outcome}" ] = partial(return_correction, corrections=rec, zero_state=z_stabs)
            self.protocol.add_edge(f"COR_{layer}_{outcome}", f"REC_{layer}", check=f"return_correction_circuit_{layer}_{outcome}(COR_{layer}_{outcome}[-1])")
            self.protocol.add_edge(f"REC_{layer}", end_node, check='True')

        # hooks
        hooks_idx = 0
        for hook_correction, flagged in zip(verification.hook_corrections, verification.are_flagged()):
            if not flagged:
                continue
            hook_stabs, rec = hook_correction[1]
            outcome = int('0' * hooks_idx + '1' + '0' * (verification.num_ancillae_hooks() - hooks_idx - 1),2)
            self.protocol.add_node(name=f"COR_{layer}_H_{hooks_idx}", circuit=self._create_stab_measurement_circuit(hook_stabs,z_stabs=not z_stabs, noisy=True))
            self.protocol.check_functions[f"ndv_outcome_check_h_{layer}_{hooks_idx}"] = partial(ndv_outcome_check, correction_index=outcome, num_measurements=num_measurements, cutoff=num_nd_measurements, flag=True)
            self.protocol.add_edge(f"NDV_{layer}", f"COR_{layer}_H_{hooks_idx}", check=f"ndv_outcome_check_h_{layer}_{hooks_idx}(NDV_{layer}[-1])")

            # recovery
            self.protocol.check_functions[f"return_correction_circuit_{layer}_hook_{hooks_idx}" ] = partial(return_correction, corrections=rec, zero_state=not z_stabs)
            self.protocol.add_edge(f"COR_{layer}_H_{hooks_idx}", f"REC_{layer}", check=f"return_correction_circuit_{layer}_hook_{hooks_idx}(COR_{layer}_H_{hooks_idx}[-1])")
            self.protocol.add_edge(f"REC_{layer}", end_node, check='True')
            hooks_idx += 1

    def _append_decoding(self, code: CSSCode) -> None:
        """
        Append the decoding circuit to the end of the protocol. 
        This consists of measuring the stabilizers of the code and using the LUT decoder to correct the errors.
        If the correction is not successful the protocol terminates with FAIL.
        """
        assert code.Hx is not None
        assert code.Hz is not None

        num_measurements = self.code.n
        num_measurements += len(code.Hz) if self.zero_state else len(code.Hx)

        self.protocol.check_functions["decode_failure"] = partial(decode_failure, code=code, decoder=self.decoder, zero_state=self.zero_state, 
                                                                  num_measurements=num_measurements)
        self.protocol.add_edge("DEC", "FAIL", check='decode_failure(DEC[-1])')


    def dss_logical_error_rates(
            self,
            err_params: dict[str, list[float]],
            p_max: dict[str, float],
            L : int = 3,
            shots: int = 1000,
            callbacks: list[callable] = []) -> tuple[float, int]:
        """
        Calculate the logical error rate of the deterministic state preparation circuit using direct Monte Carlo sampling.
        """
        sampler = qs.SubsetSampler(protocol=self.protocol, 
                                   simulator=qs.StabilizerSimulator,
                                   p_max=p_max,
                                   err_model=self.err_model,
                                   err_params=err_params)
        sampler.run(n_shots=shots, callbacks=callbacks)
        return sampler.stats() / self.code.k

    def mc_logical_error_rates(
            self,
            err_params: dict[str, list[float]],
            shots: int = 10000,
            callbacks: list[callable] = []
    ) -> tuple[float, int]:
        """
        Calculate the logical error rate of the deterministic state preparation circuit using direct Monte Carlo sampling.
        """
        sampler = qs.DirectSampler(protocol=self.protocol, 
                                   simulator=qs.StabilizerSimulator,
                                   err_model=self.err_model,
                                   err_params=err_params)
        sampler.run(n_shots=shots, callbacks=callbacks)
        return sampler.stats() / self.code.k

    def _create_stab_measurement_circuit(self, verification_stabilizers: list[npt.NDArray[np.int8]],  z_stabs: bool, are_flagged: list[bool]| None = None, noisy: bool = True) -> qs.Circuit:
        """
        Create the deterministic verification circuit for the given verification stabilizers using CNOT gates starting from the ancilla_index.
        """
        num_stabs = len(verification_stabilizers)
        if num_stabs == 0:
            return qs.Circuit([], noisy=False)

        if are_flagged is None:
            are_flagged = [False]*len(verification_stabilizers)
        circuit = []
        # init new ancillae
        num_ancllae = len(verification_stabilizers) + sum(are_flagged)
        circuit.append({"init": set(range(self._ancilla_index, self._ancilla_index + num_ancllae))})

        flag_ancilla_index = self._ancilla_index + len(verification_stabilizers)
        for stabilizer, flagged in zip(verification_stabilizers,are_flagged):
            stabilizer = _support_int(stabilizer)
            if not z_stabs:
                circuit.append({"H": {self._ancilla_index}})
            for qubit_idx, qubit in enumerate(stabilizer):
                # add flag
                if flagged and qubit_idx == 1:
                    if z_stabs:
                        circuit.append({"H": {flag_ancilla_index}})
                        circuit.append({"CNOT": {(flag_ancilla_index, self._ancilla_index)}})
                    else:
                        circuit.append({"CNOT": {(self._ancilla_index, flag_ancilla_index)}})
                if flagged and qubit_idx == len(stabilizer) - 1:
                    if z_stabs:
                        circuit.append({"CNOT": {(flag_ancilla_index, self._ancilla_index)}})
                        circuit.append({"H": {flag_ancilla_index}})
                    else:
                        circuit.append({"CNOT": {(self._ancilla_index, flag_ancilla_index)}})
                    flag_ancilla_index += 1

                # add stab cnots
                if z_stabs:
                    circuit.append({"CNOT": {(qubit, self._ancilla_index)}})
                else:
                    circuit.append({"CNOT": {(self._ancilla_index, qubit)}})
            if not z_stabs:
                circuit.append({"H": {self._ancilla_index}})
            # circuit.append({"measure": {self._ancilla_index}})
            self._ancilla_index += 1
        circuit.append({"measure": set(range(self._ancilla_index - num_stabs, flag_ancilla_index))})
        self._ancilla_index = flag_ancilla_index
        return qs.Circuit(circuit, noisy=noisy)

def qiskit_to_qsample(qiskit_circuit : QuantumCircuit) -> qs.Circuit:
    """
    Convert a Qiskit circuit to a qsample circuit. Only supports H, X, Y, Z, CX, and MEASURE gates.
    """

    custom_circuit = [{"init" : set(range(qiskit_circuit.num_qubits))}]
    for instruction, qargs, _ in qiskit_circuit.data:
        gate_name = instruction.name.upper()
        if gate_name == 'CX':
            gate_name = 'CNOT'
        # collect measurements to combine them at the end
        if gate_name == 'MEASURE':
            gate_name = 'measure'
        if len(qargs) == 1:
            qubits = qiskit_circuit.qubits.index(qargs[0])
        else:
            qubits = tuple(qiskit_circuit.qubits.index(q) for q in qargs)
        custom_circuit.append({gate_name: {qubits}})
    return qs.Circuit(custom_circuit, noisy=True)

