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
    from .state_prep_det import DeterministicVerification

    from ..codes import CSSCode


class UnsupportedCode(ValueError):
    """Raised when the code is not supported by the simulator."""

def _support_int(array: npt.NDArray[np.int8]) -> set[int]:
    """
    Return the indices of the non-zero elements of the array.
    """
    return set([int(i) for i in np.where(array)[0]])

def return_correction(outcome: int, corrections: dict, zero_state: bool) -> qs.Circuit:
    """
    Return the det correction circuit for the given outcome of the D-verification.
    """
    correction = _support_int(corrections[outcome])
    correction_circuit = []
    if zero_state:
        correction_circuit.append({"X": set(correction)})
    else:
        correction_circuit.append({"Z": set(correction)})
    return qs.Circuit(correction_circuit, noisy=False)

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
        nd_state_prep_circuit: QuantumCircuit,
        det_verification: DeterministicVerification, 
        code: CSSCode,
        err_model: qs.ErrorModel = qs.noise.E1_1,
        zero_state: bool = True) -> None:

        if code.Hx is None or code.Hz is None:
            msg = "The code must have both X and Z checks."
            raise InvalidCSSCodeError(msg)
        
        if code.distance > 3:
            msg = "Only distance 3 CSS codes are supported."
            raise UnsupportedCode(msg) 

        self.code = code
        self.err_model = err_model
        self.zero_state = zero_state
        self.det_verification = det_verification

        self.num_qubits = code.n
        self.num_nd_verify_qubits = nd_state_prep_circuit.num_qubits - code.n
        self.num_d_verify_qubits = []

        # create LUT
        self.decoder = LutDecoder(code)
        if zero_state:
            self.decoder.generate_x_lut()
        else:
            self.decoder.generate_z_lut()

        # create protocol
        self.protocol = qs.Protocol()
        self.create_det_protocol(nd_state_prep_circuit, det_verification, code)

    def create_det_protocol(self, nd_state_prep_circuit: QuantumCircuit, det_verification: DeterministicVerification, code: CSSCode) -> qs.Protocol:
        """
        Create the protocol for the noisy deterministic state preparation circuit.
        """
        self._append_nd_verification(nd_state_prep_circuit)
        self._append_det_verification(det_verification)
        self._append_decoding(code)

    def _append_nd_verification(self, nd_state_prep_circuit: QuantumCircuit) -> None:
        """
        Append the non-deterministic verification circuit to the start of the protocol.
        """
        circ = qiskit_to_qsample(nd_state_prep_circuit)
        self.protocol.add_node(name="NDV", circuit=circ)
        self.protocol.add_edge("START", "NDV", check='True')

    def _append_det_verification(self, det_verification: DeterministicVerification) -> None:
        """
        Append the different deterministic verification circuits to the protocol and connect them to the ND-verification 
        depending on the outcome of the ND-verification.
        """




        def _create_det_verification_circuit(verification_stabilizers: list[npt.NDArray[np.int8]], ancilla_index: int) -> qs.Circuit:
            """
            Create the deterministic verification circuit for the given verification stabilizers using CNOT gates starting from the ancilla_index.
            """
            circuit = []
            for stabilizer in verification_stabilizers:
                stabilizer = _support_int(stabilizer)
                if not self.zero_state:
                    circuit.append({"H": {ancilla_index}})
                for qubit in stabilizer:
                    if self.zero_state:
                        circuit.append({"CNOT": {(qubit, ancilla_index)}})
                    else:
                        circuit.append({"CNOT": {(ancilla_index, qubit)}})
                if not self.zero_state:
                    circuit.append({"H": {ancilla_index}})
                circuit.append({"measure": {ancilla_index}})
                ancilla_index += 1
            return qs.Circuit(circuit, noisy=True)
        
        ancilla_index = self.num_nd_verify_qubits + self.num_qubits
        for outcome, (verification_stabilizers, corrections) in det_verification.items():
            self.num_d_verify_qubits.append(len(verification_stabilizers))
            if outcome == 0:
                continue
            self.protocol.add_node(name=f"DV{outcome}", circuit=_create_det_verification_circuit(verification_stabilizers, ancilla_index))
            ancilla_index += len(verification_stabilizers)
            self.protocol.add_edge("NDV", f"DV{outcome}", check=f"NDV[-1] == {outcome}")
            self.protocol.check_functions[f"return_correction{outcome}" ] = partial(return_correction, corrections=corrections, zero_state=self.zero_state)
            # self.protocol.check_functions["return_correction" ] = return_correction
            self.protocol.add_edge(f"DV{outcome}", "COR", check=f"return_correction{outcome}(DV{outcome}[-1])")
            

    def _append_decoding(self, code: CSSCode) -> None:
        """
        Append the decoding circuit to the end of the protocol. 
        This consists of measuring the stabilizers of the code and using the LUT decoder to correct the errors.
        If the correction is not successful the protocol terminates with FAIL.
        """
        assert code.Hx is not None
        assert code.Hz is not None

        # create decoding circuit
        decoding_circuit = qs.Circuit(noisy=False)
        ancilla_index = max(self.num_d_verify_qubits) + self.num_nd_verify_qubits + self.num_qubits
        start_index = ancilla_index
        if self.zero_state:
            for check in code.Hz:
                supp = _support_int(check)
                for qubit in supp:
                    decoding_circuit.append({"CNOT": {(qubit, ancilla_index)}})
                ancilla_index += 1
        else:
            for check in code.Hx:
                supp = _support_int(check)
                decoding_circuit.append({"H": {ancilla_index}})
                for qubit in supp:
                    decoding_circuit.append({"CNOT": {(ancilla_index, qubit)}})
                decoding_circuit.append({"H": {ancilla_index}})
                ancilla_index += 1
            decoding_circuit.append({"H": set(range(self.num_qubits))})
        decoding_circuit.append({"measure": set(range(start_index, ancilla_index))})
        decoding_circuit.append({"measure": set(range(self.num_qubits))})

        self.protocol.add_node(name="DEC", circuit=decoding_circuit)
        self.protocol.add_edge("COR", "DEC", check='True')


# def decode_failure(syndrome: int, state: int, code: CSSCode, decoder: LutDecoder, zero_state: bool) -> bool:
        num_measurements = self.num_qubits
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
        return sampler.stats()

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
        return sampler.stats()


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