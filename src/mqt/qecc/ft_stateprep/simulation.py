"""Simulation of Non-deterministic fault tolerant state preparation."""

import stim
from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag, dag_to_circuit
import numpy as np
from collections import defaultdict

from state_prep import NDFTStatePrepCircuit
from code import CSSCode

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt
    
class NoisyNDFTStatePrepSimulator:
    """A noisy state preparation circuit."""

    def __init__(self, state_prep_circ: QuantumCircuit, code: CSSCode, p: float, zero_state: bool = False):
        self.circ = state_prep_circ
        self.num_qubits = state_prep_circ.num_qubits
        self.code = code
        self.p = p
        self.zero_state = zero_state
        self.verification_measurements = []
        self.x_measurements = []
        self.z_measurements = []
        self.data_measurements = []
        self.n_measurements = 0
        self.stim_circ = None
        self.set_p(p)
        self.decoder = LUTDecoder(code)

    def set_p(self, p: float) -> None:
        """Set the error rate."""
        self.verification_measurements = []
        self.x_measurements = []
        self.z_measurements = []
        self.data_measurements = []
        self.n_measurements = 0
        self.p = p
        self.stim_circ = self.to_stim_circ()
        self.measure_stabilizers()
        
    def to_stim_circ(self) -> stim.Circuit:
        """Convert a QuantumCircuit to a noisy STIM circuit.

        A depolarizing error model is used:
        - Single-qubit gates and idling qubits are followed by a single-qubit Pauli error with probability 2/9 p. This reflects the fact that two-qubit gates are more likely to fail.
        - Two-qubit gates are followed by a two-qubit Pauli error with probability p/15.
        - Measurements flip with a probability of 2/3 p.
        - Qubit are initialized in the -1 Eigenstate with probability 2/3 p.

        Args:
            circ: The QuantumCircuit to convert.
            p: The error rate.
        """
        initialized = [False for _ in self.circ.qubits]
        stim_circuit = stim.Circuit()
        ctrls = []

        def idle_error(used_qubits):
            for q in self.circ.qubits:
                qubit = self.circ.find_bit(q)[0]
                if initialized[qubit] and qubit not in used_qubits:
                    stim_circuit.append_operation("DEPOLARIZE1", [self.circ.find_bit(q)[0]], [2*self.p/3])

        dag = circuit_to_dag(self.circ)
        layers = dag.layers()
        used_qubits = []

        for layer in layers:
            layer_circ = dag_to_circuit(layer["graph"])

            # Apply idling errors to all qubits that were unused in the previous layer
            if len(used_qubits) > 0:
                idle_error(used_qubits)

            used_qubits = []
            for gate in layer_circ.data:
                if gate[0].name == "h":
                    ctrls.append(self.circ.find_bit(gate[1][0])[0])

                elif gate[0].name == "cx":
                    ctrl = self.circ.find_bit(gate[1][0])[0]
                    target = self.circ.find_bit(gate[1][1])[0]
                    if not initialized[ctrl]:
                        if ctrl in ctrls:
                            stim_circuit.append_operation("H", [ctrl])
                            stim_circuit.append_operation("Z_ERROR", [ctrl], [2*self.p/3])  # Wrong initialization
                        else:
                            stim_circuit.append_operation("X_ERROR", [ctrl], [2*self.p/3])  # Wrong initialization
                        initialized[ctrl] = True
                    if not initialized[target]:
                        stim_circuit.append_operation("X_ERROR", [target], [2*self.p/3])  # Wrong initialization
                        initialized[target] = True

                    stim_circuit.append_operation("CX", [ctrl, target])
                    stim_circuit.append_operation("DEPOLARIZE2", [ctrl, target], [self.p])
                    used_qubits.extend([ctrl, target])

                elif gate[0].name == "measure":
                    stim_circuit.append_operation("X_ERROR", [self.circ.find_bit(q)[0] for q in gate[1]], [2*self.p/3])
                    stim_circuit.append_operation("M", [self.circ.find_bit(q)[0] for q in gate[1]])
                    self.verification_measurements.append(self.n_measurements)
                    self.n_measurements += 1
                    used_qubits.extend([self.circ.find_bit(q)[0] for q in gate[1]])

        return stim_circuit

    def measure_stabilizers(self) -> stim.Circuit:
        """Measure the stabilizers of the code.

        An ancilla is used for each measurement.
        """
        for check in self.code.Hx:
            supp = _support(check)
            anc = self.stim_circ.num_qubits()
            self.stim_circ.append_operation("H", [anc])
            for q in supp:
                self.stim_circ.append_operation("CX", [anc, q])
            self.stim_circ.append_operation("MRX", [anc])
            self.x_measurements.append(self.n_measurements)
            self.n_measurements += 1
            
        for check in self.code.z_checks:
            supp = _support(check)
            anc = self.stim_circ.num_qubits()
            for q in supp:
                self.stim_circ.append_operation("CX", [q, anc])
            self.stim_circ.append_operation("MRZ", [anc])
            self.z_measurements.append(self.n_measurements)
            self.n_measurements += 1

    def measure_z(self) -> None:
        """Measure all data qubits in the Z basis."""
        self.data_measurements = [self.n_measurements + i for i in range(self.num_qubits)]
        self.n_measurements += self.num_qubits
        self.circuit.append_operation("MRZ", [q for q in range(self.num_qubits)])

    def measure_x(self) -> None:
        """Measure all data qubits in the X basis."""
        self.data_measurements = [self.n_measurements + i for i in range(self.num_qubits)]
        self.n_measurements += self.num_qubits
        self.circuit.append_operation("MRX", [q for q in range(self.num_qubits)])

    def logical_error_rate(self, shots=100000, shots_per_batch=100000, at_least_min_errors=True, min_errors=500):
        batch = min(shots_per_batch, shots)
        p_l = 0
        r_a = 0

        num_logical_errors = 0

        if self.zero_state:
            self.decoder.generate_x_lut()
        else:
            self.decoder.generate_z_lut()

        i = 1
        while i <= int(np.ceil(shots/batch)) or at_least_min_errors:
            num_logical_errors_batch, discarded_batch = self._simulate_batch(batch)

            p_l_batch = num_logical_errors-batch/(batch-discarded_batch)
            r_a_batch = 1-discarded_batch/batch

            # Update statistics
            num_logical_errors += num_logical_errors-batch
            p_l = ((i-1)*p_l + p_l_batch) / i
            r_a = ((i-1)*r_a + r_a_batch)/i

            if at_least_min_errors and num_logical_errors >= min_errors:
                break
            i += 1

        return p_l, r_a, num_logical_errors, i*batch

    def _simulate_batch(self, shots=1024):
        sampler = self.stim_circuit.compile_sampler()
        detection_events = sampler.sample(shots)

        # Filter events where the verification circuit flagged
        index_array = np.where(np.all(np.logical_not(detection_events[:, self.verification_measurements]), axis=1))[0]
        filtered_events = detection_events[index_array][:, self.x_measurements + self.z_measurements].astype(int)

        if len(filtered_events) == 0:  # All events were discarded
            return 0, shots
        
        state = filtered_events[:, self.data_measurements]

        if self.zero_state:
            checks = filtered_events[:, self.x_measurements]
            observables = self.code.Lz
            estimates = self.lut.batch_decode_x(checks)
        else:
            checks = filtered_events[:, self.z_measurements]
            observables = self.code.Lx
            estimates = self.lut.batch_decode_z(checks)

        corrected = state + estimates

        num_discarded = detection_events.shape[0]-filtered_events.shape[0]
        num_logical_errors = np.sum(np.any(corrected @ observables.T % 2!=0, axis=1))  # number of non-commuting corrected states
        return num_logical_errors, num_discarded

    
class LUTDecoder:
    """Lookup table decoder for a CSSState"""

    def __init__(self, code: CSSCode, init_LUTs: bool = True):
        self.code = code
        self.x_LUT = None
        self.z_LUT = None
        if init_LUTs:
            self.generate_x_LUT()
            self.generate_z_LUT()

    def batch_decode_x(self, syndromes: np.array) -> np.array:
        """Decode the X errors given a batch of syndromes."""
        return np.apply_along_axis(self.decode_x, 1, syndromes)

    def batch_decode_z(self, syndromes: np.array) -> np.array:
        """Decode the Z errors given a batch of syndromes."""
        return np.apply_along_axis(self.decode_z, 1, syndromes)
    
    def decode_x(self, syndrome: np.array) -> np.array:
        """Decode the X errors given a syndrome."""
        if self.x_LUT is None:
            self.generate_x_LUT()
        return self.x_LUT[syndrome.tobytes()]

    def decode_z(self, syndrome: np.array) -> np.array:
        """Decode the Z errors given a syndrome."""
        if self.z_LUT is None:
            self.generate_z_LUT()
        return self.z_LUT[syndrome.tobytes()]
        
    def generate_x_LUT(self) -> None:
        """Generate the lookup table for the X errors."""
        if self.x_LUT is not None:
            return

        self.x_LUT = self._generate_LUT(self.code.x_checks)
        if self.code.is_self_dual:
            self.z_LUT = self.x_LUT

    def generate_z_LUT(self) -> None:
        """Generate the lookup table for the Z errors."""
        if self.z_LUT is not None:
            return
        self.z_LUT = self._generate_LUT(self.code.z_checks)
        if self.code.is_self_dual:
            self.z_LUT = self.x_LUT
        
    def _generate_LUT(self, checks: np.array) -> dict:
        """Generate a lookup table for the stabilizer state.

        The lookup table maps error syndromes to their errors.
        """
        n_qubits = checks.shape[1]

        lut = defaultdict(list)
        for i in range(0, 2**n_qubits):
            state = np.array(list(np.binary_repr(i).zfill(n_qubits))).astype(np.int8)
            syndrome = checks @ state % 2
            lut[syndrome.tobytes()].append(state)

        # Sort according to weight
        for key, v in lut.items():
            lut[key] = np.array(sorted(v, key=lambda x: np.sum(x))[0])

        return lut


def _support(v: npt.NDArray[np.int_]) -> npt.NDArray[np.int_]:
    """Return the support of a vector."""
    return np.where(v != 0)[0]
