"""Simulation of Non-deterministic fault tolerant state preparation."""

from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
import stim
from qiskit.converters import circuit_to_dag, dag_to_circuit

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt
    from qiskit import QuantumCircuit

    from ..code import CSSCode


class NoisyNDFTStatePrepSimulator:
    """A noisy state preparation circuit."""

    def __init__(self, state_prep_circ: QuantumCircuit, code: CSSCode, p: float, zero_state: bool = True) -> None:
        """Initialize the simulator.

        Args:
            state_prep_circ: The state preparation circuit.
            code: The code to simulate.
            p: The error rate.
        zero_state: Whether thezero state is prepared or nor.
        """
        self.circ = state_prep_circ
        self.num_qubits = state_prep_circ.num_qubits
        self.code = code
        self.p = p
        self.zero_state = zero_state
        # Store which measurements are X, Z or data measurements.
        # The indices refer to the indices of the measurements in the stim circuit.
        self.verification_measurements = []  # type: list[int]
        self.x_measurements = []  # type: list[int]
        self.z_measurements = []  # type: list[int]
        self.data_measurements = []  # type: list[int]
        self.n_measurements = 0
        self.stim_circ = stim.Circuit()
        self.decoder = LUTDecoder(code)
        self.set_p(p)

    def set_p(self, p: float) -> None:
        """Set the error rate.

        This reinitializes the stim circuit.

        Args:
        p: The error rate.
        """
        self.verification_measurements = []
        self.x_measurements = []
        self.z_measurements = []
        self.data_measurements = []
        self.n_measurements = 0
        self.p = p
        self.stim_circ = self.to_stim_circ()
        self.num_qubits = self.stim_circ.num_qubits - len(self.verification_measurements)
        self.measure_stabilizers()
        if self.zero_state:
            self.measure_z()
        else:
            self.measure_x()

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

        def idle_error(used_qubits: list[int]) -> None:
            for q in self.circ.qubits:
                qubit = self.circ.find_bit(q)[0]
                if initialized[qubit] and qubit not in used_qubits:
                    stim_circuit.append_operation("DEPOLARIZE1", [self.circ.find_bit(q)[0]], [2 * self.p / 3])

        dag = circuit_to_dag(self.circ)
        layers = dag.layers()
        used_qubits = []  # type: list[int]

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
                            stim_circuit.append_operation("Z_ERROR", [ctrl], [2 * self.p / 3])  # Wrong initialization
                        else:
                            stim_circuit.append_operation("X_ERROR", [ctrl], [2 * self.p / 3])  # Wrong initialization
                        initialized[ctrl] = True
                    if not initialized[target]:
                        stim_circuit.append_operation("X_ERROR", [target], [2 * self.p / 3])  # Wrong initialization
                        initialized[target] = True

                    stim_circuit.append_operation("CX", [ctrl, target])
                    stim_circuit.append_operation("DEPOLARIZE2", [ctrl, target], [self.p])
                    used_qubits.extend([ctrl, target])

                elif gate[0].name == "measure":
                    stim_circuit.append_operation(
                        "X_ERROR", [self.circ.find_bit(q)[0] for q in gate[1]], [2 * self.p / 3]
                    )
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
            anc = self.stim_circ.num_qubits
            self.stim_circ.append_operation("H", [anc])
            for q in supp:
                self.stim_circ.append_operation("CX", [anc, q])
            self.stim_circ.append_operation("MRX", [anc])
            self.x_measurements.append(self.n_measurements)
            self.n_measurements += 1

        for check in self.code.Hz:
            supp = _support(check)
            anc = self.stim_circ.num_qubits
            for q in supp:
                self.stim_circ.append_operation("CX", [q, anc])
            self.stim_circ.append_operation("MRZ", [anc])
            self.z_measurements.append(self.n_measurements)
            self.n_measurements += 1

    def measure_z(self) -> None:
        """Measure all data qubits in the Z basis."""
        self.data_measurements = [self.n_measurements + i for i in range(self.num_qubits)]
        self.n_measurements += self.num_qubits
        self.stim_circ.append_operation("MRZ", list(range(self.num_qubits)))

    def measure_x(self) -> None:
        """Measure all data qubits in the X basis."""
        self.data_measurements = [self.n_measurements + i for i in range(self.num_qubits)]
        self.n_measurements += self.num_qubits
        self.stim_circ.append_operation("MRX", list(range(self.num_qubits)))

    def logical_error_rate(
        self,
        shots: int = 100000,
        shots_per_batch: int = 100000,
        at_least_min_errors: bool = True,
        min_errors: int = 500,
    ) -> tuple[float, float, int, int]:
        """Estimate the logical error rate of the code.

        Args:
            shots: The number of shots to use.
            shots_per_batch: The number of shots per batch.
            at_least_min_errors: Whether to continue simulating until at least min_errors are found.
            min_errors: The minimum number of errors to find before stopping.
        """
        batch = min(shots_per_batch, shots)
        p_l = 0.0
        r_a = 0.0

        num_logical_errors = 0

        if self.zero_state:
            self.decoder.generate_x_LUT()
        else:
            self.decoder.generate_z_LUT()

        i = 1
        while i <= int(np.ceil(shots / batch)) or at_least_min_errors:
            num_logical_errors_batch, discarded_batch = self._simulate_batch(batch)

            if discarded_batch != batch:
                p_l_batch = num_logical_errors_batch / (batch - discarded_batch)
                p_l = ((i - 1) * p_l + p_l_batch) / i
                
            r_a_batch = 1 - discarded_batch / batch

            # Update statistics
            num_logical_errors += num_logical_errors_batch
            r_a = ((i - 1) * r_a + r_a_batch) / i

            if at_least_min_errors and num_logical_errors >= min_errors:
                break
            i += 1

        return p_l, r_a, num_logical_errors, i * batch

    def _simulate_batch(self, shots: int = 1024) -> tuple[int, int]:
        sampler = self.stim_circ.compile_sampler()
        detection_events = sampler.sample(shots)

        # Filter events where the verification circuit flagged
        index_array = np.where(np.all(np.logical_not(detection_events[:, self.verification_measurements]), axis=1))[0]
        filtered_events = detection_events[index_array].astype(np.int8)

        if len(filtered_events) == 0:  # All events were discarded
            return 0, shots

        state = filtered_events[:, self.data_measurements]

        if self.zero_state:
            checks = filtered_events[:, self.z_measurements]
            observables = self.code.Lz
            estimates = self.decoder.batch_decode_x(checks)
        else:
            checks = filtered_events[:, self.x_measurements]
            observables = self.code.Lx
            estimates = self.decoder.batch_decode_z(checks)

        corrected = state + estimates

        num_discarded = detection_events.shape[0] - filtered_events.shape[0]
        num_logical_errors = np.sum(
            np.any(corrected @ observables.T % 2 != 0, axis=1)
        )  # number of non-commuting corrected states
        return num_logical_errors, num_discarded


class LUTDecoder:
    """Lookup table decoder for a CSSState."""

    def __init__(self, code: CSSCode, init_LUTs: bool = True) -> None:
        """Initialize the decoder.

        Args:
            code: The code to decode.
            init_LUTs: Whether to initialize the lookup tables at object creation.
        """
        self.code = code
        self.x_LUT = {}  # type: dict[bytes, npt.NDArray[np.int8]]
        self.z_LUT = {}  # type: dict[bytes, npt.NDArray[np.int8]]
        if init_LUTs:
            self.generate_x_LUT()
            self.generate_z_LUT()

    def batch_decode_x(self, syndromes: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
        """Decode the X errors given a batch of syndromes."""
        return np.apply_along_axis(self.decode_x, 1, syndromes)

    def batch_decode_z(self, syndromes: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
        """Decode the Z errors given a batch of syndromes."""
        return np.apply_along_axis(self.decode_z, 1, syndromes)

    def decode_x(self, syndrome: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
        """Decode the X errors given a syndrome."""
        if len(self.x_LUT) == 0:
            self.generate_x_LUT()
        return self.x_LUT[syndrome.tobytes()]

    def decode_z(self, syndrome: npt.NDArray[np.int8]) -> npt.NDArray[np.int8]:
        """Decode the Z errors given a syndrome."""
        if len(self.z_LUT) == 0:
            self.generate_z_LUT()
        return self.z_LUT[syndrome.tobytes()]

    def generate_x_LUT(self) -> None:
        """Generate the lookup table for the X errors."""
        if len(self.x_LUT) != 0:
            return

        self.x_LUT = LUTDecoder._generate_LUT(self.code.Hz)
        if self.code.is_self_dual():
            self.z_LUT = self.x_LUT

    def generate_z_LUT(self) -> None:
        """Generate the lookup table for the Z errors."""
        if len(self.z_LUT) != 0:
            return
        self.z_LUT = LUTDecoder._generate_LUT(self.code.Hx)
        if self.code.is_self_dual():
            self.z_LUT = self.x_LUT

    @staticmethod
    def _generate_LUT(checks: npt.NDArray[np.int_]) -> dict[bytes, npt.NDArray[np.int_]]:
        """Generate a lookup table for the stabilizer state.

        The lookup table maps error syndromes to their errors.
        """
        n_qubits = checks.shape[1]

        syndromes = defaultdict(list)
        lut = {}  # type: dict[bytes, npt.NDArray[np.int8]]
        for i in range(2**n_qubits):
            state = np.array(list(np.binary_repr(i).zfill(n_qubits))).astype(np.int8)  # type: npt.NDArray[np.int_]
            syndrome = checks @ state % 2
            syndromes[syndrome.astype(np.int8).tobytes()].append(state)

        # Sort according to weight
        for key, v in lut.items():
            lut[key] = np.array(min(v, key=np.sum))

        return lut


def _support(v: npt.NDArray[np.int_]) -> npt.NDArray[np.int_]:
    """Return the support of a vector."""
    return np.where(v != 0)[0]
