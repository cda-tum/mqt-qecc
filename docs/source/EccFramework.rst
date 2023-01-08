ECC Framework
======================================================================================================

The QECC library offers means for automatic implementation and evaluation of error-correcting codes for quantum
computing. More precisely, the library allows to automatically apply different error correction schemes to quantum
circuits provided as OpenQASM files or Qiskit QuantumCircuit objects. The "protected" quantum circuits can then be exported again in the form of OpenQASM files or can be directly used for noise-aware quantum circuit simulation. For the latter case, we also provide a wrapper script which makes use of the provided framework to apply error correction schemes to circuits and directly simulate those circuits using qiskit.

Usage
#####

Having the Python module installed, error correcting codes can be applied using apply_ecc of module qecc, like so

.. code-block:: python

  from mqt import qecc

  file = "path/to/qasm/file.qasm"  # Path to the OpenQASM file the quantum circuit shall be loaded from
  ecc = "Q7Steane"  # Error correction code that shall be applied to the quantum circuit
  ecc_frequency = 100  # After how many times a qubit is used, error correction is applied

  result = qecc.apply_ecc(file, ecc, ecc_frequency)

  print(result["circ"])

Currently, the error correction schemes Q3Shor, Q5Laflamme, Q7Steane, Q9Shor, Q9Surface, and Q18Surface are supported.

We provide a wrapper script for applying error correction to quantum circuits (provided as OpenQASM) and followed by a noise-aware quantum circuit simulation (using qiskit). The script can be used like this:

.. code-block:: file

  $ /venv/ecc_qiskit_wrapper -ecc Q7Steane -fq 100 -m D -p 0.0001 -n 2000 -fs aer_simulator_stabilizer -s 0 -f  ent_simple1000_n2.qasm
  _____Trying to simulate with D(prob=0.0001, shots=2000, n_qubits=17) Error______
  State |00> probability 0.515
  State |01> probability 0.0055
  State |10> probability 0.0025
  State |11> probability 0.477

The script offers a help function, which displays available parameters:

.. code-block:: console

  $ /venv/ecc_qiskit_wrapper --help
  usage: ecc_qiskit_wrapper [-h] [-m M] [-p P] [-n N] [-s S] -f F [-e E] [-fs FS] [-ecc ECC] [-fq FQ] [-mc MC] [-cf CF]

  QiskitWrapper interface with error correction support!

  optional arguments:
    -h, --help  show this help message and exit
    -m M        Define the error_channels (e.g., -m APD), available errors channels are amplitude damping (A), phase flip (P), bit flip (B), and depolarization (D) (Default="D")
    -p P        Set the noise probability (Default=0.001)
    -n N        Set the number of shots. 0 for deterministic simulation (Default=2000)
    -s S        Set a seed (Default=0)
    -f F        Path to a OpenQASM file
    -e E        Export circuit, with error correcting code applied, as OpenQASM circuit instead of simulation it (e.g., -e "/path/to/new/openqasm_file") (Default=None)
    -fs FS      Specify a simulator (Default: "statevector_simulator" for simulation without noise, "aer_simulator_density_matrix", for deterministic noise-aware simulation"aer_simulator_statevector", for stochastic noise-
                aware simulation). Available: [AerSimulator('aer_simulator'), AerSimulator('aer_simulator_statevector'), AerSimulator('aer_simulator_density_matrix'), AerSimulator('aer_simulator_stabilizer'),
                AerSimulator('aer_simulator_matrix_product_state'), AerSimulator('aer_simulator_extended_stabilizer'), AerSimulator('aer_simulator_unitary'), AerSimulator('aer_simulator_superop'),
                QasmSimulator('qasm_simulator'), StatevectorSimulator('statevector_simulator'), UnitarySimulator('unitary_simulator'), PulseSimulator('pulse_simulator')]
    -ecc ECC    Specify a ecc to be applied to the circuit. Currently available are Q3Shor, Q5Laflamme, Q7Steane, Q9Shor, Q9Surface, and Q18Surface (Default=none)
    -fq FQ      Specify after how many qubit usages error correction is applied to it (Default=100)

Available error-correcting codes
#####################################################

Properties
----------

.. list-table:: Properties of available error-correcting codes
  :widths: 22 13 13 13 13 13 13
  :header-rows: 1

  * - Feature
    - Q3Shor :cite:labelpar:`ShorCodes`
    - Q5Laflamme :cite:labelpar:`LaflammeCodes`
    - Q7Steane :cite:labelpar:`SteaneCode`
    - Q9Shor :cite:labelpar:`ShorCodes`
    - Q9Surface :cite:labelpar:`WoottonMinimalSurfaceCode`
    - Q18Surface (based on :cite:labelpar:`FowlerSurfaceCodes`)
  * - able to detect bit flips
    - ✔️
    - ✔️
    - ✔️
    - ✔️
    - ✔️
    - ✔️
  * - able to detect phase flips
    - ✖️
    - ✔️
    - ✔️
    - ✔️
    - ✔️
    - ✖️*
  * - #qubits for n logical qubits
    - 3n+2
    - 5n+4
    - 7n+3
    - 9n+8
    - 9n+8
    - 36n
  * - #classical bits (total)
    - 2
    - 5
    - 3
    - 8
    - 8
    - 16

\* Planned to work, but not fully implemented yet

Available logical operations
----------------------------

.. list-table:: Available operations for each error-correcting code
  :widths: 22 13 13 13 13 13 13
  :header-rows: 1

  * - Operation
    - Q3Shor
    - Q5Laflamme
    - Q7Steane
    - Q9Shor
    - Q9Surface
    - Q18Surface
  * - Pauli (X, Y, Z)
    - ✔️
    - ✔️
    - ✔️
    - ✔️
    - ✔️
    - ✔️
  * - controlled Pauli (CX,CY,CZ)
    - ✔️
    - ✖️
    - ✔️
    - ✔️
    - ✔️
    - ✖️
  * - Hadamard
    - ⚠️
    - ✖️
    - ✔️
    - ✖️
    - ✔️
    - ✔️
  * - S, S†, T, T†
    - ✔️
    - ✖️
    - ✔️
    - ✖️
    - ✖️
    - ✖️

⚠️ = operation is applied without the scheme of the error-correcting code (i.e. decoding and encoding is performed before/afterwards, respectively, and the operation is encoded as-is)
