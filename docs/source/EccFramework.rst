ECC Framework: Automatic Implementation and Evaluation of Error-Correcting Codes for Quantum Computing
======================================================================================================

The QFR library offers means for automatic implementation and evaluation of error-correcting codes for quantum computing. More precisely, the library allows to automatically apply different error correction schemes to quantum circuits provided as openqasm files. The "protected" quantum circuits can then be exported again in the form of openqasm files or can be directly used for noise-aware quantum circuit simulation. For the latter case, we also provide a wrapper script which makes use of the provided framework to apply error correction schemes to circuits and directly simulate those circuits using qiskit.

**Note: The ECC framework is only available within the current branch and can only be installed directly from source**

Installation
############

If you have not done so already, clone the repository using:

.. code-block:: console

  git clone --recurse-submodules -j8 https://github.com/pichristoph/qfr.git

Make sure you are in the main project directory for the next steps. Switch to the branch feature/ecc,

.. code-block:: console

  cd qfr
  git switch feature/ecc

and (if necessary), update the submodules.

.. code-block:: console

  git submodule update --init --recursive

Then, the ECC framework can be installed using pip

.. code-block:: console

  (venv) pip install --editable .

If you want to use Qiskit for quantum circuit simulation, you need to install it as well

.. code-block:: console

  (venv) pip install qiskit

Usage
#####

Having the Python module installed, error correcting codes can be applied using apply_ecc of module qfr, like so

.. code-block:: file

  from mqt import qfr

  file = "path/to/qasm/file.qasm" # Path to the openqasm file the quantum circuit shall be loaded from
  ecc = "Q7Steane" # Error correction code that shall be applied to the quantum circuit
  ecc_frequency = 100 # After how many times a qubit is used, error correction is applied
  ecc_mc = False # Only allow single controlled gates in the created quantum circuit
  ecc_cf = False # Only allow single clifford gates in the created quantum circuit

  result = qfr.apply_ecc(file, ecc, ecc_frequency, ecc_mc, ecc_cf)

  # print the resulting circuit
  print(result["circ"])

Currently, the error correction schemes Q3Shor, Q5Laflamme, Q7Steane, Q9Shor, Q9Surface, and Q18Surface are supported.

We provide a wrapper script for applying error correction to quantum circuits (provided as openQasm) and followed by a noise-aware quantum circuit simulation (using qiskit). The script can be used like this:

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
    -f F        Path to a openqasm file
    -e E        Export circuit, with error correcting code applied, as openqasm circuit instead of simulation it (e.g., -e "/path/to/new/openqasm_file") (Default=None)
    -fs FS      Specify a simulator (Default: "statevector_simulator" for simulation without noise, "aer_simulator_density_matrix", for deterministic noise-aware simulation"aer_simulator_statevector", for stochastic noise-
                aware simulation). Available: [AerSimulator('aer_simulator'), AerSimulator('aer_simulator_statevector'), AerSimulator('aer_simulator_density_matrix'), AerSimulator('aer_simulator_stabilizer'),
                AerSimulator('aer_simulator_matrix_product_state'), AerSimulator('aer_simulator_extended_stabilizer'), AerSimulator('aer_simulator_unitary'), AerSimulator('aer_simulator_superop'),
                QasmSimulator('qasm_simulator'), StatevectorSimulator('statevector_simulator'), UnitarySimulator('unitary_simulator'), PulseSimulator('pulse_simulator')]
    -ecc ECC    Specify a ecc to be applied to the circuit. Currently available are Q3Shor, Q5Laflamme, Q7Steane, Q9Shor, Q9Surface, and Q18Surface (Default=none)
    -fq FQ      Specify after how many qubit usages error correction is applied to it (Default=100)
    -mc MC      Only allow single controlled gates (Default=False)
    -cf CF      Only allow clifford operations (Default=False)

Available error-correcting codes and operations
###############################################

.. list-table:: Available error-correcting codes and operations
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

Properties of the implemented error-correcting codes
####################################################

.. list-table:: Available error-correcting codes and operations
  :widths: 22 13 13 13 13 13 13
  :header-rows: 1

  * - Feature
    - Q3Shor
    - Q5Laflamme
    - Q7Steane
    - Q9Shor
    - Q9Surface
    - Q18Surface
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
    - ✔️
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

More-detailed information about the error-correcting codes can be found in the README information TODO.
