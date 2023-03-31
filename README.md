[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![CI](https://img.shields.io/github/actions/workflow/status/cda-tum/qecc/ci.yml?branch=main&style=flat-square&logo=github&label=c%2B%2B)](https://github.com/cda-tum/qecc/actions/workflows/ci.yml)
[![Python CI](https://img.shields.io/github/actions/workflow/status/cda-tum/qecc/python-ci.yml?branch=main&style=flat-square&logo=github&label=python)](https://github.com/cda-tum/qecc/actions/workflows/python-ci.yml)
[![Bindings](https://img.shields.io/github/actions/workflow/status/cda-tum/qecc/deploy.yml?branch=main&style=flat-square&logo=github&label=packaging)](https://github.com/cda-tum/qecc/actions/workflows/deploy.yml)
[![codecov](https://img.shields.io/codecov/c/github/cda-tum/qecc?style=flat-square&logo=codecov)](https://codecov.io/gh/cda-tum/qecc)

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/cda-tum/qecc/main/docs/source/_static/mqt_light.png" width="60%">
    <img src="https://raw.githubusercontent.com/cda-tum/qecc/main/docs/source/_static/mqt_dark.png" width="60%">
  </picture>
  </p>

# QECC: An MQT tool for Quantum Error Correcting Codes written in C++

:warning: **This project is still in early development and breaking changes might happen frequently.**

(Additionally to the basic numerical results already provided, further data will be published continually)

A tool for quantum error correcting codes and numerical simulations developed by the
[Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/)
based on methods proposed in [[1]](https://arxiv.org/abs/2209.01180). QECC is part of the Munich Quantum Toolkit (MQT).

The tool can be used to:

- Decode quantum LDPC codes and conduct respective numerical simulations.
  - At the moment the general QLDPC
    decoder [[2]](https://ieeexplore.ieee.org/abstract/document/9682738)
    and a heuristic (which improves the runtime of the algorithm) [[1]](https://arxiv.org/abs/2209.01180) are
    implemented.
    Currently, open-source software by Joschka Roffe et
    al.: [[3]](https://github.com/quantumgizmos/bias_tailored_qldpc) is used to construct codes (toric, lifted product
    and
    hypergraph product).
- Decode (triangular) color codes and conduct respective numerical simulations.
  - The decoder is based on an analogy to the classical LightsOut puzzle and formulated as a MaxSAT problem. The SMT solver
    Z3 is used to determine minimal solutions of the MaxSAT problem, resulting in minimum-weight decoding estimates.
- Apply error correction to quantum circuits.
  - The framework allows to apply different QECC schemes to quantum circuits and either exports the resulting
    circuits or simulates them using Qiskit [[4]](https://qiskit.org/). Currently, six different ECCs are supported
    with varying extent of functionality.

<p align="center">
  <a href="https://qecc.readthedocs.io/en/latest/">
  <img width=30% src="https://img.shields.io/badge/documentation-blue?style=for-the-badge&logo=read%20the%20docs" alt="Documentation" />
  </a>
</p>

If you have any questions, feel free to contact us via [quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de) or by
creating an issue on [GitHub](https://github.com/cda-tum/qecc/issues).

## Getting Started

QECC is available via [PyPI](https://pypi.org/project/mqt.qecc/) for Linux, macOS, and Windows.

```console
(venv) $ pip install mqt.qecc
```

The following code gives an example on the usage:

### Example for decoding quantum LDPC codes

```python3
from mqt.qecc import *
import numpy as np

H = [
    [1, 0, 0, 1, 0, 1, 1],
    [0, 1, 0, 1, 1, 0, 1],
    [0, 0, 1, 0, 1, 1, 1]
]
code = Code(H, H)
decoder = UFHeuristic()
decoder.set_code(code)
x_err = sample_iid_pauli_err(code.N, 0.05)
decoder.decode(code.get_x_syndrome(x_err))
result = decoder.result
print(result)
residual_err = np.array(x_err) ^ np.array(result.estimate)
print(code.is_x_stabilizer(residual_err))
```

### Example for decoding color codes

Simply running the following code will perform a numerical analysis of the MaxSAT color code decoder for an instance of
the distance-21 triangular color code with a bit-flip error rate of 0.01 and 1000 simulations.

```python3
from mqt.qecc.cc_decoder import decoder

d = 21  # distance of the triangular code to simulate
p = 0.01  # (bit-flip) error rate
nr_sims = 1000  # number of simulations to run
decoder.run(distance=d, error_rate=p, nr_sims=nr_sims)
```

The dataset used in the paper evaluation on decoding quantum color codes is available on Zenodo:
[![a](https://img.shields.io/static/v1?label=DOI&message=10.5281/zenodo.7760135&color=inactive&style=flat-square)](https://doi.org/10.5281/zenodo.7760135)

### Example for applying error correction to a circuit

```python3
from mqt import qecc

file = "path/to/qasm/file.qasm"  # Path to the OpenQASM file the quantum circuit shall be loaded from
ecc = "Q7Steane"  # Error correction code that shall be applied to the quantum circuit
ecc_frequency = 100  # After how many times a qubit is used, error correction is applied

result = qecc.apply_ecc(file, ecc, ecc_frequency)

# print the resulting circuit as OpenQASM string
print(result["circ"])
```

A wrapper script for applying error correction to quantum circuits (provided as OpenQASM) and performing a
noise-aware quantum circuit simulation (using Qiskit) is provided. The script can be used like this:

```bash
$ (venv) ecc_qiskit_wrapper -ecc Q7Steane -fq 100 -m D -p 0.0001 -n 2000 -fs aer_simulator_stabilizer -s 0 -f  ent_simple1000_n2.qasm
_____Trying to simulate with D (prob=0.0001, shots=2000, n_qubits=17, error correction=Q7Steane) Error______
State |00> probability 0.515
State |01> probability 0.0055
State |10> probability 0.0025
State |11> probability 0.477
```

**Detailed documentation on all available methods, options, and input formats is available
at [ReadTheDocs](https://qecc.readthedocs.io/en/latest/).**

## System Requirements and Building

The implementation is compatible with any C++17 compiler and a minimum CMake version of 3.19.
Please refer to the [documentation](https://qecc.readthedocs.io/en/latest/) on how to build the project.

Building (and running) is continuously tested under Linux and macOS using the
[latest available system versions for GitHub Actions](https://github.com/actions/virtual-environments).
Windows support is currently experimental.

## Reference

If you use our tool for your research, we will be thankful if you refer to it by citing the appropriate publication:

- [![a](https://img.shields.io/static/v1?label=arXiv&message=2303.14237&color=inactive&style=flat-square)](https://arxiv.org/abs/2303.14237)
  L. Berent, L. Burgholzer, P.J. Derks, J. Eisert, and R. Wille, "Decoding quantum color codes with MaxSAT".

- [![a](https://img.shields.io/static/v1?label=arXiv&message=2301.05731&color=inactive&style=flat-square)](https://arxiv.org/pdf/2301.05731)
  T. Grurl, C. Pichler, J. Fuss and R. Wille, "Automatic Implementation and Evaluation of Error-Correcting Codes for
  Quantum Computing: An Open-Source Framework for Quantum Error-Correction," in International Conference on VLSI
  Design and International Conference on Embedded Systems (VLSID), 2023

- [![a](https://img.shields.io/static/v1?label=arXiv&message=2011.07288&color=inactive&style=flat-square)](https://arxiv.org/abs/2209.01180)
  L. Berent, L. Burgholzer, and R.
  Wille, "[Software Tools for Decoding Quantum Low-Density Parity Check Codes](https://arxiv.org/abs/2209.01180),"
  in Asia and South Pacific Design Automation Conference (ASP-DAC), 2023
