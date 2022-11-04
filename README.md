[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![CI](https://img.shields.io/github/workflow/status/cda-tum/qecc/C++?style=flat-square&logo=github&label=c%2B%2B)](https://github.com/cda-tum/qecc/actions/workflows/ci.yml)
[![Python CI](https://img.shields.io/github/workflow/status/cda-tum/qecc/Python?style=flat-square&logo=github&label=python)](https://github.com/cda-tum/qecc/actions/workflows/python-ci.yml)
[![Bindings](https://img.shields.io/github/workflow/status/cda-tum/qecc/Python%20Packaging?style=flat-square&logo=github&label=packaging)](https://github.com/cda-tum/qecc/actions/workflows/deploy.yml)
[![codecov](https://img.shields.io/codecov/c/github/cda-tum/qecc?style=flat-square&logo=codecov)](https://codecov.io/gh/cda-tum/qecc)

# QECC: An MQT tool for Quantum Error Correcting Codes written in C++

:warning: **This project is still in early development and breaking changes might happen frequently.**

(Additionally to the basic numerical results already provided, further data will be published continually)

A tool for quantum error correcting codes and numerical simulations developed by the
[Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/)
based on methods proposed in [[1]](https://arxiv.org/abs/2209.01180). QECC is part of the Munich Quantum Toolkit (MQT).

The tool can be used to decode quantum LDPC codes and conduct respective numerical simulations.

At the moment the general QLDPC
decoder [[2]](https://ieeexplore.ieee.org/abstract/document/9682738)
and a heuristic (which improves the runtime of the algorithm) [[1]](https://arxiv.org/abs/2209.01180) are implemented.
Currently, open-source software by Joschka Roffe et
al.: [[3]](https://github.com/quantumgizmos/bias_tailored_qldpc) is used to construct codes (toric, lifted product and
hypergraph product).

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

**Detailed documentation on all available methods, options, and input formats is available
at [ReadTheDocs](https://qecc.readthedocs.io/en/latest/).**

## System Requirements and Building

The implementation is compatible with any C++17 compiler and a minimum CMake version of 3.14.
Please refer to the [documentation](https://qecc.readthedocs.io/en/latest/) on how to build the project.

Building (and running) is continuously tested under Linux, macOS, and Windows using the
[latest available system versions for GitHub Actions](https://github.com/actions/virtual-environments).

## Reference

If you use our tool for your research, we will be thankful if you refer to it by citing the appropriate publication:

[![a](https://img.shields.io/static/v1?label=arXiv&message=2011.07288&color=inactive&style=flat-square)](https://arxiv.org/abs/2209.01180)
L. Berent, L. Burgholzer, and R.
Wille, "[Software Tools for Decoding Quantum Low-Density Parity Check Codes](https://arxiv.org/abs/2209.01180),"
in Asia and South Pacific Design Automation Conference (ASP-DAC), 2023
