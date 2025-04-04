# Installation

QECC is mainly developed as a C++ library.
In order to make the tool as accessible as possible, it comes with an easy-to-use Python interface.

We encourage installing QECC via pip (preferably in a [virtual environment](https://docs.python.org/3/library/venv.html)):

> ```console
> (venv) $ pip install mqt.qecc
> ```

:::{note}
In order to set up a virtual environment, you can use the following commands:

```console
$ python3 -m venv venv
$ source venv/bin/activate
```

If you are using Windows, you can use the following commands instead:

```console
$ python3 -m venv venv
$ venv\Scripts\activate.bat
```

:::

% A Detailed Walk Through
% #######################
%
% First, save the following lines as :code:`steane_example.py` in a folder where you want to install QECC and run the example:
%
% .. code-block:: python
%
% from mqt.qecc import CSSCode
% from mqt.qecc.circuit_synthesis import \*
% import numpy as np
%
% H = np.array([[1, 0, 0, 1, 0, 1, 1], [0, 1, 0, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1, 1]], dtype=np.int8)
% code = CSSCode(H, H)
% non_ft_circ = heuristic_prep_circuit(code)
% ft_circ = gate_optimal_verification_circuit(non_ft_circ)
% sim = NoisyNDFTStatePrepSimulator(ft_circ, code)
% sim.set_p(0.01, 0.0001)
% p_l, r_a, n_err, n_total = logical_error_rate(shots=100000, at_least_min_error=False)
%
% print(p_l, r_a, n_err, n_total)
%
%
% Then, the following snippet shows the installation process from setting up the virtual environment to running a small example program.
%
% .. code-block:: console
%
% $ python3 -m venv venv
% $ . venv/bin/activate
% (venv) $ pip install -U pip setuptools wheel
% (venv) $ pip install mqt.qecc
% (venv) $ python3 steane_example.py
% {
% "decodingTime(ms)": 0,
% "estimate": "[0,1,0,0,0,0,0]"
% }
% True
% [0 1 0 0 0 0 0]
