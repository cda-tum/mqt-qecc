Installation
============

QECC is mainly developed as a C++ library.
In order to make the tool as accessible as possible, it comes with an easy-to-use Python interface.

We encourage installing QECC via pip (preferably in a `virtual environment <https://docs.python.org/3/library/venv.html>`_):

    .. code-block:: console

        (venv) $ pip install mqt.qecc

In most practical cases (under 64-bit Linux, MacOS incl. Apple Silicon, and Windows), this requires no compilation and merely downloads and installs a platform-specific pre-built wheel.

.. note::
    In order to set up a virtual environment, you can use the following commands:

    .. code-block:: console

        $ python3 -m venv venv
        $ source venv/bin/activate

    If you are using Windows, you can use the following commands instead:

    .. code-block:: console

        $ python3 -m venv venv
        $ venv\Scripts\activate.bat

    It is recommended to make sure that you are using the latest version of pip, setuptools, and wheel before trying to install the project:

    .. code-block:: console

        (venv) $ pip install --upgrade pip setuptools wheel

A Detailed Walk Through
#######################

First, save the following lines as :code:`steane_example.py` in a folder where you want to install QECC and run the example:

    .. code-block:: python

        from mqt.qecc import *
        import numpy as np

        H = [[1, 0, 0, 1, 0, 1, 1], [0, 1, 0, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1, 1]]
        code = Code(H, H)
        decoder = UFHeuristic()
        decoder.set_code(code)
        x_err = sample_iid_pauli_err(code.n, 0.05)
        decoder.decode(code.get_x_syndrome(x_err))
        result = decoder.result
        residual_err = np.array(x_err) ^ np.array(result.estimate)

        print(result)
        print(code.is_x_stabilizer(residual_err))
        print(np.array(x_err).astype(int))

Then, the following snippet shows the installation process from setting up the virtual environment to running a small example program.

    .. code-block:: console

        $ python3 -m venv venv
        $ . venv/bin/activate
        (venv) $ pip install -U pip setuptools wheel
        (venv) $ pip install mqt.qecc
        (venv) $ python3 steane_example.py
        {
            "decodingTime(ms)": 0,
            "estimate": "[0,1,0,0,0,0,0]"
        }
        True
        [0 1 0 0 0 0 0]


Building from Source for Performance
####################################

In order to get the best performance out of QECC and enable platform-specific compiler optimizations that cannot be enabled on portable wheels, it is recommended to build the package from source via:

    .. code-block:: console

        (venv) $ pip install mqt.qecc --no-binary mqt.qecc

This requires a `C++ compiler <https://en.wikipedia.org/wiki/List_of_compilers#C++_compilers>`_ compiler supporting *C++17* and a minimum `CMake <https://cmake.org/>`_ version of *3.19*.

The library is continuously tested under Linux, MacOS, and Windows using the `latest available system versions for GitHub Actions <https://github.com/actions/virtual-environments>`_.
In order to access the latest build logs, visit `qecc/actions/workflows/ci.yml <https://github.com/cda-tum/mqt-qecc/actions/workflows/ci.yml>`_.
