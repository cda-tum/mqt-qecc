# MQT QECC - A tool for Quantum Error Correcting Codes written in C++

A tool for quantum error correcting codes and numerical simulations developed by the
[Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/)
based on methods proposed in [[1]](todo)

QECC is part of the Munich Quantum Toolkit (MQT; formerly known as JKQ and developed by the
[Institute for Integrated Circuits](https://iic.jku.at/eda/) at the
[Johannes Kepler University Linz](https://jku.at)).

The tool can be used to decode quantum LDPC codes and conduct respective numerical evaluations.

At the moment we have implemented the general Union-Find
decoder [[2]](https://ieeexplore.ieee.org/abstract/document/9682738)
and an improved version of the algorithm as proposed in [[1]](todo). To,
to construct codes we use the open-source software by Joshka Roffe et
al: [[3]](https://github.com/quantumgizmos/bias_tailored_qldpc).

For more information, please visit [cda.cit.tum.de/](https://www.cda.cit.tum.de/).

If you have any questions, feel free to contact us via [quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de) or by
creating an issue on [GitHub](https://github.com/cda-tum/qecc/issues).

## Usage

MQT QECC is developed as a C++ library with an easy to use Python interface.

- In order to make the library as easy to use as possible (without compilation), we provide pre-built wheels for most
  common platforms (64-bit Linux, MacOS, Windows). These can be installed using
    ```bash
    pip install mqt.qecc
    ```
  However, in order to get the best performance out of QECC, it is recommended to build it locally from the source
  distribution (see [system requirements](#system-requirements)) via
    ```bash
    pip install  mqt.qecc --no-binary mqt.qecc
    ```
  This enables platform specific compiler optimizations that cannot be enabled on portable wheels.
- Once installed, start using it in Python:
  ```python
  from mqt.qecc import *
  results = decode(code)
  decoding_simulator.simulate_wer(...)
  ```

### System Requirements

Building (and running) is continuously tested under Linux, MacOS, and Windows using the
[latest available system versions for GitHub Actions](https://github.com/actions/virtual-environments).
However, the implementation should be compatible with any current C++ compiler supporting C++17 and a minimum CMake
version of 3.14.

`flint` [flint2](https://github.com/wbhart/flint2) is used for matrix operations in GF(2). For an installation
guide please refer to the official [flint documentation](https://flintlib.org/doc/building.html). The easiest way
is to download and build locally and then use the environment variable `LD_LIBRARY_PATH`
to specify the location of flint (e.g. `export LD_LIBRARY_PATH=/usr/local/lib` on UNIX).

### Configure, Build, and Install

To start off, clone this repository using

```shell
git clone --recurse-submodules -j8 https://github.com/cda-tum/qecc 
```

Note the `--recurse-submodules` flag. It is required to also clone all the required submodules.
If you happen to forget passing the flag on your initial clone, you can initialize all the submodules by
executing `git submodule update --init --recursive` in the main project directory.

Our projects use CMake as the main build configuration tool. Building a project using CMake is a two-stage process.
First, CMake needs to be *configured* by calling

```shell 
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

This tells CMake to search the current directory `.` (passed via `-S`) for a *CMakeLists.txt* file and process it into a
directory `build` (passed via `-B`).
The flag `-DCMAKE_BUILD_TYPE=Release` tells CMake to configure a *Release* build (as opposed to, e.g., a *Debug* build).

After configuring with CMake, the project can be built by calling

```shell
cmake --build build --config Release
```

This tries to build the project in the `build` directory (passed via `--build`).
Some operating systems and developer environments explicitly require a configuration to be set, which is why
the `--config` flag is also passed to the build command. The flag `--parallel <NUMBER_OF_THREADS>` may be added to
trigger a parallel build.

Building the project this way generates

- the main library `libqec.a` (Unix) / `libqec_lib.lib` (Windows) in the `build/src` directory
- a test executable `qec_test` containing a small set of unit tests in the `build/test` directory

### Extending the Python Bindings

To extend the Python bindings you can locally install the package in edit mode, so that changes in the Python code are
instantly available.
The following example assumes you have a [virtual environment](https://docs.python.org/3/library/venv.html) set up and
activated.

```commandline
(venv) $ pip install cmake
(venv) $ pip install --editable .
```

If you change parts of the C++ code, you have to run the second line to make the changes visible in Python.

## Reference

If you use our tool for your research, we will be thankful if you refer to it by citing the appropriate publications.

```bibtex
@article{DBLP:conf/to/do,
  author    = {todo},
  title     = {todo},
  journal   = {todo},
  year      = {2022}
}
```
