# Analog Information Decoding

This submodule provides means to conduct decoding simulations for quantum CSS codes with
analog syndrome information as proposed in the corresponding paper [1].
Proper integration and setup is a work in progress.
The main functionality is provided in the `simulators` module, which contains the main classes
that can be used to conduct several types of simulations:

- `ATD_Simulator`: Analog Tanner graph decoding,
- `Single-Shot Simulator`: Analog Single-Shot decoding with meta checks.
- `QSS_Simulator`: Quasi-Single-Shot decoding, and

Moreover, `memory_experiment` contains methods for analog overlapping window decoding, in
particular the `decode_multiround` method.

## Results

The `results` directory contains the results used in the paper [1].

## Codes

The `codes` directory contains the parity-check matrices of the codes used in the paper [1].
Three dimensional toric codes can either be constructed with the hypergraph product construction
or with a library, e.g., panqec [3].

### Code construction

The `code_construction` directory contains the code used to construct higher-dimensional hypergraph
product codes and used the `compute_distances.sh` script to automatically compute bounds on the
distances of the constructed codes using the GAP library QDistRnd.

## Utils

### Plotting

The `plotting` directory contains the code used to plot the results in the paper [1].

### Data Utils

We have implemented several utility functions as part of this package that might be of independent
interest for decoding simulations and data analysis. These are provided in the `data_utils` module
in the `utils` package.

## Dependencides

The used BP+OSD implementation, as well as our implementation of the SSMSA decoder are provided
in the LDPC2 package a preliminary beta version of which is available on Github [2].

## References

- [1]: L. Berent, T. Hillmann et al.: https://arxiv.org/abs/2311.01328
- [2]: J. Roffe et al. https://github.com/quantumgizmos/ldpc/tree/ldpc_v2
- [3]: E. Huang et al. https://github.com/panqec/panqec
