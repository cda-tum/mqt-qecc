MaxSAT Decoder
=============

The MaxSAT package provides a MaxSAT decoder for decoding quantum codes. Currently the automatic simulation of
triangular color codes is supported. The decoder is based on the Z3 MaxSAT solver.

Usage
#####

The decoder is implemented in :func:`~mqt.qecc.cc_decoder` and simulations can be run automatically with the run method:

.. code-block:: python

  from mqt.qecc.cc_decoder import decoder

  d = 21  # distance of the triangular code to simulate
  p = 0.01  # (bit-flip) error rate
  nr_sims = 1000  # number of simulations to run
  decoder.run(distance=d, error_rate=p, nr_sims=nr_sims)

Currently the only supported error model is bit-flip errors with perfect syndrome measurements.

Furthermore, a command-line interface is provided to run simulations in parallel (e.g. using GNU Parallel)
.. code-block:: console

  $ parallel -j 8 mqt.qecc.cc-decoder {1} {2} --nr_sims 10000 --results_dir ./results/maxsat --decoder maxsat ::: $(seq 3 2 23) ::: $(seq 0.001 0.001 0.175)
