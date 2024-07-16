#!/usr/bin/env bash

# set number of parallel processes
N=16
export export LC_NUMERIC=C
# run the simulations
 parallel -j $N mqt.qecc.cc-decoder {1} {2} --nr_sims 10000 --results_dir ./results/maxsat --decoder maxsat ::: $(seq 3 2 21) ::: $(seq 0.001 0.001 0.175)
parallel -j $N mqt.qecc.cc-decoder {1} {2} --nr_sims 10000 --results_dir ./results/tn --decoder tn ::: $(seq 3 2 21) ::: $(seq 0.001 0.001 0.175)
