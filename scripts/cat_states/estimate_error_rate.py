"""Estimate Cat State Preparation Logical Error Rate."""

from __future__ import annotations

import argparse

from mqt.qecc.circuit_synthesis.cat_states import CatStatePreparationExperiment, cat_state_balanced_tree

permutations = {4: None, 8: [7, 1, 2, 3, 4, 5, 6, 0], 16: [0, 1, 6, 10, 13, 3, 5, 15, 2, 8, 11, 14, 4, 9, 12, 7]}


def main() -> None:
    """Run the logical error rate estimation for cat state preparation circuits."""
    parser = argparse.ArgumentParser(description="Estimate logical error rate for cat state preparation circuits")

    parser.add_argument("-p", "--p_error", type=float, help="Physical error rate")
    parser.add_argument("-w", "--weight", type=int, help="Weight of the cat state. Can be 4, 8 or 16")
    parser.add_argument(
        "-n", "--n_samples", type=int, default=int(100e9), help="Number of samples to estimate logical error rate"
    )

    args = parser.parse_args()
    if args.weight not in {4, 8, 16}:
        msg = "Weight must be 4, 8 or 16"
        raise ValueError(msg)

    w = args.weight
    p = args.p_error
    n = args.n_samples

    circ = cat_state_balanced_tree(w)
    circ_anc = cat_state_balanced_tree(w)
    perm = permutations[w]

    experiment = CatStatePreparationExperiment(circ, circ_anc, perm)
    ra, ra_error, error_rates, error_rates_error = experiment.sample_cat_state(p, n)

    print(";".join([str(ra), str(ra_error), str(error_rates), str(error_rates_error)]))


if __name__ == "__main__":
    main()
