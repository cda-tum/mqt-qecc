"""Test the Hill Climbing."""

from __future__ import annotations

import mqt.qecc.co3.utils.hill_climber as co


def test_neighborhood():
    """Checks an example neighborhood."""
    layout_type = "row"
    m = 4
    n = 4
    metric = "crossing"
    circuit = [(5, 4), (3, 1), (2, 0)]
    max_restarts = None
    max_iterations = None
    hc = co.HillClimbing(max_restarts, max_iterations, circuit, layout_type, m, n, metric)
    layout = {1: (1, 2), 5: (1, 3), 3: (2, 2), 0: (2, 3), 4: (3, 2), 2: (3, 3)}

    neighborhood = hc.gen_neighborhood(layout)
    aim_neighborhood = [
        {1: (1, 2), 5: (3, 2), 3: (2, 2), 0: (2, 3), 4: (1, 3), 2: (3, 3)},
        {1: (2, 2), 5: (1, 3), 3: (1, 2), 0: (2, 3), 4: (3, 2), 2: (3, 3)},
        {1: (1, 2), 5: (1, 3), 3: (2, 2), 0: (3, 3), 4: (3, 2), 2: (2, 3)},
    ]

    error_message = "Generates unexpected neighborhood."
    assert neighborhood == aim_neighborhood, error_message


def test_crossing_metric():
    """Checks an example cost."""
    layout_type = "row"
    m = 4
    n = 4
    metric = "crossing"
    circuit = [(5, 4), (3, 1), (2, 0)]
    max_restarts = None
    max_iterations = None
    hc = co.HillClimbing(max_restarts, max_iterations, circuit, layout_type, m, n, metric)
    layout = {1: (1, 2), 5: (1, 3), 3: (2, 2), 0: (2, 3), 4: (3, 2), 2: (3, 3), "factory_positions": []}

    cost = hc.evaluate_solution(layout)
    expected_cost = 2

    error_message = "Generates unexpected crossing cost value."
    assert cost == expected_cost, error_message
