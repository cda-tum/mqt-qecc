from __future__ import annotations

import numpy as np
import pytest

from mqt.qecc.circuit_synthesis.synthesis_utils import CandidateAction, GaussianElimination

SIMPLE_MATRIX = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype=np.int8)


@pytest.fixture
def get_instance():
    """A basic fixture to provide a fresh GaussianElimination instance for each test."""
    return GaussianElimination(matrix=SIMPLE_MATRIX.copy())


# Use parametrize to test all logical branches
@pytest.mark.parametrize(
    ("cost", "backtrack_required", "used_cnots", "failed_cnots", "expected_action"),
    [
        # Case 1: Greedy step, should evaluate
        (-1, False, [], [], CandidateAction.EVALUATE),
        # Case 2: Already used, should skip
        (-1, False, [(0, 1)], [], CandidateAction.SKIP),
        # Case 3: Already failed, should skip
        (-1, False, [], [(0, 1)], CandidateAction.SKIP),
        # Case 4: First time hitting non-negative cost, should restart search
        (0, False, [], [], CandidateAction.RESTART_SEARCH),
        (1, False, [], [], CandidateAction.RESTART_SEARCH),
        # Case 5: Already in backtrack mode and hitting non-negative cost, should trigger backtrack
        (0, True, [], [], CandidateAction.TRIGGER_BACKTRACK),
        (1, True, [], [], CandidateAction.TRIGGER_BACKTRACK),
    ],
)
def test_get_candidate_action(get_instance, cost, backtrack_required, used_cnots, failed_cnots, expected_action):
    """Verify that _get_candidate_action returns the correct action based on state."""
    # Arrange: Set up the instance state for the test case
    get_instance.backtrack_required = backtrack_required
    get_instance.used_cnots = used_cnots
    get_instance.failed_cnots = failed_cnots

    # Act: Call the method
    action = get_instance._get_candidate_action(i=0, j=1, costs=cost)

    # Assert: Check the result
    assert action == expected_action
