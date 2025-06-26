from __future__ import annotations

import numpy as np
import pytest

from mqt.qecc.circuit_synthesis.synthesis_utils import CandidateAction, GaussianElimination

STEANE_CHECK_MATRIX = np.array([[1, 1, 0, 0, 1, 1, 0], [1, 0, 1, 0, 1, 0, 1], [0, 0, 0, 1, 1, 1, 1]], dtype=np.int8)


@pytest.fixture
def get_instance():
    """A basic fixture to provide a fresh GaussianElimination instance for each test."""
    return GaussianElimination(matrix=STEANE_CHECK_MATRIX.copy())


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


def test_apply_cnot_to_matrix_updates_matrix_and_costs_correctly(get_instance):
    """Verify the target column of the matrix is correctly updated after a CNOT."""
    expected_matrix = np.array([[1, 1, 0, 0, 0, 1, 0], [1, 0, 1, 0, 0, 0, 1], [0, 0, 0, 1, 1, 1, 1]], dtype=np.int8)
    get_instance._apply_cnot_to_matrix(i=0, j=4)
    np.testing.assert_array_equal(get_instance.matrix, expected_matrix)
    expected_instance = GaussianElimination(matrix=expected_matrix)
    expected_costs = expected_instance.costs
    np.testing.assert_array_equal(get_instance.costs, expected_costs)


# This is our new test
@pytest.mark.parametrize(
    (
        "scenario",
        "initial_used_cols",
        "costs_are_positive",
        "expected_return",
        "expect_reset_called",
        "expected_final_used_cols",
    ),
    [
        # Scenario 1: Stagnated and no used columns -> should call triangular_reset
        ("local_minimum", [], True, True, True, []),
        # Scenario 2: Stagnated with used columns -> should clear used_columns
        ("parallel_deadend", [1, 2], True, True, False, []),
        # Scenario 3: Not stagnated (negative cost exists) -> should do nothing
        ("not_stagnated", [1, 2], False, False, False, [1, 2]),
    ],
)
def test_handle_stagnation(
    get_instance,
    monkeypatch,
    scenario,
    initial_used_cols,
    costs_are_positive,
    expected_return,
    expect_reset_called,
    expected_final_used_cols,
):
    """Verify that _handle_stagnation behaves correctly in all scenarios."""
    costs_unused = np.array([[0, 1], [1, 0]]) if costs_are_positive else np.array([[0, -1], [1, 0]])
    get_instance.used_columns = initial_used_cols

    reset_call_log = []  # This list will act as our spy.
    # We replace the real, complex _triangular_reset with a fake function
    # that just appends to our log list when it's called.
    monkeypatch.setattr(get_instance, "_triangular_reset", lambda: reset_call_log.append(True))

    return_value = get_instance._handle_stagnation(costs_unused)

    assert return_value == expected_return
    assert get_instance.used_columns == expected_final_used_cols

    if expect_reset_called:
        assert len(reset_call_log) == 1
    else:
        assert len(reset_call_log) == 0


def test_basic_elimination_integration_with_known_matrix(get_instance):
    """Tests the full basic_elimination method from start to finish.

    It verifies that the correct sequence of eliminations is produced
    and that the matrix is fully reduced at the end.
    """
    expected_final_matrix = np.array(
        [[0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0]], dtype=np.int8
    )

    expected_eliminations = [(0, 4), (1, 5), (2, 6), (1, 0), (3, 4), (5, 6), (0, 2), (3, 5)]

    get_instance.basic_elimination()

    # 1. Did it produce the correct sequence of operations?
    assert get_instance.eliminations == expected_eliminations

    # 2. Is the algorithm correctly reporting that it's finished?
    assert get_instance.is_reduced() is True

    # 3. Is the final matrix what we expect?
    np.testing.assert_array_equal(get_instance.matrix, expected_final_matrix)
