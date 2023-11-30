"""Tests for the memory experiment simulator."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from mqt.qecc.analog_information_decoding.simulators.memory_experiment_v2 import build_multiround_pcm, move_syndrome

if TYPE_CHECKING:
    from numpy._typing import NDArray


def test_build_mr_pcm() -> None:
    """Test build_multiround_pcm function."""
    pcm: NDArray[np.int32] = np.array([[1, 1, 0], [0, 1, 1]]).astype(np.int32)
    mr_pcm = build_multiround_pcm(pcm, 1)
    np.zeros((2, 3))
    np.identity(2)
    r1 = np.hstack([pcm, np.zeros(pcm.shape), np.identity(2), np.zeros((2, 2))])
    r2 = np.hstack([np.zeros(pcm.shape), pcm, np.identity(2), np.identity(2)])
    expected = np.vstack((r1, r2))

    assert np.array_equal(mr_pcm.toarray(), expected)


def test_move_syndrome() -> None:
    """Test move_syndrome function."""
    # three bit syndrome over 4 rounds
    syndr = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    res = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
    assert np.array_equal(res, move_syndrome(syndr))

    syndr = np.array([[0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]])
    res = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [1, 1, 0, 0]])
    assert np.array_equal(res, move_syndrome(syndr))

    syndr = np.array([[0, 0, 1, 1, 0, 0], [0, 0, 1, 1, 0, 0], [0, 0, 1, 1, 0, 0]])
    res = np.array([[1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0]])
    assert np.array_equal(res, move_syndrome(syndr))


def test_decode_multiround() -> None:
    """Test decoding of multiround syndrome."""
    # three bit syndrome over 4 rounds
    syndr = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    res = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
    assert np.array_equal(res, move_syndrome(syndr))

    syndr = np.array([[0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]])
    res = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [1, 1, 0, 0]])
    assert np.array_equal(res, move_syndrome(syndr))

    syndr = np.array([[0, 0, 1, 1, 0, 0], [0, 0, 1, 1, 0, 0], [0, 0, 1, 1, 0, 0]])
    res = np.array([[1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0]])
    assert np.array_equal(res, move_syndrome(syndr))
