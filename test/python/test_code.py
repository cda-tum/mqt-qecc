"""Test the CSSCode class."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

from mqt.qecc import CSSCode, InvalidCSSCodeError

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt


@pytest.fixture
def rep_code() -> tuple[npt.NDArray[np.int8] | None, npt.NDArray[np.int8] | None]:
    """Return the parity check matrices for the repetition code."""
    hx = np.array([[1, 1, 0], [0, 0, 1]])
    hz = None
    return hx, hz


@pytest.fixture
def steane_code() -> tuple[npt.NDArray[np.int8], npt.NDArray[np.int8]]:
    """Return the check matrices for the Steane code."""
    hx = np.array([[1, 1, 1, 1, 0, 0, 0], [1, 0, 1, 0, 1, 0, 1], [0, 1, 1, 0, 1, 1, 0]])
    hz = hx
    return hx, hz


def test_invalid_css_codes() -> None:
    """Test that an invalid CSS code raises an error."""
    # Violates CSS condition
    hx = np.array([[1, 1, 1]])
    hz = np.array([[1, 0, 0]])
    with pytest.raises(InvalidCSSCodeError):
        CSSCode(distance=3, Hx=hx, Hz=hz)

    # Distances don't match
    hz = np.array([[1, 1, 0]])
    with pytest.raises(InvalidCSSCodeError):
        CSSCode(distance=3, Hx=hx, Hz=hz, x_distance=4, z_distance=1)

    # Checks not over the same number of qubits
    hz = np.array([[1, 1]])
    with pytest.raises(InvalidCSSCodeError):
        CSSCode(distance=3, Hx=hx, Hz=hz)

    # Invalid distance
    with pytest.raises(InvalidCSSCodeError):
        CSSCode(distance=-1, Hx=hx)

    # Checks not provided
    with pytest.raises(InvalidCSSCodeError):
        CSSCode(distance=3)


@pytest.mark.parametrize("checks", ["steane_code", "rep_code"])
def test_logicals(checks: tuple[npt.NDArray[np.int8] | None, npt.NDArray[np.int8] | None], request) -> None:  # type: ignore[no-untyped-def]
    """Test the logical operators of the CSSCode class."""
    hx, hz = request.getfixturevalue(checks)
    code = CSSCode(distance=3, Hx=hx, Hz=hz)
    assert code.Lx is not None
    assert code.Lz is not None
    assert code.Lx.shape[1] == code.Lz.shape[1] == hx.shape[1]
    assert code.Lx.shape[0] == code.Lz.shape[0]

    # assert that logicals anticommute
    assert code.Lx @ code.Lz.T % 2 != 0

    # assert that logicals commute with stabilizers
    if code.Hz is not None:
        assert np.all(code.Lx @ code.Hz.T % 2 == 0)
    if code.Hx is not None:
        assert np.all(code.Lz @ code.Hx.T % 2 == 0)


def test_errors(steane_code: tuple[npt.NDArray[np.int8], npt.NDArray[np.int8]]) -> None:
    """Test error detection and symdromes."""
    hx, hz = steane_code
    code = CSSCode(distance=3, Hx=hx, Hz=hz)
    e1 = np.array([1, 0, 0, 0, 0, 0, 0])
    e2 = np.array([0, 1, 0, 0, 1, 0, 0])
    e3 = np.array([0, 0, 0, 0, 0, 1, 1])
    e4 = np.array([0, 1, 1, 1, 0, 0, 0])

    assert np.array_equal(code.get_x_syndrome(e1), code.get_z_syndrome(e2))
    assert np.array_equal(code.get_x_syndrome(e2), code.get_z_syndrome(e2))

    x_syndrome_1 = code.get_x_syndrome(e1)
    x_syndrome_2 = code.get_x_syndrome(e2)
    x_syndrome_3 = code.get_x_syndrome(e3)
    x_syndrome_4 = code.get_x_syndrome(e4)

    assert np.array_equal(x_syndrome_1, x_syndrome_2)
    assert not np.array_equal(x_syndrome_1, x_syndrome_3)
    assert np.array_equal(x_syndrome_1, x_syndrome_4)

    # e1 and e2 have same syndrome but if we add them we get a logical error
    assert code.check_if_logical_x_error((e1 + e2) % 2)
    assert code.check_if_logical_z_error((e1 + e2) % 2)
    assert not code.stabilizer_eq_x_error(e1, e2)
    assert not code.stabilizer_eq_z_error(e1, e2)

    # e1 and e4 on the other hand do not induce a logical error because they are stabilizer equivalent
    assert not code.check_if_logical_x_error((e1 + e4) % 2)
    assert not code.check_if_logical_z_error((e1 + e4) % 2)
    assert code.stabilizer_eq_x_error(e1, e4)
    assert code.stabilizer_eq_z_error(e1, e4)


def test_steane(steane_code: tuple[npt.NDArray[np.int8], npt.NDArray[np.int8]]) -> None:
    """Test utility functions and correctness of the Steane code."""
    hx, hz = steane_code
    code = CSSCode(distance=3, Hx=hx, Hz=hz)
    assert code.n == 7
    assert code.k == 1
    assert code.distance == 3
    assert code.is_self_dual()

    x_paulis, z_paulis = code.stabs_as_pauli_strings()
    assert x_paulis is not None
    assert z_paulis is not None
    assert len(x_paulis) == len(z_paulis) == 3
    assert x_paulis == ["XXXXIII", "XIXIXIX", "IXXIXXI"]
    assert z_paulis == ["ZZZZIII", "ZIZIZIZ", "IZZIZZI"]

    x_log = code.x_logicals_as_pauli_string()
    z_log = code.z_logicals_as_pauli_string()
    assert x_log.count("X") == 3
    assert x_log.count("I") == 4
    assert z_log.count("Z") == 3
    assert z_log.count("I") == 4

    hx_reordered = hx[::-1, :]
    code_reordered = CSSCode(distance=3, Hx=hx_reordered, Hz=hz)
    assert code == code_reordered
