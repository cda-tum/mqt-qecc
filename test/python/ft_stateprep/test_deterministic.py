"""Test synthesis and simulation of deterministic FT state preparation circuits."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from ldpc import mod2
from qsample import noise

from mqt.qecc import CSSCode
from mqt.qecc.ft_stateprep import DeterministicVerificationHelper, NoisyDFTStatePrepSimulator, heuristic_prep_circuit

if TYPE_CHECKING:
    import numpy.typing as npt

    from mqt.qecc.ft_stateprep import DeterministicVerification, StatePrepCircuit

# Simulation parameters
err_params = {"q": [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1]}
err_model = noise.E1_1
shots_dss = 2000
p_max = {"q": 0.01}
L = 3


@pytest.fixture
def steane_code_sp_plus() -> StatePrepCircuit:
    """Return a non-ft state preparation circuit for the Steane code."""
    steane_code = CSSCode.from_code_name("Steane")
    sp_circ = heuristic_prep_circuit(steane_code, zero_state=False)
    sp_circ.compute_fault_sets()
    return sp_circ


@pytest.fixture
def css_11_1_3_code_sp() -> StatePrepCircuit:
    """Return a non-ft state preparation circuit for the 11_1_3 code."""
    check_matrix = np.array([
        [1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0],
        [0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
        [0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0],
        [0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1],
    ])
    code = CSSCode(distance=3, Hx=check_matrix, Hz=check_matrix)
    sp_circ = heuristic_prep_circuit(code)
    sp_circ.compute_fault_sets()
    return sp_circ


def in_span(m: npt.NDArray[np.int_], v: npt.NDArray[np.int_]) -> bool:
    """Check if a vector is in the row space of a matrix."""
    return bool(mod2.rank(np.vstack((m, v))) == mod2.rank(m))


def assert_statistics(
    verify: DeterministicVerification,
    num_ancillas_verification: int,
    num_cnots_verification: int,
    num_ancillas_correction: int,
    num_cnots_correction: int,
    num_ancillas_hooks: int = 0,
    num_cnots_hooks: int = 0,
    num_ancillas_hook_corrections: int = 0,
    num_cnots_hook_corrections: int = 0,
) -> None:
    """Assert that the statistics of a deterministic verification are correct."""
    assert verify.num_ancillas_verification() == num_ancillas_verification
    assert verify.num_cnots_verification() == num_cnots_verification
    assert verify.num_ancillas_correction() <= num_ancillas_correction
    assert verify.num_cnots_correction() <= num_cnots_correction
    assert verify.num_ancillas_hooks() == num_ancillas_hooks
    assert verify.num_cnots_hooks() == num_cnots_hooks
    assert verify.num_ancillas_hook_corrections() == num_ancillas_hook_corrections
    assert verify.num_cnots_hook_corrections() == num_cnots_hook_corrections


def assert_stabs(verify: DeterministicVerification, code: CSSCode, z_stabs: bool) -> None:
    """Assert that the measurement stabs of a deterministic verification are correct."""
    checks = np.vstack((code.Hz, code.Lz))
    checks_other = np.vstack((code.Hx, code.Lx))
    if not z_stabs:
        checks, checks_other = checks_other, checks

    for stab in verify.stabs:
        assert in_span(checks, stab)
    for correction in verify.det_correction.values():
        stabs, _ = correction
        for stab in stabs:
            assert in_span(checks, stab)
    for hook in verify.hook_corrections:
        if not hook:
            continue
        for correction in hook.values():
            stabs, _ = correction
            for stab in stabs:
                assert in_span(checks_other, stab)


def assert_scaling(simulation_results: list[npt.NDArray[np.float64]]) -> None:
    """Assert that the logical error rates scales approximately quadratically."""
    dss_upper_bound = simulation_results[-2]
    x = np.log10(err_params["q"])
    y = np.log10(dss_upper_bound)
    m = np.diff(y) / np.diff(x)
    assert np.average(m[:4]) > 1.7


def test_11_1_3_det_verification(css_11_1_3_code_sp: StatePrepCircuit) -> None:
    """Test deterministic verification of the 11_1_3 code state preparation circuit."""
    verify_helper = DeterministicVerificationHelper(css_11_1_3_code_sp)
    verify_x, verify_z = verify_helper.get_solution()

    # as this is not optimal it might be possible that different verification result in different corrections
    assert_statistics(verify_x, 2, 8, 4, 14, 0, 0)
    assert_stabs(verify_x, css_11_1_3_code_sp.code, z_stabs=True)

    assert_statistics(verify_z, 1, 4, 1, 4, 1, 2, 1, 3)
    assert_stabs(verify_z, css_11_1_3_code_sp.code, z_stabs=False)

    # perform simulation
    simulator = NoisyDFTStatePrepSimulator(
        css_11_1_3_code_sp.circ, (verify_x, verify_z), css_11_1_3_code_sp.code, err_model
    )
    simulation_results = simulator.dss_logical_error_rates(err_params, p_max, L, shots_dss)
    assert_scaling(simulation_results)


def test_steane_det_verification(steane_code_sp_plus: StatePrepCircuit) -> None:
    """Test deterministic verification of the Steane code state preparation circuit."""
    verify_helper = DeterministicVerificationHelper(steane_code_sp_plus)
    verify_z_opt, verify_x_opt = verify_helper.get_solution()
    verify_z_global, verify_x_global = verify_helper.get_global_solution()

    # Check right statistics
    for verify_x, verify_z in zip((verify_x_opt, verify_x_global), (verify_z_opt, verify_z_global)):
        assert_statistics(verify_z, 1, 3, 1, 3, 0, 0)
        assert_stabs(verify_z, steane_code_sp_plus.code, z_stabs=False)

        # second verification is trivial
        assert verify_x.num_ancillas_total() == 0
        assert verify_x.num_cnots_total() == 0

        # perform simulation
        simulator = NoisyDFTStatePrepSimulator(
            steane_code_sp_plus.circ, (verify_z, verify_x), steane_code_sp_plus.code, err_model, False
        )
        simulation_results = simulator.dss_logical_error_rates(err_params, p_max, L, shots_dss)
        assert_scaling(simulation_results)
