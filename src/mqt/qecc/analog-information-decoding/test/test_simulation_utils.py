import math

import numpy as np
from scipy.sparse import csr_matrix, coo_matrix

import utils
from utils import simulation_utils


def test_check_logical_err_h() -> None:
    H = np.array(
        [
            [1, 0, 0, 1, 0, 1, 1],
            [0, 1, 0, 1, 1, 0, 1],
            [0, 0, 1, 0, 1, 1, 1]
        ])
    # check with logical
    estimate = np.array([1, 0, 0, 0, 0, 0, 1])
    assert simulation_utils.check_logical_err_h(H, np.array([0, 0, 0, 0, 1, 1, 0]), estimate) == True
    #
    # check with stabilizer
    estimate2 = np.array([0, 0, 0, 0, 0, 0, 1])
    assert simulation_utils.check_logical_err_h(H, np.array([1, 1, 1, 0, 0, 0, 0]), estimate2) == False

    # check with all zeros
    estimate3 = np.array([0, 0, 0, 0, 0, 0, 0])
    assert simulation_utils.check_logical_err_h(H, np.array([0, 0, 0, 0, 0, 0, 0]), estimate3) == False


def test_is_logical_err() -> None:
    # check with logical
    Lsc = np.array(
        [[1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    residual = np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    assert simulation_utils.is_logical_err(Lsc, residual) == True

    # check with stabilizer
    residual2 = np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    assert simulation_utils.is_logical_err(Lsc, residual2) == False

    # check with all zeros
    residual2 = np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    assert simulation_utils.is_logical_err(Lsc, residual2) == False

    # check with non-min weight logical
    residual3 = np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    assert simulation_utils.is_logical_err(Lsc, residual3) == True


def test_get_analog_llr() -> None:
    analog_syndr = np.array([0.5, 0, 0, -1, 0, 1])
    sigma = 0.8

    assert np.allclose(simulation_utils.get_analog_llr(analog_syndr, sigma),
                       np.array([1.5625, 0., 0., -3.125, 0., 3.125]))

    sigma = 0.1
    assert np.allclose(simulation_utils.get_analog_llr(analog_syndr, sigma), np.array([100, 0., 0., -200., 0., 200.]))

def test_generate_err() -> None:
    # no errors
    p = 0.0
    n = 10
    ch = np.ones(n) * p
    channel = [np.copy(ch), np.copy(ch), np.copy(ch)]
    residual = [np.zeros(n).astype(np.int32), np.zeros(n).astype(np.int32)]

    expected = [np.zeros(n).astype(np.int32), np.zeros(n).astype(np.int32)]
    assert np.array_equal(simulation_utils.generate_err(n, channel, residual), expected)

    residual[0][0] = 1
    residual[1][0] = 1

    expected = [np.copy(residual[0]), np.copy(residual[1])]
    res = simulation_utils.generate_err(n, channel, residual)
    assert np.array_equal(res[0], expected[0]) and np.array_equal(res[1], expected[1])


def test_get_sigma_from_syndr_er() -> None:
    ser = 0.1

    assert math.ceil(simulation_utils.get_sigma_from_syndr_er(ser)) == math.ceil(0.780304146072379)
    ser = 1.0
    assert math.ceil(simulation_utils.get_sigma_from_syndr_er(ser)) == -0.0
    ser = 0.0
    assert math.ceil(simulation_utils.get_sigma_from_syndr_er(ser)) == 0.0


def test_get_error_rate_from_sigma() -> None:
    sigma = 0.3
    assert np.isclose([simulation_utils.get_error_rate_from_sigma(sigma)], [0.00042906])
    sigma = 0.5
    assert np.isclose([simulation_utils.get_error_rate_from_sigma(sigma)], [0.02275])
    sigma = 0.0
    assert simulation_utils.get_error_rate_from_sigma(sigma) == 0.0


def test_get_virtual_check_init_vals() -> None:
    noisy_syndr = np.array([0.5, 0, 0, -1, 0, 100])
    sigma = 0.8

    assert np.allclose(simulation_utils.get_virtual_check_init_vals(noisy_syndr, sigma),
                       np.array([1.73288206e-001, 5.00000000e-001, 5.00000000e-001, 4.20877279e-002, 5.00000000e-001,
                                 1.91855567e-136]))

    sigma = 0.1
    assert np.allclose(simulation_utils.get_virtual_check_init_vals(noisy_syndr, sigma),
                       np.array([3.72007598e-44, 5.00000000e-01, 5.00000000e-01, 1.38389653e-87, 5.00000000e-01,
                                 0.00000000e+00]))
    sigma = 0.0
    res = simulation_utils.get_virtual_check_init_vals(noisy_syndr, sigma)

    assert math.isnan(res[1]) and res[0] == 0.


def test_generate_syndr_err() -> None:
    channel = np.array([0.0, 1.0])

    assert np.array_equal(simulation_utils.generate_syndr_err(channel), np.array([0., 1.]))


def test_get_noisy_analog_syndr() -> None:
    perfect_s = np.array([1, 0])
    sigma = 0.0

    assert np.array_equal(simulation_utils.get_noisy_analog_syndrome(perfect_s, sigma), np.array([-1., 1.]))


def test_err_chnl_setup() -> None:
    p = 0.1
    bias = [1., 1., 1.]
    n = 10
    ar = np.ones(n) * p / 3
    exp = [np.copy(ar), np.copy(ar), np.copy(ar)]
    res = simulation_utils.error_channel_setup(p, bias, n)

    assert np.array_equal(res, exp)

    bias = [1., 0., 0.]
    ar = np.ones(n) * p
    exp = [np.copy(ar), np.zeros(n), np.zeros(n)]
    res = simulation_utils.error_channel_setup(p, bias, n)

    assert np.array_equal(res, exp)

    bias = [1., 1., 0.]
    ar = np.ones(n) * p / 2
    exp = [np.copy(ar), np.copy(ar), np.zeros(n)]
    res = simulation_utils.error_channel_setup(p, bias, n)

    assert np.array_equal(res, exp)

    bias = [np.inf, 0., 0.]
    ar = np.ones(n) * p
    exp = [np.copy(ar), np.zeros(n), np.zeros(n)]
    res = simulation_utils.error_channel_setup(p, bias, n)

    assert np.array_equal(res, exp)


def test_build_ss_pcm() -> None:
    H = np.array([
        [1, 1, 0, 0],
        [0, 1, 1, 0],
        [0, 0, 1, 1],
        [1, 0, 0, 1]]
    )
    M = np.array([
        [1, 1, 0, 0],
        [0, 1, 1, 0],
        [0, 0, 1, 1],
        [1, 0, 0, 1]]
    )
    id_r = np.identity(M.shape[1])
    zeros = np.zeros((M.shape[0], H.shape[1]))
    exp = np.block([[H, id_r], [zeros, M]]).astype(np.int32)
    assert np.array_equal(utils.build_single_stage_pcm(H, M), exp)


def test_get_signed_from_binary() -> None:
    binary = np.array([1, 0, 0, 1, 0, 1])
    exp = np.array([-1, 1, 1, -1, 1, -1])

    assert np.array_equal(utils.get_signed_from_binary(binary), exp)


def test_get_binary_from_analog() -> None:
    exp = np.array([1, 0, 0, 1, 0, 1])
    analog = np.array([-1.0, 3., 1., -1., 1., -2])

    assert np.array_equal(utils.get_binary_from_analog(analog), exp)
