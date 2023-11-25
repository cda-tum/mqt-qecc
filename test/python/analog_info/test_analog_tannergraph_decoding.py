from unittest import mock

import numpy as np
import pytest
from mqt.qecc.analog_information_decoding.simulators.analog_tannergraph_decoding import AnalogTannergraphDecoder
from mqt.qecc.analog_information_decoding.utils.data_utils import BpParams
from typing import TYPE_CHECKING

from mqt.qecc import AtdSimulator
from mqt.qecc.analog_information_decoding.utils import simulation_utils

if TYPE_CHECKING:
    from numpy._typing import NDArray


@pytest.fixture()
def error_channel():
    """Return error channel."""
    return np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])


@pytest.fixture()
def pcm():
    """Return parity check matrix."""
    return np.array([
        [1, 0, 0, 1, 0, 1, 1],
        [0, 1, 0, 1, 1, 1, 1],
        [0, 0, 1, 0, 1, 0, 1]
    ]).astype(np.int32)


@pytest.fixture()
def atd(error_channel, pcm) -> AnalogTannergraphDecoder:
    """Return distance of the hexagonal color code."""
    return AnalogTannergraphDecoder(
        pcm=pcm,
        bp_params=BpParams(osd_method="osd0"),
        error_channel=error_channel,
        sigma=0.1,
        ser=None,
    )


def test_set_analog_syndrome(atd: AnalogTannergraphDecoder, error_channel) -> None:
    analog_syndr = np.array([0.1, 0.1, 0.1])
    expected_virtual_inits = np.array([2.0611536181902108e-09, 2.0611536181902108e-09, 2.0611536181902108e-09])
    expected = np.hstack([error_channel, expected_virtual_inits])

    res = atd.decode(analog_syndr)
    assert np.array_equal(atd.bposd_decoder.channel_probs, expected)
    assert res is not None and len(res) == 10


def test_atd_simulator_data_error_channels_setup(pcm) -> None:
    per = 0.1
    sigma = 0.1
    sim = AtdSimulator(
        hx=pcm,
        hz=pcm,
        lx=np.array([]),
        lz=np.array([]),
        codename='test',
        data_err_rate=per,
        sigma=sigma,
        seed=666,
        bp_params=BpParams(osd_method="osd0"),
        decoding_method="atd",
    )
    expected_err_chnl = simulation_utils.error_channel_setup(
        error_rate=per,
        xyz_error_bias=[0.1, 0.1, 0.1],
        nr_qubits=pcm.shape[1],
    )
    assert np.array_equal(sim.x_decoder.error_channel, expected_err_chnl[0] + expected_err_chnl[1])
    assert np.array_equal(sim.z_decoder.error_channel, expected_err_chnl[2] + expected_err_chnl[1])


def test_atd_simulator_syndrome_error_channels_setup(pcm) -> None:
    per = 0.1
    sigma = 0.1
    ser = simulation_utils.get_error_rate_from_sigma(sigma)
    sim = AtdSimulator(
        hx=pcm,
        hz=pcm,
        lx=np.array([]),
        lz=np.array([]),
        codename='test',
        data_err_rate=per,
        sigma=sigma,
        seed=666,
        bp_params=BpParams(osd_method="osd0"),
        decoding_method="atd",
    )
    expec_chnl = simulation_utils.error_channel_setup(
        error_rate=ser,
        xyz_error_bias=[0.1, 0.1, 0.1],
        nr_qubits=1,
    )
    assert sim.syndr_err_rate == simulation_utils.get_error_rate_from_sigma(sigma=sigma)
    assert sim.x_sigma == simulation_utils.get_sigma_from_syndr_er(expec_chnl[0][0] + expec_chnl[1][0])
    assert sim.z_sigma == simulation_utils.get_sigma_from_syndr_er(expec_chnl[2][0] + expec_chnl[1][0])


def test_atd_simulator_syndrome_error_channels_setup(pcm) -> None:
    per = 0.1
    ser = 0.1
    sim = AtdSimulator(
        hx=pcm,
        hz=pcm,
        lx=np.array([]),
        lz=np.array([]),
        codename='test',
        data_err_rate=per,
        syndr_err_rate=ser,
        seed=666,
        bp_params=BpParams(osd_method="osd0"),
        decoding_method="atd",
    )
    expec_chnl = simulation_utils.error_channel_setup(
        error_rate=ser,
        xyz_error_bias=[0.1, 0.1, 0.1],
        nr_qubits=1,
    )
    assert sim.x_sigma == simulation_utils.get_sigma_from_syndr_er(expec_chnl[0][0] + expec_chnl[1][0])
    assert sim.z_sigma == simulation_utils.get_sigma_from_syndr_er(expec_chnl[2][0] + expec_chnl[1][0])

def test_single_sample(pcm)-> None:
    per = 0.1
    ser = 0.1
    sim = AtdSimulator(
        hx=pcm,
        hz=pcm,
        lx=np.array([1,1,1,1,1,1,1]),
        lz=np.array([1,1,1,1,1,1,1]),
        codename='test',
        data_err_rate=per,
        syndr_err_rate=ser,
        seed=666,
        bp_params=BpParams(osd_method="osd0"),
        decoding_method="atd",
    )

    assert sim.single_sample() is not None
    assert sim.x_bp_iterations is not None
    assert sim.z_bp_iterations is not None

def test_safe_results(pcm)-> None:
    per = 0.1
    ser = 0.1
    sim = AtdSimulator(
        hx=pcm,
        hz=pcm,
        lx=np.array([1,1,1,1,1,1,1]),
        lz=np.array([1,1,1,1,1,1,1]),
        codename='test',
        data_err_rate=per,
        syndr_err_rate=ser,
        seed=666,
        bp_params=BpParams(osd_method="osd0"),
        decoding_method="atd",
        output_path="./results"
    )
    with mock.patch('json.dump', return_value=True):
            res = sim.save_results(1,1,1)
    assert res is not None
    assert res["code_K"] == 3
    assert res["code_N"] == 7
