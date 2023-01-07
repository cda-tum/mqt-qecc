import mqt.qecc.ecc_framework_qiskit_wrapper
import pytest
from typing import TYPE_CHECKING, Any, cast


@pytest.fixture(params=['none',
                        'Id',
                        'Q3Shor',
                        'Q7Steane',
                        ])
def simulator(request: Any) -> str:
    return cast(str, request.param)


@pytest.fixture(params=['B', 'P', 'D'])
def noise_models(request: Any) -> str:
    return cast(str, request.param)


def test_with_stab_simulator(simulator, noise_models, script_runner):
    """Testing the script with different parameters"""
    ret = script_runner.run('ecc_framework_qiskit_wrapper',
                            '-m', noise_models,
                            '-p', '0.001',
                            '-n', '1000',
                            '-s', '1',
                            '-fs', 'aer_simulator_extended_stabilizer',
                            '-ecc', simulator,
                            '-f', 'ExampleCircuit.qasm',
                            )
    assert ret.success

    def test_amplitude_damping_simulator(simulator, noise_models, script_runner):
        """Testing the script with different parameters"""
        ret = script_runner.run('ecc_framework_qiskit_wrapper',
                                '-m', 'APD',
                                '-p', '0.001',
                                '-n', '1000',
                                '-s', '1',
                                '-fs', 'aer_simulat_extended_stabilizer',
                                '-ecc', 'Id',
                                '-f', 'ExampleCircuit.qasm',
                                )
        assert ret.success


def test_failing_simulators(script_runner):
    """Testing the script with unsupported ecc."""
    ret = script_runner.run('ecc_framework_qiskit_wrapper',
                            '-m', 'D',
                            '-p', '0.001',
                            '-n', '1000',
                            '-s', '1',
                            '-fs', 'aer_simulator_extended_stabilizer',
                            '-ecc', 'UnsupportedEcc',
                            '-f', 'ExampleCircuit.qasm',
                            )

    assert not ret.success
    assert 'No ECC found for' in ret.stderr
