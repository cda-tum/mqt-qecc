from __future__ import annotations

import argparse
from unittest import mock

import mqt.qecc.ecc_framework_qiskit_wrapper


@mock.patch(
    "argparse.ArgumentParser.parse_args",
    return_value=argparse.Namespace(
        m="D",
        p=0.001,
        n=1000,
        s=1,
        f="ghz02.qasm",
        e=None,
        fs="aer_simulator_extended_stabilizer",
        ecc="Q7Steane",
        fq=0,
    ),
)
def test_command(mock_args):
    mqt.qecc.ecc_framework_qiskit_wrapper.main()


if __name__ == "__main__":
    test_command()
