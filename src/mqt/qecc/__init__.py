"""MQT QECC library.

This file is part of the MQT QECC library released under the MIT license.
See README.md or go to https://github.com/cda-tum/qecc for more information.
"""

from __future__ import annotations

from ._version import version as __version__
from .analog_information_decoding.simulators.analog_tannergraph_decoding import AnalogTannergraphDecoder, AtdSimulator
from .analog_information_decoding.simulators.quasi_single_shot_v2 import QssSimulator
from .pyqecc import (
    Code,
    Decoder,
    DecodingResult,
    DecodingResultStatus,
    DecodingRunInformation,
    GrowthVariant,
    UFDecoder,
    UFHeuristic,
    apply_ecc,
    sample_iid_pauli_err,
)

__all__ = [
    "__version__",
    "Code",
    "Decoder",
    "UFHeuristic",
    "UFDecoder",
    "GrowthVariant",
    "DecodingResult",
    "DecodingResultStatus",
    "DecodingRunInformation",
    "sample_iid_pauli_err",
    "apply_ecc",
    "AnalogTannergraphDecoder",
    "AtdSimulator",
    # "SoftInfoDecoder",
    "QssSimulator",
]
