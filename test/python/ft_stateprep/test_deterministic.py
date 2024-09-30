"""Test synthesis and simulation of deterministic FT state preparation circuits."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from mqt.qecc import CSSCode
from mqt.qecc.ft_stateprep import DeterministicVerification, DeterministicVerificationHelper

if TYPE_CHECKING:
    import numpy.typing as npt
    from qiskit import QuantumCircuit
