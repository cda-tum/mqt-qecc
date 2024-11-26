Encoder Circuit Synthesis for CSS Codes
=======================================

QECC provides functionality to synthesize encoder circuits for arbitrary :math:`[[n, k, d]]` quantum CSS codes.

    .. currentmodule:: mqt.qecc.circuit_synthesis

Non-fault tolerant state preparation circuits can be synthesized using :func:`depth_optimal_encoding_circuit`, :func:`gate_optimal_encoding_circuit` and :func:`heuristic_encoding_circuit`.

    .. autofunction:: depth_optimal_encoding_circuit

    .. autofunction:: gate_optimal_encoding_circuit

    .. autofunction:: heuristic_encoding_circuit
