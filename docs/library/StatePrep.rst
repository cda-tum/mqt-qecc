Fault tolerant state preparation
================================

QECC provides functionality to synthesize and simulate state preparation circuits for logical basis states for arbitrary :math:`[[n, k, d]]` quantum CSS codes.

    .. currentmodule:: mqt.qecc.ft_stateprep

Non-fault tolerant state preparation circuits can be synthesized using :func:`depth_optimal_prep_circuit`, :func:`gate_optimal_prep_circuit` and :func:`heuristic_prep_circuit`.

    .. autofunction:: depth_optimal_prep_circuit

    .. autofunction:: gate_optimal_prep_circuit

    .. autofunction:: heuristic_prep_circuit

These methods return a :class:`StatePrepCircuit` from which the circuit can be obtained as a qiskit :code:`QuantumCircuit` object via the :code:`circ` member. The :class:`StatePrepCircuit` class contains methods for computing the state preparation circuit's fault set. :class:`StatePrepCircuit` are given as input to the verification circuit synthesis methods which add verification measurements to the circuit such that postselection on +1 measurement results of these circuit outputs a state with a logical error rate on the order of :math:`O(p^{\frac{d}{2}})`.

Gate-optimal verification circuits can be generated using the :func:`gate_optimal_verification_circuit` method. This method uses the SMT solver `Z3 <https://github.com/Z3Prover/z3>`_ for iteratively searching for better verification circuits. The search is guided by a :code:`min_timeout` and :code:`max_timeout` parameter. Initially, the search is only allowed to continue for the minimal amount of time. This time budget is successively increased until a solution is found. At this point the maximum number of CNOTs is reduced until the synthesis takes :code:`max_timeout` time or if the SAT solver manages to prove that no smaller circuit exists.

    .. autofunction:: gate_optimal_verification_circuit

If the optimal synthesis takes too long, the :func:`heuristic_verification_circuit` method can be used. This method reduces the synthesis of state preparation circuits to a set cover problem. Quality of solution can be traded with performance via the :code:`max_covering_sets` parameter. The smaller this parameter is set the lower the number of sets from which a covering is obtained.

    .. autofunction:: heuristic_verification_circuit

State preparation circuits can be simulated using the :class:`NoisyNDFTStatePrepSimulator` class.

    .. autoclass:: NoisyNDFTStatePrepSimulator
        :members:

:math:`d<5` codes
---------------------------------------------------------
For small distance codes QECC provides functionality to not only synthesize state preparation circuits based on post-selection but also preparation protocols that yield the logical Pauli state deterministically.
Such a deterministic protocol consists of three parts: (1) a (non-deterministic) verification, (2) a deterministic correction if one of the verification measurements yields a non-trivial result, and (3) a hook correction for the case that one of the verification measurement flags indicates a hook error.
To facilitate the handling (e.g. collecting statistics regarding ancilla and CNOT count) QECC provides the :class:`DeterministicVerification` class.

    .. autoclass:: DeterministicVerification
        :members:

The :class:`DeterministicVerification` for a certain :class:`StatePrepCircuit` can be obtained using the wrapper class :class:`DeterministicVerificationHelper` which provides the two functions :func:`get_solution` and :func:`get_global_solution` where the latter optimizes over all possible(*potentially exponentially many*)  non-deterministic verifications to find the global optimum (recommended only for codes with small qubit number).
The two :class:`DeterministicVerification` objects returned correspond to the two layers of verification (X and Z errors).

    .. autoclass:: DeterministicVerificationHelper
        :members:
The resulting protocol consisting of the two-layered verification can be converted to a `Qsample <https://github.com/dpwinter/qsample>`_ protocol and simulated using Dynamic Subset Sampling :cite:labelpar:`heussenDynamicalSubsetSampling2024`.
The protocol construction and simulation is done automatically by the :class:`NoisyDFTStatePrepSimulator` class where the returned numbers correspond to the upper and lower bounds with corresponding errors as described at the `Qsample Documentation <https://dpwinter.github.io/qsample/>`_.

    .. autoclass:: NoisyDFTStatePrepSimulator
        :members:
