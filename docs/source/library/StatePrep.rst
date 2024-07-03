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
