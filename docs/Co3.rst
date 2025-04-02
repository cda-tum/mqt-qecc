Compilation beyond the Surface Code ``co3``
===========================

This submodule contains an elementary routing routine for CNOT + T compilation on a hexagonal routing graph.
Moreover the routing assumes that each accessible boundary can host both Z and X operators for lattice surgery.
Thus, if a substrate cannot use all boundaries for lattice surgery (e.g. folded surface code), one has to add adaptions to the routing etc and choose layouts wisely.
Hence, this code is only viable for color codes only. Adaptions needed to incorporate valid paths for e.g. the folded surface code substrate.

Layouts
#######

Basic layouts (sparse, pair, row, hex) can be generated automatically in the ``HexagonalLattice`` class.
However, to place factories in a suitable and controlled manner, it is advisable to construct layouts manually as shown in ``plots/construct_layouts.ipynb``.
A layout describes which nodes on the routing graph are used as logical data qubits and factory locations. The remainder is the routing ancilla space.
The mapping of logical qubit labels onto those chosen data qubit locations on the graph is another task.


Compilation of given layout and qubit label allocation
######################################################

The higher level compilation follows a simple greedy routine for solving the VDP problem. We greedily extended this to include paths to factories as well.
Note that the class ``ShortestFirstRouterTGatesDyn`` and particularly the method ``find_total_vdp_layers_dyn`` should be used to perform routing.


Optimization of qubit label allocation by Hill Climbing
#######################################################

Once chosen a layout, one can optimize the qubit label allocation. This is crucial to exploit parallelism of the original circuit.
The class ``HillClimbing`` performs a simple hill climbing routine to optimize the qubit label mapping based on a heuristic metric which computes the initial crossing of shortest paths as well as a more reliable (yet expensive) metric which computes the routing for each Hill climbing iteration and directly aims to reduce the resulting layers.


Microscopic Details
###################

We consider two microscopic substrates, both leading to a hexagonal routing graph.
First, the class ``SnakeBuilderSTDW`` builds stabilizers and subsets of stabilizers to perform logical meausrements for the color code connected by semi transparent domain walls.
The class ``SnakeBuilderSC`` builds the surface code snakes required to perform lattice surgery between logical folded surface codes. However, this can only display snakes where you can embed the snake in 2d.
A notebook with example constructions can be found in ``/microscopic/snake_examples.ipynb``.
