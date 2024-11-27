"""Methods for working with graphical code representations based on the ZX-calculus."""

from __future__ import annotations

from itertools import starmap
from typing import TYPE_CHECKING

import numpy as np
import pyzx as zx
from ldpc import mod2

if TYPE_CHECKING:  # pragma: no cover
    import numpy.typing as npt

    from ..codes import CSSCode


class ZXCSSEncoder:
    """Class representing the encoding isometry of a CSS code with the ZX-calculus."""

    def __init__(self, diag: zx.Graph, logical_in: list[zx.VT], stab_in: list[zx.VT]) -> None:
        """Initialize the encoder from a ZX-diagram.

        Args:
            diag: The ZX-diagram representing the encoding isometry.
            logical_in: The input verices for logical operators. Must be of type zx.VertexType.Boundary.
            stab_in: The input vertices for stabilizer operators. Must be of type zx.VertexType.Boundary.
        """
        self.diag = diag
        self.logicals = logical_in
        self.stabs = stab_in

    @classmethod
    def from_code(cls, code: CSSCode, z_x_form: bool = True) -> ZXCSSEncoder:
        """Create a ZX-diagram from a CSS code.

        Args:
            code: The CSS code.
            z_x_form: Whether to use the Z-X form of the code.
        """
        if z_x_form:
            stabs = code.Hx
            logicals = code.Lx
        else:
            stabs = code.Hz
            logicals = code.Lz

        return cls.from_stabilizers_and_logicals(stabs, logicals, z_x_form)

    @classmethod
    def from_stabilizers_and_logicals(
        cls, stabs: npt.NDArray[np.int8], logicals: npt.NDArray[np.int8], z_x_form: bool = True
    ) -> ZXCSSEncoder:
        """Create a ZX-diagram from the stabilizers and logicals of a CSS code.

        Args:
            stabs: The stabilizers of the code.
            logicals: The logical operators of the code.
            z_x_form: Whether to use the Z-X form of the code or the X-Z form.
        """
        g = zx.Graph(backends="simple")
        rank = mod2.rank(stabs)
        if rank != stabs.shape[0]:
            msg = "Stabilizers must be linearly independent."
            raise ValueError(msg)

        n_in_stabs = stabs.shape[0]
        n_in_logicals = logicals.shape[0]
        n_out_stabs = logicals.shape[0] + stabs.shape[1]
        n_out_logicals = n_in_logicals

        in_type = zx.VertexType.Z if z_x_form else zx.VertexType.X
        out_type = zx.VertexType.X if z_x_form else zx.VertexType.Z

        v_in_stabs = g.add_vertices(n_in_stabs, vtype=in_type)
        v_out_stabs = g.add_vertices(n_out_stabs, vtype=out_type)
        v_in_logicals = g.add_vertices(n_in_logicals, vtype=in_type)
        v_out_logicals = g.add_vertices(n_out_logicals, vtype=out_type)
        v_boundary_stabs = g.add_vertices(n_out_stabs, vtype=zx.VertexType.BOUNDARY)
        v_boundary_logicals = g.add_vertices(n_out_logicals, vtype=zx.VertexType.BOUNDARY)

        # connect outputs to boundary
        g.add_edges(starmap(g.edge, zip(v_out_stabs, v_boundary_stabs)))
        g.add_edges(starmap(g.edge, zip(v_out_logicals, v_boundary_logicals)))

        # Add interior edges for stabilizers
        g.add_edges(g.edge(v_in_stabs[row], v_out_stabs[col]) for row, col in cls._edge_list_from_adj_matrix(stabs))

        # Add interior edges for logicals
        g.add_edges(
            g.edge(v_in_logicals[row], v_out_logicals[col]) for row, col in cls._edge_list_from_adj_matrix(logicals)
        )

        return cls(g, v_in_logicals, v_in_stabs)

    @staticmethod
    def _edge_list_from_adj_matrix(adj: npt.NDArray[np.int8]) -> list[tuple[int, int]]:
        """Create a list of edges from an adjacency matrix."""
        rows, cols = np.triu_indices_from(adj, k=1)
        mask = adj[rows, cols] != 0
        return list(zip(rows[mask], cols[mask]))
