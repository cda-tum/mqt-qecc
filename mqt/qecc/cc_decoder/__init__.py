"""Contains the implementation of the tensor network decoder for the hexagonal color code."""
from mqt.qecc.cc_decoder.color_code import ColorCode, LatticeType
from mqt.qecc.cc_decoder.comparison import tn_decoder
from mqt.qecc.cc_decoder.hexagonal_color_code import HexagonalColorCode
from mqt.qecc.cc_decoder.square_octagon_color_code import SquareOctagonColorCode


def code_from_string(lattice_type: str, distance: int) -> ColorCode:
    """Construct a color code from a string defining the lattice and a distance."""
    if lattice_type == LatticeType.HEXAGON:
        return HexagonalColorCode(distance)
    if lattice_type == LatticeType.SQUARE_OCTAGON:
        return SquareOctagonColorCode(distance)
    msg = f"Unknown lattice type {lattice_type}. Please choose either hexagon or square_octagon."
    raise ValueError(msg)


__all__ = [
    "tn_decoder",
    "HexagonalColorCode",
    "ColorCode",
    "LatticeType",
    "SquareOctagonColorCode",
    "code_from_string",
]
