# Created by Peter-Jan Derks
# Hexagonal Color Code layout construction adapted from https://github.com/peter-janderks/restriction_decoder_domain_wall_colour_code
from __future__ import annotations

import numpy as np
from ldpc import mod2


class HexagonalColorCode:
    def __init__(self, distance: int):
        self.distance = distance
        self.ancilla_qubits: set[tuple[int, int]] = set()
        self.data_qubits: set[tuple[int, int]] = set()

        colour = ["r", "b", "g"]
        y = 0
        x_max = distance + distance // 2
        while x_max > 0:
            ancilla_colour = colour[y % 3]
            if ancilla_colour == "r":
                self.red_row(x_max, y)
            elif ancilla_colour == "b":
                self.blue_row(x_max, y)
            else:
                self.green_row(x_max, y)
            x_max -= 1
            y += 1

        self.qubits_to_faces: dict[int, list[int]] = {}
        self.faces_to_qubits: dict[int, list[int]] = {}

        self.H = np.zeros((len(self.ancilla_qubits), len(self.data_qubits)), dtype=int)
        self.L = np.zeros((len(self.ancilla_qubits), len(self.data_qubits)), dtype=int)

        self.construct_layout()
        self.n = len(self.qubits_to_faces)

    def red_row(self, x_max: int, y: int) -> None:
        """
        create red ancilla qubits
        """
        i = 0
        x_row = y
        while i < x_max:
            data_or_ancilla = i % 3
            if data_or_ancilla == 0 or data_or_ancilla == 2:
                self.data_qubits.add((x_row, y))
            else:
                self.ancilla_qubits.add((x_row, y))
            i += 1
            x_row += 2

    def blue_row(self, x_max: int, y: int) -> None:
        """
        create blue ancilla qubits
        """
        i = 0
        x_row = y
        while i < x_max:
            data_or_ancilla = i % 3
            if data_or_ancilla == 0 or data_or_ancilla == 1:
                self.data_qubits.add((x_row, y))
            else:
                self.ancilla_qubits.add((x_row, y))
            i += 1
            x_row += 2

    def green_row(self, x_max: int, y: int) -> None:
        """
        create green ancilla qubits
        """
        i = 0
        x_row = y
        while i < x_max:
            data_or_ancilla = i % 3
            if data_or_ancilla == 1 or data_or_ancilla == 2:
                self.data_qubits.add((x_row, y))
            else:
                self.ancilla_qubits.add((x_row, y))
            x_row += 2
            i += 1

    def compute_logical(self) -> None:
        # lz logical operators
        # lz\in ker{hx} AND \notin Im(Hz.T)
        ker_hx = mod2.nullspace(self.H)  # compute the kernel basis of hx
        im_hz_transp = mod2.row_basis(self.H)  # compute the image basis of hz.T
        log_stack = np.vstack([im_hz_transp, ker_hx])
        pivots = mod2.row_echelon(log_stack.T)[3]
        log_op_indices = [i for i in range(im_hz_transp.shape[0], log_stack.shape[0]) if i in pivots]
        self.L = log_stack[log_op_indices]

    def construct_layout(self) -> None:
        coords_to_idx: dict[tuple[int, int], int] = {}
        # builds a map: {(x,y): index} for each qubit with coordinates (x,y)
        # initializes the {qubit_index: [faces]} adjacency list
        for idx, q in enumerate(self.data_qubits):
            coords_to_idx[q] = idx
            self.qubits_to_faces[idx] = []

        for idx, (x, y) in enumerate(self.ancilla_qubits):
            qbts: list[int] = []
            for coord in [(x - 1, y + 1), (x - 2, y), (x - 1, y - 1), (x + 1, y - 1), (x + 2, y), (x + 1, y + 1)]:
                if coord in coords_to_idx:
                    qubit_idx = coords_to_idx[coord]
                    qbts.append(qubit_idx)
                    self.qubits_to_faces[qubit_idx].append(idx)
            self.faces_to_qubits[idx] = qbts
            for qb in qbts:
                self.H[idx, qb] = 1

        # L is the matrix of logicals of the code
        self.compute_logical()

    def get_syndrome(self, error: np.ndarray) -> np.ndarray:
        return self.H @ error % 2

    def check_if_logical_error(self, residual: np.ndarray) -> bool:
        return (self.L @ residual % 2).any()
