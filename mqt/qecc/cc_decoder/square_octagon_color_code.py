import numpy as np


class SquareOctagonColorCode:
    def __init__(self, distance: int) -> None:
        """4.8.8 triangular colour code.

        This class can be used to generate the parity check matrix of 4.8.8 triangular colour code.
        This code has parameters [n, k, d] = [1/2*d**2 + d - 1/2, 1, d].

        Args:
        ----
            distance: Distance of the code to generate. Must be an odd integer.
        """
        self.distance = distance
        self.square_ancilla_qubits: set[tuple[int, int]] = set()
        self.octagon_ancilla_qubits: set[tuple[int, int]] = set()
        self.data_qubits: set[tuple[int, int]] = set()
        self.qubits_to_faces: dict[int, list[int]] = {}
        self.faces_to_qubits: dict[int, list[int]] = {}
        self.add_qubits()
        self.H = np.zeros((len(self.ancilla_qubits), len(self.data_qubits)), dtype=int)
        self.construct_layout()

    def add_qubits(self):
        self.bottom_row_ancillas()
        y = 1
        x_max = self.distance

        while y <= (self.distance + self.distance // 2):
            self.even_data_qubit_row(x_max, y)
            if y == 1:
                x_max -= 1
            else:
                x_max -= 2
            y += 1
            if y <= (self.distance + self.distance // 2):
                self.odd_data_qubit_row(x_max, y)
                y += 1
            if y <= (self.distance + self.distance // 2):
                if y % 2 == 0:
                    self.even_ancilla_qubit_row(x_max, y)
                else:
                    self.odd_ancilla_qubit_row(x_max, y)
                y += 1

        self.ancilla_qubits = self.square_ancilla_qubits.union(self.octagon_ancilla_qubits)

    def even_ancilla_qubit_row(self, x_max, y):
        x = y + 1
        n_qubits = 0
        while True:
            if n_qubits >= x_max:
                break
            self.square_ancilla_qubits.add((x, y))
            n_qubits += 1
            x += 3
            if n_qubits >= x_max:
                break
            self.octagon_ancilla_qubits.add((x, y))
            n_qubits += 1
            x += 3

    def odd_ancilla_qubit_row(self, x_max, y):
        x = y - 2
        n_qubits = 0
        while True:
            if n_qubits >= x_max:
                break
            self.octagon_ancilla_qubits.add((x, y))
            n_qubits += 1
            x += 3
            if n_qubits >= x_max:
                break
            self.square_ancilla_qubits.add((x, y))
            n_qubits += 1
            x += 3

    def even_data_qubit_row(self, x_max, y):
        x = y - 1
        n_qubits = 0
        while True:
            if n_qubits >= x_max:
                break
            self.data_qubits.add((x, y))
            n_qubits += 1
            x += 2
            if n_qubits >= x_max:
                break
            self.data_qubits.add((x, y))
            n_qubits += 1
            x += 4

    def odd_data_qubit_row(self, x_max, y):
        x = y + 1
        n_qubits = 0
        while True:
            if n_qubits >= x_max:
                break
            self.data_qubits.add((x, y))
            n_qubits += 1
            x += 2
            if n_qubits >= x_max:
                break
            self.data_qubits.add((x, y))
            n_qubits += 1
            x += 4

    def bottom_row_ancillas(self):
        for x in range(4, self.distance // 2 * 6, 6):
            self.octagon_ancilla_qubits.add((x, 0))

    def construct_layout(self) -> None:
        coords_to_idx: dict[tuple[int, int], int] = {}
        # builds a map: {(x,y): index} for each qubit with coordinates (x,y)
        # initializes the {qubit_index: [faces]} adjacency list
        for idx, q in enumerate(self.data_qubits):
            coords_to_idx[q] = idx
            self.qubits_to_faces[idx] = []

        for idx, (x, y) in enumerate(self.square_ancilla_qubits):
            qbts: list[int] = []
            for coord in [(x - 1, y - 1), (x + 1, y - 1), (x + 1, y + 1), (x - 1, y + 1)]:
                if coord in coords_to_idx:
                    qubit_idx = coords_to_idx[coord]
                    qbts.append(qubit_idx)
                    self.qubits_to_faces[qubit_idx].append(idx)

            self.faces_to_qubits[idx] = qbts
            for qb in qbts:
                self.H[idx, qb] = 1

        for idx, (x, y) in enumerate(self.octagon_ancilla_qubits):
            idx += len(self.square_ancilla_qubits)
            qbts: list[int] = []
            for coord in [
                (x - 2, y - 1),
                (x - 1, y - 2),
                (x + 1, y - 2),
                (x + 2, y - 1),
                (x + 2, y + 1),
                (x + 1, y + 2),
                (x - 1, y + 2),
                (x - 2, y + 1),
            ]:
                if coord in coords_to_idx:
                    qubit_idx = coords_to_idx[coord]
                    qbts.append(qubit_idx)
                    self.qubits_to_faces[qubit_idx].append(idx)

            self.faces_to_qubits[idx] = qbts
            for qb in qbts:
                self.H[idx, qb] = 1
