from __future__ import annotations

import json
import os
import subprocess
from typing import TYPE_CHECKING, Any

import numpy as np
import scipy.io as sio
from bposd.hgp import hgp
from ldpc import mod2
from scipy import sparse
from scipy.sparse import coo_matrix, csr_matrix

from mqt.qecc.analog_information_decoding.code_construction import code_constructor

if TYPE_CHECKING:
    from numpy._typing import NDArray


def is_all_zeros(array: NDArray[np.int_]) -> bool:
    """Check if array is all zeros."""
    return not np.any(array)


def sparse_all_zeros(mat: csr_matrix) -> bool:
    """Check if sparse matrix is all zeros."""
    mat.data %= 2
    return bool(mat.sum() == 0)


def run_checks_scipy(d_1: csr_matrix, d_2: csr_matrix, d_3: csr_matrix, d_4: csr_matrix) -> None:
    if not (sparse_all_zeros(d_1 @ d_2) and sparse_all_zeros(d_2 @ d_3) and sparse_all_zeros(d_3 @ d_4)):
        msg = "Error generating 4D code, boundary maps do not square to zero"
        raise Exception(msg)


def generate_4D_product_code(
    a_1: csr_matrix,
    a_2: csr_matrix,
    a_3: csr_matrix,
    p: csr_matrix,
    checks: bool = True,
) -> tuple[csr_matrix, csr_matrix, csr_matrix, csr_matrix]:
    r, c = p.shape

    id_r = sparse.identity(r, dtype=int)
    id_c = sparse.identity(c, dtype=int)
    id_n0 = sparse.identity(a_1.shape[0], dtype=int)
    id_n1 = sparse.identity(a_2.shape[0], dtype=int)
    id_n2 = sparse.identity(a_3.shape[0], dtype=int)
    id_n3 = sparse.identity(a_3.shape[1], dtype=int)

    d_1 = sparse.hstack((sparse.kron(a_1, id_r), sparse.kron(id_n0, p)))

    x = sparse.hstack((sparse.kron(a_2, id_r), sparse.kron(id_n1, p)))
    y = sparse.kron(a_1, id_c)
    dims = (y.shape[0], x.shape[1] - y.shape[1])
    nmat = csr_matrix(np.zeros(dims))
    z = sparse.hstack((nmat, y))
    d_2 = sparse.vstack((x, z))

    x = sparse.hstack((sparse.kron(a_3, id_r), sparse.kron(id_n2, p)))
    y = sparse.kron(a_2, id_c)
    dims = (y.shape[0], x.shape[1] - y.shape[1])
    mat = csr_matrix(np.zeros(dims))
    z = sparse.hstack([mat, y])
    d_3 = sparse.vstack((x, z))

    d_4 = sparse.vstack((sparse.kron(id_n3, p), sparse.kron(a_3, id_c)))

    if checks:
        run_checks_scipy(d_1, d_2, d_3, d_4)

    return d_1, d_2, d_3, d_4


def generate_3D_product_code(
    a_1: csr_matrix, a_2: csr_matrix, p: csr_matrix
) -> tuple[csr_matrix, csr_matrix, csr_matrix]:
    r, c = p.shape

    id_r = sparse.identity(r, dtype=int)
    id_c = sparse.identity(c, dtype=int)
    id_n0 = sparse.identity(a_1.shape[0], dtype=int)
    id_n1 = sparse.identity(a_2.shape[0], dtype=int)
    id_n2 = sparse.identity(a_2.shape[1], dtype=int)

    d_1 = sparse.hstack((sparse.kron(a_1, id_r), sparse.kron(id_n0, p)))

    x = sparse.hstack((sparse.kron(a_2, id_r), sparse.kron(id_n1, p)))
    y = sparse.kron(a_1, id_c)
    dims = (y.shape[0], x.shape[1] - y.shape[1])
    z = sparse.hstack((csr_matrix(np.zeros(dims), dtype=int), y))
    d_2 = sparse.vstack((x, z))

    d_3 = sparse.vstack((sparse.kron(id_n2, p), sparse.kron(a_2, id_c)))

    if not (sparse_all_zeros(d_1 @ d_2) and sparse_all_zeros(d_2 @ d_3)):
        msg = "Error generating 3D code, boundary maps do not square to zero"
        raise Exception(msg)

    return d_1, d_2, d_3  # mx, hx, hzT # hx, hzT, mzT


def create_outpath(codename: str) -> str:
    """Create output path for code."""
    path = f"/codes/generated_codes/{codename}/"

    if not os.path.exists(path):
        os.makedirs(path)

    return path


def save_code(
    hx: csr_matrix,
    hz: csr_matrix,
    mz: csr_matrix,
    codename: str,
    lx: csr_matrix = None,
    lz: csr_matrix = None,
) -> None:
    path = create_outpath(codename)

    Ms: list[csr_matrix] = [hx, hz, mz, lx, lz]
    names = ["hx", "hz", "mz", "lx", "lz"]
    for mat, name in zip(Ms, names):
        if type(mat) != type(None):
            path_str = path + name
            try:
                np.savetxt(path_str + ".txt", mat.todense(), fmt="%i")
            except:
                np.savetxt(path_str + ".txt", mat, fmt="%i")
            sio.mmwrite(
                path_str + ".mtx",
                coo_matrix(mat),
                comment="Field: GF(2)",
                field="integer",
            )


def run_compute_distances(codename: str) -> None:
    path = "/codes/generated_codes/" + codename
    subprocess.run(["bash", "compute_distances_3D.sh", path])


def _compute_distances(hx: NDArray[np.int32], hz: NDArray[np.int32], codename: str) -> None:
    run_compute_distances(codename)
    _, n = hx.shape
    codeK = n - mod2.rank(hx) - mod2.rank(hz)
    with open(f"/codes/generated_codes/{codename}/info.txt") as f:
        code_dict: dict[str, Any] = dict(
            line[: line.rfind("#")].split(" = ") for line in f if not line.startswith("#") and line.strip()
        )

    code_dict["n"] = n
    code_dict["k"] = codeK
    code_dict["dX"] = int(code_dict["dX"])
    code_dict["dZ"] = int(code_dict["dZ"])

    print("Code properties:", code_dict)
    with open(f"/codes/generated_codes/{codename}/code_params.txt", "w") as file:
        file.write(json.dumps(code_dict))


def _store_code_params(hx: csr_matrix, hz: csr_matrix, codename: str) -> None:
    code_dict = {}
    hx, hz = hx.todense(), hz.todense()
    m, n = hx.shape
    codeK = n - mod2.rank(hx) - mod2.rank(hz)
    code_dict["n"] = n
    code_dict["k"] = codeK
    with open(f"/codes/generated_codes/{codename}/code_params.txt", "w") as file:
        file.write(json.dumps(code_dict))


def create_code(
    constructor: str,
    seed_codes: list[csr_matrix],
    codename: str,
    compute_distance: bool = False,
    compute_logicals: bool = False,
    checks: bool = False,
) -> None:
    # Construct initial 2 dim code
    if constructor == "hgp":
        code = hgp(seed_codes[0], seed_codes[1])
    else:
        msg = f"No constructor specified or the specified constructor {constructor} not implemented."
        raise ValueError(msg)

    # Extend to 3D HGP
    A1 = sparse.csr_matrix(code.hx)
    A2 = sparse.csr_matrix(code.hz.T)
    res = generate_3D_product_code(A1, A2, sparse.csr_matrix(seed_codes[2]))
    hx, hzT, mzT = res

    hz = hzT.transpose()
    mz = mzT.transpose()
    if compute_logicals:
        lx, lz = code_constructor._compute_logicals(hx.todense(), hz.todense())
        save_code(hx=hx, hz=hz, mz=mz, codename=codename, lx=lx, lz=lz)

    if compute_distance:
        _compute_distances(hx.todense(), hz.todense(), codename)

    else:
        _store_code_params(hx.todense(), hz.todense(), codename)
    return
