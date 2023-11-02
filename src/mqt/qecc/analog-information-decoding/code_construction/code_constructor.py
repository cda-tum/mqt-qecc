from __future__ import annotations

import json
import os
import subprocess

import ldpc.code_util
import numpy as np
import scipy.io as sio
import scipy.sparse as scs
from bposd.hgp import hgp
from ldpc import mod2
from scipy import sparse


class HD_HGP:
    def __init__(self, boundaries) -> None:
        self.boundaries = boundaries


def is_all_zeros(array):
    return not np.any(array)


def run_checks_scipy(d_1, d_2, d_3, d_4):
    sd_1 = scs.csr_matrix(d_1)
    sd_2 = scs.csr_matrix(d_2)
    sd_3 = scs.csr_matrix(d_3)
    sd_4 = scs.csr_matrix(d_4)

    if not (
        is_all_zeros((sd_1 * sd_2).todense() % 2)
        and is_all_zeros((sd_2 * sd_3).todense() % 2)
        and is_all_zeros((sd_3 * sd_4).todense() % 2)
    ):
        msg = "Error generating 4D code, boundary maps do not square to zero"
        raise Exception(msg)


def generate_4D_product_code(A_1, A_2, A_3, P, checks=True):
    r, c = P.shape
    id_r = np.identity(r, dtype=np.int32)
    id_c = np.identity(c, dtype=np.int32)
    id_n0 = np.identity(A_1.shape[0], dtype=np.int32)
    id_n1 = np.identity(A_2.shape[0], dtype=np.int32)
    id_n2 = np.identity(A_3.shape[0], dtype=np.int32)
    id_n3 = np.identity(A_3.shape[1], dtype=np.int32)

    d_1 = np.hstack((np.kron(A_1, id_r), np.kron(id_n0, P)))

    x = np.hstack((np.kron(A_2, id_r), np.kron(id_n1, P)))
    y = np.kron(A_1, id_c)
    dims = (y.shape[0], x.shape[1] - y.shape[1])
    z = np.hstack((np.zeros(dims, dtype=np.int32), y))
    d_2 = np.vstack((x, z))

    x = np.hstack((np.kron(A_3, id_r), np.kron(id_n2, P)))
    y = np.kron(A_2, id_c)
    dims = (y.shape[0], x.shape[1] - y.shape[1])
    z = np.hstack((np.zeros(dims, dtype=np.int32), y))
    d_3 = np.vstack((x, z))

    d_4 = np.vstack((np.kron(id_n3, P), np.kron(A_3, id_c)))

    if checks:
        run_checks_scipy(d_1, d_2, d_3, d_4)

    return d_1, d_2, d_3, d_4


def generate_3D_product_code(A_1, A_2, P):
    r, c = P.shape
    id_r = np.identity(r, dtype=np.int32)
    id_c = np.identity(c, dtype=np.int32)
    id_n0 = np.identity(A_1.shape[0], dtype=np.int32)
    id_n1 = np.identity(A_2.shape[0], dtype=np.int32)
    id_n2 = np.identity(A_2.shape[1], dtype=np.int32)

    d_1 = np.hstack((np.kron(A_1, id_r), np.kron(id_n0, P)))

    x = np.hstack((np.kron(A_2, id_r), np.kron(id_n1, P)))
    y = np.kron(A_1, id_c)
    dims = (y.shape[0], x.shape[1] - y.shape[1])
    z = np.hstack((np.zeros(dims, dtype=np.int32), y))
    d_2 = np.vstack((x, z))

    d_3 = np.vstack((np.kron(id_n2, P), np.kron(A_2, id_c)))

    if not (is_all_zeros(d_1 @ d_2 % 2) and is_all_zeros(d_2 @ d_3 % 2)):
        msg = "Error generating 3D code, boundary maps do not square to zero"
        raise Exception(msg)

    return d_1, d_2, d_3


def create_outpath(codename: str) -> str:
    path = f"generated_codes/{codename}/"

    if not os.path.exists(path):
        os.makedirs(path)

    return path


def save_code(hx, hz, mx, mz, codename, lx=None, lz=None):
    path = create_outpath(codename)

    Ms = [hx, hz, mx, mz, lx, lz]
    names = ["hx", "hz", "mx", "mz", "lx", "lz"]
    for mat, name in zip(Ms, names):
        if type(mat) != type(None):
            np.savetxt(path + name + ".txt", mat, fmt="%i")
            sio.mmwrite(
                path + name + ".mtx",
                sparse.coo_matrix(mat),
                comment="Field: GF(2)",
            )


def run_compute_distances(codename):
    path = "generated_codes/" + codename
    subprocess.run(["bash", "compute_distances.sh", path])


def _compute_distances(hx, hz, codename):
    run_compute_distances(codename)
    code_dict = {}
    m, n = hx.shape
    codeK = n - mod2.rank(hx) - mod2.rank(hz)
    with open(f"generated_codes/{codename}/info.txt") as f:
        code_dict = dict(
            line[: line.rfind("#")].split(" = ") for line in f if not line.startswith("#") and line.strip()
        )

    code_dict["n"] = n
    code_dict["k"] = codeK
    code_dict["dX"] = int(code_dict["dX"])
    code_dict["dZ"] = int(code_dict["dZ"])
    code_dict["dMX"] = int(code_dict["dMX"])
    code_dict["dMZ"] = int(code_dict["dMZ"])

    print("Code properties:", code_dict)
    with open(f"generated_codes/{codename}/code_params.txt", "w") as file:
        file.write(json.dumps(code_dict))

    return


def _compute_logicals(hx, hz):
    def compute_lz(hx, hz):
        # lz logical operators
        # lz\in ker{hx} AND \notin Im(Hz.T)

        ker_hx = mod2.nullspace(hx)  # compute the kernel basis of hx
        im_hzT = mod2.row_basis(hz)  # compute the image basis of hz.T

        # in the below we row reduce to find vectors in kx that are not in the image of hz.T.
        log_stack = np.vstack([im_hzT, ker_hx])
        pivots = mod2.row_echelon(log_stack.T)[3]
        log_op_indices = [i for i in range(im_hzT.shape[0], log_stack.shape[0]) if i in pivots]
        return log_stack[log_op_indices]

    lx = compute_lz(hz, hx)
    lz = compute_lz(hx, hz)
    return lx, lz


def create_code(
    constructor: str,
    seed_codes: list,
    codename: str,
    compute_distance: bool = False,
    compute_logicals: bool = False,
    lift_parameter=None,
    checks: bool = False,
):
    # Construct initial 2 dim code
    if constructor == "hgp":
        code = hgp(seed_codes[0], seed_codes[1])
    else:
        msg = f"No constructor specified or the specified constructor {constructor} not implemented."
        raise ValueError(msg)

    # Extend to 3D HGP
    A1 = code.hx
    A2 = code.hz.T
    res = generate_3D_product_code(A1, A2, seed_codes[2])

    # Build 4D HGP code
    mx, hx, hzT, mzT = generate_4D_product_code(*res, seed_codes[3], checks=checks)

    hz = hzT.T
    mz = mzT.T

    # Perform checks
    if np.any(hzT @ mzT % 2) or np.any(hx @ hzT % 2) or np.any(mx @ hx % 2):
        msg = "err"
        raise Exception(msg)
    save_code(hx, hz, mx, mz, codename)

    if compute_logicals:
        lx, lz = _compute_logicals(hx, hz)
        save_code(hx, hz, mx, mz, codename, lx=lx, lz=lz)

    else:
        save_code(hx, hz, mx, mz, codename)

    if compute_distance:
        _compute_distances(hx, hz, codename)
    return


if __name__ == "__main__":
    for d in range(3, 8):
        seed_codes = [
            ldpc.codes.ring_code(d),
            ldpc.codes.ring_code(d),
            ldpc.codes.ring_code(d),
            ldpc.codes.ring_code(d),
        ]

        constructor = "hgp"
        codename = f"4D_toric_{d:d}"
        compute_distance = False
        compute_logicals = True
        create_code(
            constructor,
            seed_codes,
            codename,
            compute_distance,
            compute_logicals,
        )
