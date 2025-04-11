"""Estimate logical error rate for CSS state preparation circuits for a given code and physical error rate."""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

from qiskit import QuantumCircuit

from mqt.qecc import CSSCode
from mqt.qecc.circuit_synthesis.simulation import SteaneNDFTStatePrepSimulator
from mqt.qecc.circuit_synthesis.state_prep import heuristic_prep_circuit
from mqt.qecc.codes import HexagonalColorCode, SquareOctagonColorCode


def main() -> None:
    """Run the logical error rate estimation for a given code and physical error rate."""
    available_codes = ["steane", "tetrahedral", "shor", "surface", "cc_4_8_8", "cc_6_6_6", "hamming", "carbon"]
    parser = argparse.ArgumentParser(description="Estimate logical error rate for CSS state preparation circuits")
    parser.add_argument(
        "code",
        type=str,
        help="Code for which to estimate logical error rate. Available codes: " + ", ".join(available_codes),
    )
    parser.add_argument("-p", "--p_error", type=float, help="Physical error rate")
    parser.add_argument("-p_idle_factor", "--p_idle_factor", type=float, default=0.01, help="Idling error rate")
    parser.add_argument("--zero_state", default=True, action="store_true", help="Synthesize logical |0> state.")
    parser.add_argument(
        "--plus_state", default=False, dest="zero_state", action="store_false", help="Synthesize logical |+> state."
    )
    parser.add_argument("--x_errors", default=True, action="store_true", help="Calculate error rates for X-errors")
    parser.add_argument(
        "--z_errors", default=False, dest="x_errors", action="store_false", help="Calculate error rates for Z errors"
    )
    parser.add_argument("-n", "--n_errors", type=int, default=500, help="Number of errors to sample")
    parser.add_argument(
        "-d", "--distance", type=int, default=3, help="Code Distance (only required for surface and color codes)"
    )

    args = parser.parse_args()
    code_name = args.code
    if "surface" in code_name:
        d = args.distance
        code = CSSCode.from_code_name("surface", d)
        code_name = f"rotated_surface_d{d}"
    elif code_name == "cc_4_8_8_d7":
        d = 7
        code = SquareOctagonColorCode(d)
        lut_path = (Path("__file__") / "../../eval/luts/decoder_488_7.pickle").resolve()
        if lut_path.exists():
            with lut_path.open("rb") as f:
                lut = pickle.load(f)
        else:
            msg = "LUT file not found."
            raise ValueError(msg)
    elif code_name == "cc_6_6_6_d7":
        d = 7
        code = HexagonalColorCode(d)
    elif code_name == "cc_4_8_8_d5":
        d = 5
        code = SquareOctagonColorCode(d)
    elif code_name == "cc_6_6_6_d5":
        d = 5
        code = HexagonalColorCode(d)
    elif code_name in available_codes:
        code = CSSCode.from_code_name(code_name)
    else:
        raise ValueError("Code " + code_name + " not available. Available codes: " + ", ".join(available_codes))

    prefix = (Path(__file__) / "../circuits/").resolve()
    circ_file_core = f"{code_name}_heuristic_"

    # check if file exists
    # if not (prefix / code_name / circ_file).exists():
    #     # create circuit
    #     # NOTE: error message for missing circuits
    #     pass
    # else:
    circuits = []
    # load circuit from file
    for _id in [0, 1, 2, 3]:
        circ_file = circ_file_core + str(_id)
        circuits.append(QuantumCircuit.from_qasm_file(prefix / code_name / circ_file))

    sim = SteaneNDFTStatePrepSimulator(
        circ1=circuits[0],
        circ2=circuits[1],
        code=code,
        circ3=circuits[2],
        circ4=circuits[3],
        check_circuit=None if args.x_errors else heuristic_prep_circuit(code, zero_state=False).circ,
        p=args.p_error,
        p_idle=args.p_idle_factor * args.p_error,
        decoder=lut if code_name == "cc_4_8_8_d7" else None,
    )
    if args.x_errors:
        res = sim.logical_error_rate(min_errors=args.n_errors)
    else:
        sim.set_p(args.p_error, args.p_idle_factor * args.p_error)
        res = sim.secondary_logical_error_rate(min_errors=args.n_errors)

    print(",".join([str(x) for x in res]))


if __name__ == "__main__":
    main()
