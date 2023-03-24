from __future__ import annotations

import datetime
import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from mqt.qecc.cc_decoder.hexagonal_color_code import HexagonalColorCode
from z3 import Bool, Not, Optimize, Xor, simplify

if TYPE_CHECKING:
    from z3 import ModelRef


@dataclass
class LightsOut:
    lights_to_switches: dict[int, list[int]]
    switches_to_lights: dict[int, list[int]]
    switch_vars: list[Bool] | None = None
    helper_vars: dict[int, list[Bool]] = field(default_factory=dict)
    optimizer: Optimize = field(default_factory=Optimize)

    def preconstruct_parity_constraint(self, light: int, indices: list[int]) -> None:
        """
        Preconstruct the parity constraints for a light.
        Adds all constraint to the optimizer that are independent of the value of the light.
        """
        helper_vars = self.helper_vars[light]

        for i in range(1, len(indices) - 1):
            constraint = Xor(self.switch_vars[indices[i]], helper_vars[i]) == helper_vars[i - 1]
            self.optimizer.add(simplify(constraint))

        constraint = self.switch_vars[indices[-1]] == helper_vars[-1]
        self.optimizer.add(simplify(constraint))

    def complete_parity_constraint(self, light: int, indices: list[int], val: bool) -> None:
        """
        Completes the parity constraints for a light.
        Adds the constraint that is dependent on the value of the light.
        """
        helper_vars = self.helper_vars[light]

        constraint = Xor(self.switch_vars[indices[0]], helper_vars[0]) == val
        self.optimizer.add(simplify(constraint))

    def preconstruct_z3_instance(self) -> None:
        """
        Preconstruct the z3 instance for the lights-out problem
        Creates all necessary variables, adds the known parts of the parity constraints.
        Soft constraints are added to the optimizer with default weights.
        """
        if self.switch_vars is None:
            self.switch_vars = [Bool(f"switch_{i}") for i in range(len(self.switches_to_lights))]

        for light, switches in self.lights_to_switches.items():
            if light not in self.helper_vars:
                self.helper_vars[light] = [Bool(f"helper_{light}_{i}") for i in range(len(switches) - 1)]
            self.preconstruct_parity_constraint(light, switches)

        for switch in self.switch_vars:
            self.optimizer.add_soft(Not(switch))

    def validate_model(self, model: ModelRef, lights: list[bool]) -> bool:
        """
        Validate the model by checking if pressing the switches turns off all lights.
        """
        for i, var in enumerate(self.switch_vars):
            if model[var]:
                # flip all lights that are controlled by this switch
                for light in self.switches_to_lights[i]:
                    lights[light] = not lights[light]

        return all(not light for light in lights)

    def count_switches(self, model: ModelRef) -> int:
        """
        Count the number of switches that are set to true.
        """
        return sum(1 for var in self.switch_vars if model[var])

    def solve(
        self, lights: list[bool], solver_path: str = "z3"
    ) -> tuple[list[int], datetime.timedelta, datetime.timedelta]:
        """
        Solve the lights-out problem for a given pattern.
        Assumes that the z3 instance has already been pre-constructed.
        """
        # push a new context to the optimizer
        self.optimizer.push()

        # add the problem specific constraints
        start = datetime.datetime.now()
        for light, val in enumerate(lights):
            self.complete_parity_constraint(light, self.lights_to_switches[light], val)
        constr_time = datetime.datetime.now() - start

        if solver_path == "z3":
            # solve the problem
            start = datetime.datetime.now()
            result = self.optimizer.check()
            solve_time = datetime.datetime.now() - start
            assert str(result) == "sat", "No solution found"

            # validate the model
            model = self.optimizer.model()
            assert self.validate_model(model, lights), "Model is invalid"

            switches = [1 if model[var] else 0 for var in self.switch_vars]
        else:
            switches = []
            self.optimizer.set("pp.wcnf", True)
            wcnf = str(self.optimizer)
            # Note: This merely calls the solver. It does not interpret the output.
            #       This is just to measure the time it takes to solve the problem.
            with Path("./solver-out_" + solver_path.split("/")[-1] + ".txt").open("a+") as out:
                start = datetime.datetime.now()
                subprocess.run([solver_path, wcnf], stdout=out)
                solve_time = datetime.datetime.now() - start

        # pop the context from the optimizer
        self.optimizer.pop()

        return switches, constr_time, solve_time


def simulate_error_rate(distance: int, error_rate: float, nr_sims: int, solver_path: str = "z3") -> dict[str, Any]:
    code = HexagonalColorCode(distance=distance)
    problem = LightsOut(code.faces_to_qubits, code.qubits_to_faces)

    start = datetime.datetime.now()
    problem.preconstruct_z3_instance()
    preconstr_time = datetime.datetime.now() - start
    min_wt_logical = 0
    logical_errors = 0
    avg_constr_time = 0.0
    avg_solve_time = 0.0
    rng = np.random.default_rng()
    for i in range(nr_sims):
        # sample error
        error = rng.choice([0, 1], size=code.n, p=[1 - error_rate, error_rate])

        # get syndrome
        syndrome = code.get_syndrome(error)
        lights = [bool(b == 1) for b in syndrome]  # noqa: SIM901

        # compute estimate
        estimate, constr_time, solve_time = problem.solve(lights, solver_path=solver_path)
        if len(estimate) > 0:
            # check if the estimate is correct
            residual = (error + np.array(estimate)) % 2
            if code.check_if_logical_error(residual):
                logical_errors += 1
                wt = np.sum(residual)  # compute the min weight of a logical error
                if min_wt_logical == 0 or wt < min_wt_logical:
                    min_wt_logical = int(wt)

        # compute rolling average of the times
        avg_constr_time = (avg_constr_time * i + constr_time.microseconds) / (i + 1)
        avg_solve_time = (avg_solve_time * i + solve_time.microseconds) / (i + 1)

    avg_total_time = avg_constr_time + avg_solve_time
    logical_error_rate = logical_errors / nr_sims
    logical_error_rate_eb = np.sqrt((1 - logical_error_rate) * logical_error_rate / nr_sims)

    return {
        "distance": distance,
        "p": error_rate,
        "logical_error_rate": logical_error_rate,
        "logical_error_rate_eb": logical_error_rate_eb,
        "preconstr_time": preconstr_time.microseconds,
        "avg_constr_time": avg_constr_time,
        "avg_solve_time": avg_solve_time,
        "avg_total_time": avg_total_time,
        "min_wt_logical_err": min_wt_logical,
    }


def run(
    distance: int, error_rate: float, nr_sims: int = 10000, results_dir: str = "./results_maxsat", solver: str = "z3"
) -> None:
    data = simulate_error_rate(distance, error_rate, nr_sims, solver)
    str = solver.split("/")[-1]
    filename = f"./distance={distance},p={round(error_rate, 4)},solver={str}.json"
    path = Path(results_dir)
    path.mkdir(parents=True, exist_ok=True)
    with (path / filename).open("w") as out:
        out.write(json.dumps(data))
