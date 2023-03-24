from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from scipy.optimize import curve_fit


def plot_ler_vs_distance(code_dict: dict[float, Any], ax: Axes, pers: list[float]) -> None:
    for p in pers:
        ds = []
        lers = []
        for d, ler in sorted(code_dict[p].items()):
            ds.append(int(d))
            lers.append(ler)
        ax.plot(ds, lers, label="p=" + str(round(p, 2)), marker="o", linestyle="-")
    ax.set_yscale("log")
    ax.legend()
    ax.set_ylabel("Logical failure rate")
    ax.set_xlabel(r"Code distance $\it{d}$")


def threshold_fit(variables: tuple[float, float], B0: float, B1: float, B2: float, mu: float, pth: float) -> None:
    p, L = variables
    return B0 + B1 * (p - pth) * pow(L, 1 / mu) + B2 * pow((p - pth) * pow(L, 1 / mu), 2)


def calculate_threshold(
    code_dict: dict[int, Any],
    min_distance: int = 1,
    max_distance: int = 100,
    distances: list[int] = None,
    min_per: float = 0.06,
    max_per: float = 0.13,
    title: str = None,
    ax: Axes = None,
) -> None:
    ler_data = []
    ler_eb_data = []
    per_data = []
    distance_data = []
    for distance in code_dict:
        per_array = []
        ler_array = []
        ler_eb = []
        if distances is None:
            if not min_distance < int(distance) < max_distance:
                continue
        elif int(distance) not in distances:
            continue

        for index, per in enumerate(code_dict[distance]["p"]):
            if min_per < per < max_per:
                per_array.append(per)
                ler_array.append(code_dict[distance]["logical_error_rate"][index])
                ler_eb.append(code_dict[distance]["logical_error_rate_eb"][index])
        per_data.extend(per_array)
        ler_data.extend(ler_array)
        ler_eb_data.extend(ler_eb)
        distance_data.extend([int(distance) for _ in range(len(per_array))])

    popt, pcov = curve_fit(threshold_fit, (per_data, distance_data), ler_data, maxfev=10000)
    if ax is not None:
        ax.axvline(x=popt[-1], color="black", linestyle="dashed")
        print("threshold: ", popt[-1])
    error = np.sqrt(np.diag(pcov))[-1]
    print("error: " + str(error))

    distance_array = []
    for distance in code_dict:
        if distances is None:
            if min_distance < int(distance) < max_distance:
                distance_array.append(int(distance))
        elif int(distance) in distances:
            distance_array.append(int(distance))
    distance_array.sort()
    for distance in distance_array:
        per_array = code_dict[distance]["p"]
        ler_array = code_dict[distance]["logical_error_rate"]
        ler_eb = code_dict[distance]["logical_error_rate_eb"]

        if per_array != [] and ax is not None:
            ax.errorbar(per_array, ler_array, yerr=ler_eb, label="d = " + str(distance), fmt="|")

    if ax is not None:
        ax.legend()
        ax.set_xlabel("Physical error rate")
        ax.set_ylabel("Logical error rate")
        ax.set_title(title)
        ax.set_yscale("log")
        ax.set_xlim(min_per, max_per)


def generate_plots(results_dir: Path, results_file: Path) -> None:
    """
    Generates the plots for the paper
    """
    # read in all generated data
    data = []
    for file in results_dir.glob("*.json"):
        with file.open() as f:
            data.append(json.loads(f.read()))

    fig, ax = plt.subplots(3, 2, figsize=(12, 12))
    metrics: dict[int, Any] = {}
    per_metrics: dict[float, Any] = {}

    for result in data:
        d = result["distance"]
        p = result["p"]

        if d not in metrics:
            metrics[d] = {
                "p": [],
                "logical_error_rate": [],
                "logical_error_rate_eb": [],
                "avg_total_time": [],
                "min_wt_logical_err": -1,
            }
        if p not in per_metrics:
            per_metrics[p] = {}

        metrics[d]["p"].append(p)
        metrics[d]["logical_error_rate"].append(result["logical_error_rate"])
        metrics[d]["logical_error_rate_eb"].append(result["logical_error_rate_eb"])
        metrics[d]["avg_total_time"].append(result["avg_total_time"])

        if result["min_wt_logical_err"] > 0:
            if metrics[d]["min_wt_logical_err"] == -1:
                metrics[d]["min_wt_logical_err"] = result["min_wt_logical_err"]
            else:
                metrics[d]["min_wt_logical_err"] = min(metrics[d]["min_wt_logical_err"], result["min_wt_logical_err"])
        per_metrics[p][d] = result["logical_error_rate"]
    for d, mdata in sorted(metrics.items()):
        (
            mdata["p"],
            mdata["logical_error_rate"],
            mdata["avg_total_time"],
            mdata["logical_error_rate_eb"],
        ) = zip(
            *sorted(
                zip(mdata["p"], mdata["logical_error_rate"], mdata["avg_total_time"], mdata["logical_error_rate_eb"])
            )
        )
        ax[0][0].errorbar(mdata["p"], mdata["logical_error_rate"], mdata["logical_error_rate_eb"], label=f"d={d}")
        ax[0][0].set_xlabel("Physical error rate")
        ax[0][0].set_ylabel("Logical error rate")
        ax[0][0].legend()

        ax[1][0].plot(mdata["p"], mdata["avg_total_time"], label="d=" + str(d))
        ax[1][0].set_xlabel("Physical error rate")
        ax[1][0].set_ylabel("Average time per run (µs)")
        ax[1][0].legend()
        ax[1][0].set_ylim(0, 185000)

    # plot the average time per run over the distance
    ds = []
    p_data: dict[float, Any] = {}
    pers = [0.001, 0.02, 0.05, 0.08, 0.13]
    for d, mdata in sorted(metrics.items()):
        ds.append(d)
        for i, p in enumerate(mdata["p"]):
            if p in pers:
                if p not in p_data:
                    p_data[p] = {"d": [], "t": []}
                p_data[p]["d"].append(d)
                p_data[p]["t"].append(mdata["avg_total_time"][i])
    for p, pdata in sorted(p_data.items()):
        ax[1][1].plot(ds, pdata["t"], label="p=" + str(p))

    ax[1][1].set_xlabel("Distance")
    ax[1][1].set_ylabel("Average time per run (µs)")
    ax[1][1].legend()
    ax[1][1].set_xticks(ds)
    ax[1][1].set_ylim(0, 185000)
    calculate_threshold(metrics, ax=ax[0][1], title="Threshold")
    plot_ler_vs_distance(per_metrics, ax=ax[2][0], pers=[0.02, 0.05, 0.08, 0.11, 0.13])
    # save plot as vector graphic
    plt.savefig(results_file, bbox_inches="tight")


def generate_plots_tn(results_dir: Path, results_file: Path) -> None:
    # read in all generated data
    data = []
    for file in results_dir.glob("*.json"):
        with file.open() as f:
            data.append(json.loads(f.read()))

    # prepare code to per,ler map and print
    code_to_xys: dict[float, Any] = {}
    for run in data:
        xys = code_to_xys.setdefault(run["n_k_d"][-1], [])
        xys.append((run["physical_error_rate"], run["logical_failure_rate"]))

    for xys in code_to_xys.values():
        xys.sort(key=lambda xy: xy[0])

    fig, ax = plt.subplots(2, 2, figsize=(12, 10))
    # add data
    for code, xys in sorted(code_to_xys.items()):
        ax[0][0].plot(*zip(*xys), "x-", label=f"d={code}")
    ax[0][0].set_xlabel("Physical error rate")
    ax[0][0].set_ylabel("Logical error rate")
    ax[0][0].legend()

    # prepare code to per,average-time-per-run map and print
    code_to_xys = {}
    for run in data:
        xys = code_to_xys.setdefault(run["n_k_d"][-1], [])
        # convert from seconds to microseconds
        xys.append((run["error_probability"], (run["wall_time"] / run["n_run"]) * 1e6))

    for xys in code_to_xys.values():
        xys.sort(key=lambda xy: xy[0])

    for code, xys in sorted(code_to_xys.items()):
        ax[1][0].plot(*zip(*xys), "x-", label=f"d={code}")
    ax[1][0].set_xlabel("Physical error rate")
    ax[1][0].set_ylabel("Average time per run (µs)")
    ax[1][0].legend()
    ax[1][0].set_ylim(0, 300000)

    ds = []
    p_data: dict[float, Any] = {}
    pers = [0.001, 0.021, 0.051, 0.081, 0.111]
    for d, cdata in sorted(code_to_xys.items()):
        ds.append(d)
        for _, (p, t) in enumerate(cdata):
            if p in pers:
                if p not in p_data:
                    p_data[p] = {"d": [], "t": []}
                p_data[p]["d"].append(d)
                p_data[p]["t"].append(t)
    for p, pdata in sorted(p_data.items()):
        ax[1][1].plot(ds, pdata["t"], label="p=" + str(p))

    ax[1][1].set_xlabel("Distance")
    ax[1][1].set_ylabel("Average time per run (µs)")
    ax[1][1].legend()
    # ax[1][1].set_yscale("log")
    ax[1][1].set_xticks(ds)
    ax[1][1].set_ylim(0, 300000)
    # save plot as vector graphic
    plt.savefig(results_file, bbox_inches="tight")


def generate_plots_comp(results_dir: Path, results_file: Path) -> None:
    fig, ax = plt.subplots(2, figsize=(12, 12))
    cols = ["blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "cyan", "olive"]
    idx = 0
    solverToCol = {}
    pToCol = {}
    for subdir, _, files in os.walk(results_dir):
        if not files:
            continue
        data = []
        solver = subdir.split("/")[-1]
        if solver not in solverToCol:
            solverToCol[solver] = cols[idx]
            idx += 1
        for f in files:
            fp = subdir + "/" + f
            with Path(fp).open() as ff:
                data.append(json.loads(ff.read()))

        metrics: dict[int, dict[str, list[Any]]] = {}
        for result in data:
            d = result["distance"]
            p = result["p"]
            if d not in metrics:
                metrics[d] = {"p": [], "avg_total_time": []}

            metrics[d]["p"].append(p)
            metrics[d]["avg_total_time"].append(result["avg_total_time"])

        for d, mdata in sorted(metrics.items()):
            if d == 21 or d == 17:
                (
                    mdata["p"],
                    mdata["avg_total_time"],
                ) = zip(*sorted(zip(mdata["p"], mdata["avg_total_time"])))
                ax[0].plot(
                    mdata["p"],
                    mdata["avg_total_time"],
                    label="solver=" + solver + ", d=" + str(d),
                    color=solverToCol[solver],
                )
        ds = []
        idx = 0
        p_data: dict[float, dict[str, list[Any]]] = {}
        pers = [0.001, 0.051, 0.131]
        for d, pdata in sorted(metrics.items()):
            ds.append(d)
            for i, p in enumerate(pdata["p"]):
                if p in pers:
                    if p not in p_data:
                        p_data[p] = {"d": [], "t": []}
                        pToCol[p] = cols[idx]
                        idx += 1
                    p_data[p]["d"].append(d)
                    p_data[p]["t"].append(pdata["avg_total_time"][i])
        for p, ppdata in sorted(p_data.items()):
            ax[1].plot(
                ds,
                ppdata["t"],
                label="p=" + str(p) + ", " + (solver if solver == "z3" else "CASHWMaxSAT‑CorePlus"),
                color=pToCol[p],
                marker="x" if solver == "z3" else "o",
            )

        ax[1].set_xlabel("Distance")
        ax[1].set_ylabel("Average time per run (µs)")
        ax[1].legend()
        ax[1].set_xticks(ds)

        ax[0].set_xlabel("Physical error rate")
        ax[0].set_ylabel("Average time per run (µs)")
        ax[0].legend()
        ax[0].legend(loc="center left", bbox_to_anchor=(1, 0.5))

    # save plot as vector graphic
    plt.savefig(results_file, bbox_inches="tight")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "results_dir",
        type=str,
        help="Path to the directory containing the results",
    )
    parser.add_argument(
        "--mode",
        type=str,
        default="maxsat",
        help="Mode to use. Either 'maxsat' or 'tn'. Default: maxsat",
    )
    parser.add_argument(
        "--results_file",
        type=str,
        default="./results.pdf",
        help="File where the results should be saved. Default: ./results.pdf",
    )
    args = parser.parse_args()

    if args.mode == "maxsat":
        generate_plots(Path(args.results_dir), Path(args.results_file))
    elif args.mode == "tn":
        generate_plots_tn(Path(args.results_dir), Path(args.results_file))
    elif args.mode == "comparison":
        generate_plots_comp(Path(args.results_dir), Path(args.results_file))
    else:
        msg = "Unknown mode: " + args.mode
        raise ValueError(msg)
