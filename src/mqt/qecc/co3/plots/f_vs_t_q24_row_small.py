"""Creates results for circuits with 24 qubits."""

# from layouts import gen_layout, remove_edge_per_factory
# from evaluation import collect_data_space_time, plot_space_time, plot_ratio_vs_t, plot_f_vs_t
from __future__ import annotations

import pickle  # noqa: S403
from pathlib import Path

import mqt.qecc.co3 as co

path = "./results/f_vs_time_q24_ratio08_small_row_250321_2_8_f_t_B"

# ROW
factories_q24_row = [(0, 3), (0, 9), (0, 15), (6, 6), (6, 12), (6, 18), (5, 3), (2, 2)]

g, data_qubit_locs, factory_ring = co.plots.gen_layout("row", 24, factories_q24_row)
custom_layout_q24_row_f8_cc = [data_qubit_locs, g.copy()]

g, data_qubit_locs, factory_ring = co.plots.gen_layout("row", 24, factories_q24_row[:2])
custom_layout_q24_row_f4_cc = [data_qubit_locs, g.copy()]

hc_params = {
    "metric": "crossing",
    "max_restarts": 10,
    "max_iterations": 50,
    "routing": "dynamic",
    "optimize_factories": False,  # True, #want to count in factory crossings as well
    "free_rows": None,
    "parallel": True,
    "processes": 8,
}

instances = [
    # CC GRAPH
    # f=8
    {
        "q": 24,
        "t": 2,
        "min_depth": 24 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q24_row_f8_cc,
        "factory_locs": factories_q24_row,
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "CC",
        "circuit_type": "random",
    },
    {
        "q": 24,
        "t": 8,
        "min_depth": 24 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q24_row_f8_cc,
        "factory_locs": factories_q24_row,
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "CC",
        "circuit_type": "random",
    },
    # f=4
    {
        "q": 24,
        "t": 2,
        "min_depth": 24 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q24_row_f4_cc,
        "factory_locs": factories_q24_row[:2],
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "CC",
        "circuit_type": "random",
    },
    {
        "q": 24,
        "t": 8,
        "min_depth": 24 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q24_row_f4_cc,
        "factory_locs": factories_q24_row[:2],
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "CC",
        "circuit_type": "random",
    },
]


reps = 50
both_metric = True


res_lst = co.plots.collect_data_space_time(instances, hc_params, reps, path, both_metric)
# print(len(custom_layout_q24_row_f8_cc[1].edges()))
# print(len(custom_layout_q24_row_f8_FSC[1].edges()))


with Path(path).open("rb") as f:
    res_lst = pickle.load(f)  # noqa: S301

with Path(path).open("rb") as f:
    res_lst = pickle.load(f)  # noqa: S301

path = "./results/f_vs_time_q24_ratio08_small_row_250321_2_8_f_t_B_metricrouting"

with Path(path).open("rb") as f:
    res_lst_routing = pickle.load(f)  # noqa: S301

for i, res in enumerate(res_lst_routing):
    layout_type = res["instances"][i]["layout_name"]
    q = res["instances"][i]["q"]
    ratio = res["instances"][i]["ratio"]
    factories = len(res["instances"][i]["factory_locs"])
    num_init_list = res["num_init_lst"]
    num_final_list = res["num_final_lst"]
    improvement_lst = []
    for ni, nf in zip(num_init_list, num_final_list):
        improvement_lst.append((ni - nf) / ni)

# check improvements from hc. not too high for small ratios and sparser layouts.
for i, res in enumerate(res_lst):
    layout_type = res["instances"][i]["layout_name"]
    q = res["instances"][i]["q"]
    ratio = res["instances"][i]["ratio"]
    factories = len(res["instances"][i]["factory_locs"])
    num_init_list = res["num_init_lst"]
    num_final_list = res["num_final_lst"]
    improvement_lst = []
    for ni, nf in zip(num_init_list, num_final_list):
        improvement_lst.append((ni - nf) / ni)
