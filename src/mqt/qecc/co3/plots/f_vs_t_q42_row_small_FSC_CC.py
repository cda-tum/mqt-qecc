from __future__ import annotations

import pickle
from pathlib import Path

import mqt.qecc.co3 as co

path = "./results/f_vs_time_q42_ratio08_small_row_FSC"

# ROW
factories_q42_row = [
    (0, 3),
    (0, 9),
    (6, 12),
    (6, 18),
    (0, 15),
    (6, 6),
    (5, 3),
    (2, 2),
    (0, 21),
    (0, 27),
    (6, 24),
    (6, 30),
]  # , (1,33)]

g, data_qubit_locs, factory_ring = co.plots.gen_layout("row", 42, factories_q42_row)
custom_layout_q42_row_f12_CC = [data_qubit_locs, g.copy()]
g = co.plots.remove_edge_per_factory(g, factories_q42_row)
custom_layout_q42_row_f12_FSC = [data_qubit_locs, g]

g, data_qubit_locs, factory_ring = co.plots.gen_layout("row", 42, factories_q42_row[:4])
custom_layout_q42_row_f4_CC = [data_qubit_locs, g.copy()]
g = co.plots.remove_edge_per_factory(g, factories_q42_row[:4])
custom_layout_q42_row_f4_FSC = [data_qubit_locs, g]


hc_params = {
    "metric": "crossing",
    "max_restarts": 30,
    "max_iterations": 50,
    "routing": "dynamic",
    "optimize_factories": False,  # True, #want to count in factory crossings as well
    "free_rows": None,
    "parallel": True,
    "processes": 8,
}

instances = [
    # FSC, run fsc with removed graph edges first to ensure that the chosen circutis really work. i encountered 1 case when the crossing metric could not be computed but could not reprodcue it
    # f=12
    {
        "q": 42,
        "t": 4,
        "min_depth": 42 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q42_row_f12_FSC,
        "factory_locs": factories_q42_row,
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "FSC",
        "circuit_type": "random",
    },
    {
        "q": 42,
        "t": 12,
        "min_depth": 42 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q42_row_f12_FSC,
        "factory_locs": factories_q42_row,
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "FSC",
        "circuit_type": "random",
    },
    # f=4
    {
        "q": 42,
        "t": 4,
        "min_depth": 42 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q42_row_f4_FSC,
        "factory_locs": factories_q42_row[:4],
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "FSC",
        "circuit_type": "random",
    },
    {
        "q": 42,
        "t": 12,
        "min_depth": 42 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q42_row_f4_FSC,
        "factory_locs": factories_q42_row[:4],
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "FSC",
        "circuit_type": "random",
    },
    # CC
    # f=12
    {
        "q": 42,
        "t": 4,
        "min_depth": 42 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q42_row_f12_CC,
        "factory_locs": factories_q42_row,
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "CC",
        "circuit_type": "random",
    },
    {
        "q": 42,
        "t": 12,
        "min_depth": 42 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q42_row_f12_CC,
        "factory_locs": factories_q42_row,
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "CC",
        "circuit_type": "random",
    },
    # f=4
    {
        "q": 42,
        "t": 4,
        "min_depth": 42 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q42_row_f4_CC,
        "factory_locs": factories_q42_row[:4],
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "CC",
        "circuit_type": "random",
    },
    {
        "q": 42,
        "t": 12,
        "min_depth": 42 * 4,
        "tgate": True,
        "ratio": 0.8,
        "custom_layout": custom_layout_q42_row_f4_CC,
        "factory_locs": factories_q42_row[:4],
        "layout_type": "custom",
        "layout_name": "row",
        "graphtype": "CC",
        "circuit_type": "random",
    },
]


reps = 10
both_metric = True
res_lst = co.plots.collect_data_space_time(instances, hc_params, reps, path, both_metric)


with Path(path).open("rb") as f:
    res_lst = pickle.load(f)

with Path(path).open("rb") as f:
    res_lst = pickle.load(f)

path = "./results/f_vs_time_q42_ratio08_small_row_metricrouting"

with Path(path).open("rb") as f:
    res_lst_routing = pickle.load(f)

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
