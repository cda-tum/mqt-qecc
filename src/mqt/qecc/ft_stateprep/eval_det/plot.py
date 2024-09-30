# %%
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mqt.qecc.ft_stateprep import DeterministicVerification

# %%
data = pd.read_csv("results.csv", index_col=0)
data.reset_index(inplace=True, drop=True)
data.head()


# %% [markdown]
# # 1. Compute numbers for the table

# %%
eval_names = {"array": np.array, "int8": np.int8, "nan": np.nan}
for index, row in data.iterrows():
    verify_0 = DeterministicVerification(
        eval(row["verification_stabs_0"], eval_names),
        eval(row["recovery_stabs_0"], eval_names),
        eval(row["flags_0"], eval_names),
    )
    verify_1 = DeterministicVerification(
        eval(row["verification_stabs_1"], eval_names),
        eval(row["recovery_stabs_1"], eval_names),
        eval(row["flags_1"], eval_names),
    )
    # add column
    funs = [min, max, sum, np.mean]
    data.at[index, "verify_num_anc_0"] = verify_0.num_ancillas_verification()
    data.at[index, "verify_num_cnots_0"] = verify_0.num_cnots_verification()
    data.at[index, "verify_num_anc_1"] = verify_1.num_ancillas_verification()
    data.at[index, "verify_num_cnots_1"] = verify_1.num_cnots_verification()
    data.at[index, "correction_num_anc_0"] = str([verify_0.stat_ancillas_correction(fun) for fun in funs])
    data.at[index, "correction_num_anc_list_0"] = str(verify_0.stat_ancillas_correction(list))
    data.at[index, "correction_num_cnots_0"] = str([verify_0.stat_cnots_correction(fun) for fun in funs])
    data.at[index, "correction_num_cnots_list_0"] = str(verify_0.stat_cnots_correction(list))
    data.at[index, "correction_num_anc_1"] = str([verify_1.stat_ancillas_correction(fun) for fun in funs])
    data.at[index, "correction_num_anc_list_1"] = str(verify_1.stat_ancillas_correction(list))
    data.at[index, "correction_num_cnots_1"] = str([verify_1.stat_cnots_correction(fun) for fun in funs])
    data.at[index, "correction_num_cnots_list_1"] = str(verify_1.stat_cnots_correction(list))
    data.at[index, "flag_num_anc_0"] = verify_0.num_ancillas_hooks()
    data.at[index, "flag_num_cnots_0"] = verify_0.num_cnots_hooks()
    data.at[index, "flag_num_anc_1"] = verify_1.num_ancillas_hooks()
    data.at[index, "flag_num_cnots_1"] = verify_1.num_cnots_hooks()
    data.at[index, "flag_correction_num_anc_0"] = str([verify_0.stat_ancillas_hook_corrections(fun) for fun in funs])
    data.at[index, "flag_correction_num_anc_list_0"] = str(verify_0.stat_ancillas_hook_corrections(list))
    data.at[index, "flag_correction_num_cnots_0"] = str([verify_0.stat_cnots_hook_corrections(fun) for fun in funs])
    data.at[index, "flag_correction_num_cnots_list_0"] = str(verify_0.stat_cnots_hook_corrections(list))
    data.at[index, "flag_correction_num_anc_1"] = str([verify_1.stat_ancillas_hook_corrections(fun) for fun in funs])
    data.at[index, "flag_correction_num_anc_list_1"] = str(verify_1.stat_ancillas_hook_corrections(list))
    data.at[index, "flag_correction_num_cnots_1"] = str([verify_1.stat_cnots_hook_corrections(fun) for fun in funs])
    data.at[index, "flag_correction_num_cnots_list_1"] = str(verify_1.stat_cnots_hook_corrections(list))

# %%
# remove above columns
data_table = data.copy()

stats_columns = [
    "verification_stabs_0",
    "recovery_stabs_0",
    "flags_0",
    "verification_stabs_1",
    "recovery_stabs_1",
    "flags_1",
]
error_rates_columns = ["logical_error_rates"]
data_table = data_table.drop(columns=stats_columns)
data_table = data_table.drop(columns=error_rates_columns)
# data_table

# %% [markdown]
# # 2. Create Plot with logical error rates

# %%
fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(111)

physical_error_rates = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1]
# fmts = ['o--', 's--', 'd--', 'x--', 'v--', '^--', '<--', '>--' , 'o--']
fmts = ["o:", "s:", "d:", "x:", "v:", "^:", "<:", ">:", "o:"]
fmt_idx = 0

data_sorted = data.sort_values(by="code")

sim_col = "logical_error_rates"
for i, row in data_sorted.iterrows():
    if not row["zero_state"] and row["code"] != "tetrahedral":
        continue
    if row["zero_state"] and row["code"] == "tetrahedral":
        continue
    if row["procedure"] != "heuristic" or row["verification"] != "optimal":
        continue
    code_name = row["code"].capitalize()
    if code_name == "16_2_4":
        code_name = r"$[[16, 2, 4]]$"
    elif code_name == "11_1_3":
        code_name = r"$[[11, 1, 3]]$"
    elif code_name == "Hypercube":
        code_name = "Tesseract"
    upper_bound, std = eval(row[sim_col], eval_names)[-2:]
    ax.errorbar(physical_error_rates, upper_bound, label=code_name, fmt=fmts[fmt_idx])
    fmt_idx += 1


# add linear line as reference
ax.plot(physical_error_rates, physical_error_rates, label="Linear", linestyle="--", color="black", alpha=0.5)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylabel("Logical error rate $p_{\mathrm{L}}$")
ax.set_xlabel(r"Physical error rate $p$")
ax.set_xlim(1e-4, 0.5)
ax.set_ylim(1e-7, 1)

ax.legend(loc="lower right")
# define order of legend
handles, labels = ax.get_legend_handles_labels()
order = [0, 7, 6, 8, 1, 9, 4, 3, 2, 5]
ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc="lower right")

plt.show()
fig.savefig("logical_error_rates.pdf", bbox_inches="tight")
