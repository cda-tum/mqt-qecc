from __future__ import annotations

import itertools
import json
import os
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class BpParams:
    """Class to store parameters for BP decoding."""

    bp_method: str = "msl"
    max_bp_iter: int = 30
    osd_order: int = 10
    osd_method: str = "osd_cs"
    ms_scaling_factor: float = 0.75
    schedule: str = "parallel"
    omp_thread_count: int = 1
    random_serial_schedule: int = 0
    serial_schedule_order: list[int] | None = None
    cutoff: float = np.inf

    @classmethod
    def from_dict(cls, dict_):
        class_fields = {f.name for f in fields(cls)}
        return BpParams(**{k: v for k, v in dict_.items() if k in class_fields})


def extract_settings(filename):
    """Extracts all settings from a parameter file and returns them as a dictionary."""
    keyword_lists = {}

    with open(filename) as file:
        for line in file:
            json_data = json.loads(line.strip())
            for keyword, value in json_data.items():
                if keyword not in keyword_lists:
                    keyword_lists[keyword] = []
                if value not in keyword_lists[keyword]:
                    keyword_lists[keyword].append(value)

    return keyword_lists


def load_data(
    input_filenames: list[str],
) -> list[dict]:
    """Loads data from a list of JSON files and returns it as a list of dictionaries."""
    data = []
    for file in input_filenames:
        path = Path(file)

        try:
            ldata = json.load(path.open())
            data.append(ldata)
        except:
            merge_json_files(path.with_suffix(""))
            ldata = json.load(path.open())
            data.append(ldata)
    return data


def calculate_error_rates(success_cnt, runs, code_params):
    """Calculates logical error rate, logical error rate error bar, word error rate,
    and word error rate error bar.
    """
    logical_err_rate = 1.0 - (success_cnt / runs)
    logical_err_rate_eb = np.sqrt((1 - logical_err_rate) * logical_err_rate / runs)
    word_error_rate = 1.0 - (1 - logical_err_rate) ** (1 / code_params["k"])
    word_error_rate_eb = (
        logical_err_rate_eb * ((1 - logical_err_rate_eb) ** (1 / code_params["k"] - 1)) / code_params["k"]
    )
    return (
        logical_err_rate,
        logical_err_rate_eb,
        word_error_rate,
        word_error_rate_eb,
    )


def is_converged(x_success, z_success, runs, code_params, precission):
    x_cond = _check_convergence(x_success, runs, code_params, precission_cutoff=precission)
    z_cond = _check_convergence(z_success, runs, code_params, precission_cutoff=precission)
    return x_cond == z_cond is True


def _check_convergence(success_cnt, runs, code_params, precission_cutoff):
    _, _, ler, ler_eb = calculate_error_rates(success_cnt, runs, code_params)
    if ler_eb != 0.0:
        if ler_eb / ler < precission_cutoff:
            return True
        return None
    else:
        return False


def create_outpath(
    x_meta: bool = False,
    z_meta: bool = False,
    bias: list[float] | None = None,
    codename: str | None = None,
    single_stage: bool = True,
    sus_th_depth: int | None = None,
    rounds: int | None = None,
    id: int = 0,
    overwrite: bool = False,
    analog_info: bool = False,
    analog_tg: bool = False,
    repetitions: int | None = None,
    experiment: str = "wer_per_round",
    **kwargs,
) -> str:
    path = f"results/{experiment:s}/"
    if analog_info:
        path += "analog_info/"  # ranvendraan et al analog info decoder
    elif analog_tg:
        path += "analog_tg/"  # analog tannergraph bp
    else:
        path += "hard_syndrome/"  # standard hard syndrome decoding only
    if bias is not None:
        path += f"single_stage={single_stage}/bias={bias[0]}_{bias[1]}_{bias[2]}/"

    if sus_th_depth:
        path += f"sus_th_depth={sus_th_depth}/"
    elif rounds:
        path += f"rounds={rounds}/"
    # else:
    #     raise ValueError(
    #         "Either sus_th_depth or window_count_per_run must be specified."
    #     )

    if repetitions:
        path += f"repetitions={repetitions}/"

    if x_meta:
        path += "x-meta=true/"
    else:
        path += "x-meta=false/"

    if z_meta:
        path += "z-meta=true/"
    else:
        path += "z-meta=false/"

    path += f"{codename:s}/"

    if "syndr_err_rate" not in kwargs or kwargs["syndr_err_rate"] is None:
        if "sigma" in kwargs:
            path += f"per_{kwargs['data_err_rate']:.3e}_sigma_{kwargs['sigma']:.3e}/"
        if "z_sigma" in kwargs:
            path += f"per_{kwargs['data_err_rate']:.3e}_x_sigma_{kwargs['x_sigma']:.3e}_z_sigma_{kwargs['z_sigma']:.3e}"
    else:
        path += f"per_{kwargs['data_err_rate']:.3e}_ser_{kwargs['syndr_err_rate']:.3e}/"

    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

    if overwrite is False:
        f_loc = path + f"id_{id}.json"
        while os.path.exists(f_loc):
            id += 1
            f_loc = path + f"id_{id}.json"

    while not os.path.exists(f_loc):
        open(f_loc, "w").close()

    return f_loc


def replace_inf(lst: list):
    new_lst = []
    for item in lst:
        if np.isinf(item):
            new_lst.append("i")
        else:
            new_lst.append(item)
    return new_lst


def product_dict(**kwargs):
    """Generate a iterator of dictionaries where each dictionary is a cartesian product
    of the values associated with each key in the input dictionary.
    """
    keys = kwargs.keys()
    vals = kwargs.values()
    for instance in itertools.product(*vals):
        yield dict(zip(keys, instance))


def zip_dict(**kwargs):
    """Create a iterator of dictionaries where each dictionary contains the zip() of the
    values associated with each key in the input dictionary.
    """
    return (dict(zip(kwargs.keys(), values)) for values in zip(*kwargs.values()))


def _update_error_rates(success_cnt, runs, code_K):
    """Calculates logical error rate, logical error rate error bar, word error rate,
    and word error rate error bar.
    """
    logical_err_rate = 1.0 - (success_cnt / runs)
    logical_err_rate_eb = np.sqrt((1 - logical_err_rate) * logical_err_rate / runs)
    word_error_rate = 1.0 - (1 - logical_err_rate) ** (1 / code_K)
    word_error_rate_eb = logical_err_rate_eb * ((1 - logical_err_rate_eb) ** (1 / code_K - 1)) / code_K
    return (
        logical_err_rate,
        logical_err_rate_eb,
        word_error_rate,
        word_error_rate_eb,
    )


def merge_datasets(datasets: list[dict[str, Any]]) -> dict[str, Any]:
    """Merges a list of dictionaries into a single dictionary. The values for the fields "nr_runs",
    "x_success_cnt" and "z_success_cnt" are extracted from each dictionary and added together.

    Args:
        datasets (List[Dict[str, Any]]): A list of dictionaries to be merged.

    Returns:
        Dict[str, Any]: A dictionary containing the merged data.
    """
    if not datasets:
        return {}

    # Start with a copy of the first dictionary in the list
    merged_data = dict(datasets[0])

    # Extract and add up the values for "nr_runs", "x_success_cnt", and "z_success_cnt"
    for data in datasets[1:]:
        merged_data["nr_runs"] += data.get("nr_runs", 0)
        merged_data["x_success_cnt"] += data.get("x_success_cnt", 0)
        merged_data["z_success_cnt"] += data.get("z_success_cnt", 0)

    # Update logical and word error rates based on accumulated data.
    runs = merged_data["nr_runs"]
    x_success_cnt = merged_data["x_success_cnt"]
    z_success_cnt = merged_data["z_success_cnt"]

    x_ler, x_ler_eb, x_wer, x_wer_eb = _update_error_rates(x_success_cnt, runs, merged_data["code_K"])
    z_ler, z_ler_eb, z_wer, z_wer_eb = _update_error_rates(z_success_cnt, runs, merged_data["code_K"])

    update_dict = {
        "x_ler": x_ler,
        "x_ler_eb": x_ler_eb,
        "x_wer": x_wer,
        "x_wer_eb": x_wer_eb,
        "x_success_cnt": x_success_cnt,
        "z_ler": z_ler,
        "z_ler_eb": z_ler_eb,
        "z_wer": z_wer,
        "z_wer_eb": z_wer_eb,
        "z_success_cnt": z_success_cnt,
    }
    merged_data.update(update_dict)
    return merged_data


def _merge_datasets_x(_datasets: list[dict[str, Any]]) -> dict[str, Any]:
    """Merges a list of dictionaries into a single dictionary. The values for the fields "nr_runs",
    "x_success_cnt" and "z_success_cnt" are extracted from each dictionary and added together.

    Args:
        datasets (List[Dict[str, Any]]): A list of dictionaries to be merged.

    Returns:
        Dict[str, Any]: A dictionary containing the merged data.
    """
    datasets = _datasets.copy()

    if not datasets:
        return {}

    # remove datasets that do not contain x_success_cnt
    for i, data in enumerate(datasets):
        if "x_success_cnt" not in data:
            datasets.pop(i)

    # Start with a copy of the first dictionary in the list that contains z_success_cnt
    # and remove that dict from the list
    for i, data in enumerate(datasets):
        if "x_success_cnt" in data:
            merged_data = dict(datasets.pop(i))
            break

    # Extract and add up the values for "nr_runs", "x_success_cnt", and "z_success_cnt"
    for data in datasets:
        try:
            merged_data["nr_runs"] += data.get("nr_runs", 0)
            merged_data["x_success_cnt"] += data.get("x_success_cnt", 0)
        # merged_data["z_success_cnt"] += data.get("z_success_cnt", 0)
        except:
            pass

    # Update logical and word error rates based on accumulated data.
    runs = merged_data["nr_runs"]
    x_success_cnt = merged_data["x_success_cnt"]
    # z_success_cnt = merged_data["z_success_cnt"]

    x_ler, x_ler_eb, x_wer, x_wer_eb = _update_error_rates(x_success_cnt, runs, merged_data["code_K"])

    update_dict = {
        "x_ler": x_ler,
        "x_ler_eb": x_ler_eb,
        "x_wer": x_wer,
        "x_wer_eb": x_wer_eb,
        "x_success_cnt": x_success_cnt,
    }
    merged_data.update(update_dict)
    return merged_data


def _merge_datasets_z(_datasets: list[dict[str, Any]]) -> dict[str, Any]:
    """Merges a list of dictionaries into a single dictionary. The values for the fields "nr_runs",
    "x_success_cnt" and "z_success_cnt" are extracted from each dictionary and added together.

    Args:
        datasets (List[Dict[str, Any]]): A list of dictionaries to be merged.

    Returns:
        Dict[str, Any]: A dictionary containing the merged data.
    """
    datasets = _datasets.copy()
    # print("datasets", datasets)
    if not datasets:
        return {}

    # remove datasets that do not contain z_success_cnt
    for i, data in enumerate(datasets):
        if "z_success_cnt" not in data:
            datasets.pop(i)

    # print(datasets)

    # Start with a copy of the first dictionary in the list that contains z_success_cnt
    # and remove that dict from the list
    for i, data in enumerate(datasets):
        if "z_success_cnt" in data:
            merged_data = dict(datasets.pop(i))
            break

    # merged_data = dict(datasets[0])

    # Extract and add up the values for "nr_runs", "x_success_cnt", and "z_success_cnt"
    for data in datasets:
        try:
            merged_data["nr_runs"] += data.get("nr_runs", 0)
            # merged_data["z_success_cnt"] += data.get("z_success_cnt", 0)
            merged_data["z_success_cnt"] += data.get("z_success_cnt", 0)
        except:
            pass

    # Update logical and word error rates based on accumulated data.
    runs = merged_data["nr_runs"]
    # x_success_cnt = merged_data["x_success_cnt"]
    z_success_cnt = merged_data["z_success_cnt"]

    z_ler, z_ler_eb, z_wer, z_wer_eb = _update_error_rates(z_success_cnt, runs, merged_data["code_K"])

    update_dict = {
        "z_ler": z_ler,
        "z_ler_eb": z_ler_eb,
        "z_wer": z_wer,
        "z_wer_eb": z_wer_eb,
        "z_success_cnt": z_success_cnt,
    }
    merged_data.update(update_dict)
    return merged_data


from json.decoder import JSONDecodeError


def merge_json_files(input_path: str) -> None:
    """Iterates through all subfolders in the input folder, loads all JSON files in each subfolder
    and merges them using the `merge_datasets` function. The resulting dictionaries are
    stored in a list and then dumped as a JSON file in the input folder.

    Args:
        input_path (str): The path to the input folder.

    Returns:
        None
    """
    output_data: list[dict[str, Any]] = []
    for folder_name in os.listdir(input_path):
        folder_path = os.path.join(input_path, folder_name)
        if os.path.isdir(folder_path):
            data: list[dict[str, Any]] = []
            for filename in os.listdir(folder_path):
                if filename.endswith(".json"):
                    file_path = os.path.join(folder_path, filename)
                    with open(file_path) as file:
                        try:
                            json_data = json.load(file)
                            data.append(json_data)
                        except JSONDecodeError:
                            pass
            merged_data = merge_datasets(data)
            if merged_data:
                output_data.append(merged_data)

    # save output to parent directiory
    code_name = os.path.basename(os.path.normpath(input_path))
    parent_dir = os.path.abspath(os.path.join(input_path, os.pardir))
    with open(os.path.join(parent_dir, f"{code_name:s}.json"), "w") as output_file:
        json.dump(output_data, output_file, ensure_ascii=False, indent=4)


def merge_json_files_x(input_path: str) -> None:
    """Iterates through all subfolders in the input folder, loads all JSON files in each subfolder
    and merges them using the `merge_datasets` function. The resulting dictionaries are
    stored in a list and then dumped as a JSON file in the input folder.

    Args:
        input_path (str): The path to the input folder.

    Returns:
        None
    """
    output_data: list[dict[str, Any]] = []
    for folder_name in os.listdir(input_path):
        folder_path = os.path.join(input_path, folder_name)
        if os.path.isdir(folder_path):
            data: list[dict[str, Any]] = []
            for filename in os.listdir(folder_path):
                if filename.endswith(".json"):
                    file_path = os.path.join(folder_path, filename)
                    with open(file_path) as file:
                        try:
                            json_data = json.load(file)
                            data.append(json_data)
                        except JSONDecodeError:
                            pass
            merged_data = _merge_datasets_x(data)
            if merged_data:
                output_data.append(merged_data)

    # save output to parent directiory
    code_name = os.path.basename(os.path.normpath(input_path))
    parent_dir = os.path.abspath(os.path.join(input_path, os.pardir))
    with open(os.path.join(parent_dir, f"{code_name:s}.json"), "w") as output_file:
        json.dump(output_data, output_file, ensure_ascii=False, indent=4)


def merge_json_files_z(input_path: str) -> None:
    """Iterates through all subfolders in the input folder, loads all JSON files in each subfolder
    and merges them using the `merge_datasets` function. The resulting dictionaries are
    stored in a list and then dumped as a JSON file in the input folder.

    Args:
        input_path (str): The path to the input folder.

    Returns:
        None
    """
    output_data: list[dict[str, Any]] = []
    for folder_name in os.listdir(input_path):
        folder_path = os.path.join(input_path, folder_name)
        if os.path.isdir(folder_path):
            data: list[dict[str, Any]] = []
            for filename in os.listdir(folder_path):
                if filename.endswith(".json"):
                    file_path = os.path.join(folder_path, filename)
                    with open(file_path) as file:
                        try:
                            json_data = json.load(file)
                            data.append(json_data)
                        except JSONDecodeError:
                            pass
            merged_data = _merge_datasets_z(data)
            if merged_data:
                output_data.append(merged_data)

    # save output to parent directiory
    code_name = os.path.basename(os.path.normpath(input_path))
    parent_dir = os.path.abspath(os.path.join(input_path, os.pardir))
    with open(os.path.join(parent_dir, f"{code_name:s}.json"), "w") as output_file:
        json.dump(output_data, output_file, ensure_ascii=False, indent=4)


def merge_json_files_xz(input_path: str) -> None:
    """Iterates through all subfolders in the input folder, loads all JSON files in each subfolder
    and merges them using the `merge_datasets` function. The resulting dictionaries are
    stored in a list and then dumped as a JSON file in the input folder.

    Args:
        input_path (str): The path to the input folder.

    Returns:
        None
    """
    output_data: list[dict[str, Any]] = []
    for folder_name in os.listdir(input_path):
        folder_path = os.path.join(input_path, folder_name)
        if os.path.isdir(folder_path):
            data: list[dict[str, Any]] = []
            for filename in os.listdir(folder_path):
                if filename.endswith(".json"):
                    file_path = os.path.join(folder_path, filename)
                    with open(file_path) as file:
                        try:
                            json_data = json.load(file)
                            data.append(json_data)
                        except JSONDecodeError:
                            pass
            # print(folder_path, filename)
            # print(data)
            merged_data_x = _merge_datasets_x(data)
            merged_data_z = _merge_datasets_z(data)
            merged_data = _combine_xz_data(merged_data_x, merged_data_z)
            if merged_data:
                output_data.append(merged_data)

    # save output to parent directiory
    code_name = os.path.basename(os.path.normpath(input_path))
    parent_dir = os.path.abspath(os.path.join(input_path, os.pardir))
    with open(os.path.join(parent_dir, f"{code_name:s}.json"), "w") as output_file:
        json.dump(output_data, output_file, ensure_ascii=False, indent=4)


def _combine_xz_data(xdata: dict | None, zdata: dict | None) -> dict:
    """Combine the x and z data into a single dictionary.
    Before doing that, rename "runs" in each dictionary to "x_runs" and "z_runs" respectively.

    """
    if xdata and zdata:
        # print(xdata)
        xdata["x_runs"] = xdata.pop("nr_runs")
        zdata["z_runs"] = zdata.pop("nr_runs")
        xdata.update(zdata)
        return xdata
    elif xdata:
        xdata["x_runs"] = xdata.pop("nr_runs")
        return xdata
    elif zdata:
        zdata["z_runs"] = zdata.pop("nr_runs")
        return zdata
    else:
        return {}
