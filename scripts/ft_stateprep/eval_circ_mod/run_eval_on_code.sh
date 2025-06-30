#!/bin/bash
# Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

# Description: Estimate logical error rates for given code in parallel. Results are saved in a csv file described by the first argument. All other arguments are passed to the python script.

declare -a p=("0.5" "0.3" "0.1" "0.09" "0.07" "0.05" "0.03" "0.01" "0.009" "0.007" "0.005" "0.003" "0.001" "0.0009" "0.0007" "0.0005" "0.0003" "0.0001")

echo "p p_l acceptance errors runs p_l_error acceptance_error" > "$1.csv"

run_and_write() {
    local res=$(python estimate_logical_error_rate.py ${@:2:$#-2} "-p" "${@: -1}")
    local line="${@: -1} ${res}"
    (flock -e 200 echo $line >> "$1.csv") 200>lock
}

export -f run_and_write

parallel --load 16 --link run_and_write $@ ::: ${p[@]}
