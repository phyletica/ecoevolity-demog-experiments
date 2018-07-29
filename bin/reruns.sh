#!/bin/bash

bin_dir="$(pwd)"

while read qsub_path
do
    dir_path="$(dirname "$qsub_path")"
    file_name="$(basename "$qsub_path")"

    cd "$dir_path"
    ls "$file_name"

    prefix="${file_name/-qsub\.sh/}"
    run_number="${prefix##*run-}"
    sim_base="${prefix%-run-*}"

    op_log_file="run-${run_number}-${sim_base}-operator-run-1.log"
    state_log_file="run-${run_number}-${sim_base}-state-run-1.log"
    stdout_file="run-${run_number}-${sim_base}.yml.out"

    if [ -e "$op_log_file" ]
    then
        rm "$op_log_file"
    fi

    if [ -e "$state_log_file" ]
    then
        rm "$state_log_file"
    fi

    if [ -e "$stdout_file" ]
    then
        rm "$stdout_file"
    fi

    # nsub -r -t "05:00:00" "$file_name"
    nsub -t "02:00:00" "$file_name"

    cd "$bin_dir"

done < "reruns.txt"
