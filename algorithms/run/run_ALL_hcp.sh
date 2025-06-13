#!/bin/bash

set -m

filename="HPA"

executable="./../${filename}.out"
input_folder="../../data_processed/ALL_hcp_processed"
output_file="${filename}_ALL_hcp.txt"
time_limit_minutes=10  # set timeout minutes
time_limit_seconds=$((time_limit_minutes /10))

> "$output_file"


for file in "$input_folder"/*.hcp
do
    echo "$(basename "$file") " >> "$output_file"

    timeout "$time_limit_seconds" "$executable" "$file" >> "$output_file"

    status=$?
    if [ $status -eq 124 ]; then
        echo "Time Limit Exceeded" >> "$output_file"
    fi
done