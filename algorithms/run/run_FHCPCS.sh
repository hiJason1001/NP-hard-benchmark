#!/bin/bash

set -m

filename="HPA"

executable="./../${filename}.out"
input_folder="../../data_processed/FHCPCS_processed"
output_file="${filename}_FHCPCS_processed.txt"
time_limit_minutes=10  # set timeout minutes
time_limit_seconds=$((time_limit_minutes * 60))

> "$output_file"


for i in {1..1001}
do
	file="$input_folder/graph${i}.hcp"
    echo "$(basename "$file") " >> "$output_file"

    timeout "$time_limit_seconds" "$executable" "$file" >> "$output_file"

    status=$?
    if [ $status -eq 124 ]; then
        echo "TimeLimitExceeded" >> "$output_file"
    fi
done