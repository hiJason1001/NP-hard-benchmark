#!/bin/bash

filename="NN"

executable="./${filename}"
input_folder="NP-hard-benchmark/data_processed/examples"
output_file="${filename}_examples.txt"

> "$output_file"

for i in {1..100}
do
    file="$input_folder/example${i}.hcp"
    echo "example${i}.hcp " >> "$output_file"
    "$executable" "$file" >> "$output_file"
done