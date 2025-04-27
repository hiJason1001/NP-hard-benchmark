#!/bin/bash

executable="./heldKarp.out"
input_folder="../data_processed/examples"
output_file="heldKarp.txt"

# Clear the output file before starting
> "$output_file"

for i in {1..100}
do
    file="$input_folder/example${i}.hcp"
    echo "example${i}.hcp " >> "$output_file"
    "$executable" "$file" >> "$output_file"
done