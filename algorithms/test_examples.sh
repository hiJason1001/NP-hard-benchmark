#!/bin/bash

executable="./heldKarp.out"
input_folder="../data_processed/examples"
output_file="answer.txt"

# Clear the output file before starting
> "$output_file"

for i in {1..100}
do
    file="$input_folder/example${i}.hcp"
    echo -n "example${i}.hcp " >> "$output_file"
    "$executable" "$file" >> "$output_file"
done