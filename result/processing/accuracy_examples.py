import re

NAME = "HPA_examples"

ACTUAL_FILE = "result/exmaples_actual_results.txt"

INPUT_FILE = f"result/HPA/{NAME}.txt"
OUTPUT_FILE = f"{NAME}_accuracy.txt"

with open(INPUT_FILE, 'r') as file:
    experimental_data = file.read()

with open(ACTUAL_FILE, 'r') as file:
    actual_data = file.read()

lines = experimental_data.strip().splitlines()

experimental_results = {}
i = 0
while i < len(lines):
    line = lines[i].strip()
    
    if line.endswith(".hcp"):
        if i + 2 >= len(lines):
            break
        
        vertices_edges_line = lines[i + 1].strip()
        time_line = lines[i + 2].strip()
        yes_no_line = lines[i + 3].strip() 
        
        if yes_no_line in ["Yes", "No"]:
            experimental_results[line] = yes_no_line
        else:
            i += 4
            continue
        
        i += 4
    else:
        i += 1

actual_results = {}
lines_actual = actual_data.strip().splitlines()

for line in lines_actual:
    parts = line.split()
    if len(parts) >= 2:
        actual_results[parts[0]] = parts[1]

correct_matches = 0
incorrect_matches = 0
correct_yes = 0
incorrect_yes = 0
correct_no = 0
incorrect_no = 0

for test_case, exp_result in experimental_results.items():
    actual_result = actual_results.get(test_case)
    
    if actual_result:
        if exp_result == actual_result:
            correct_matches += 1
            if actual_result == "Yes":
                correct_yes += 1
            else:
                correct_no += 1
        else:
            incorrect_matches += 1
            if actual_result == "Yes":
                incorrect_yes += 1
            else:
                incorrect_no += 1

total_test_cases = len(experimental_results)
accuracy_percentage = (correct_matches / total_test_cases) * 100 if total_test_cases > 0 else 0

accuracy_metrics = f"""Accuracy: {accuracy_percentage:.2f}%
Total Test Cases: {total_test_cases}
Correct Matches: {correct_matches}
Incorrect Matches: {incorrect_matches}

Correct Matches with Yes: {correct_yes}
Incorrect Matches with Yes: {incorrect_yes}

Correct Matches with No: {correct_no}
Incorrect Matches with No: {incorrect_no}
"""

with open(OUTPUT_FILE, 'w') as output_file:
    output_file.write(accuracy_metrics)

print(f"Accuracy metrics saved to {OUTPUT_FILE}")
