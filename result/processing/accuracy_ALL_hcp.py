# ========================
# CONFIGURATION
# ========================
EXPERIMENTAL_FILE = "result/ALL_hcp/hybridHam_ALL_hcp.txt"
OUTPUT_FILE = "accuracy.txt"
# ========================

with open(EXPERIMENTAL_FILE, 'r') as file:
    experimental_data = file.read()

lines = experimental_data.strip().splitlines()

experimental_results = {}
i = 0
while i < len(lines):
    line = lines[i].strip()
    
    if line.endswith(".hcp"):
        if i + 3 >= len(lines):
            break
        
        yes_no_line = lines[i + 3].strip() 
        
        if yes_no_line in ["Yes", "No"]:
            experimental_results[line] = yes_no_line
            i += 4
        else:
            i += 1
    else:
        i += 1

correct_matches = 0
incorrect_matches = 0

for test_case, exp_result in experimental_results.items():
    if exp_result == "Yes":
        correct_matches += 1
    else:
        incorrect_matches += 1

total_test_cases = len(experimental_results)
accuracy_percentage = (correct_matches / total_test_cases) * 100 if total_test_cases > 0 else 0

accuracy_metrics = f"""Accuracy: {accuracy_percentage:.2f}%
Total Test Cases: {total_test_cases}
Correct (Expected "Yes"): {correct_matches}
Incorrect (Got "No" instead of "Yes"): {incorrect_matches}
"""

with open(OUTPUT_FILE, 'w') as output_file:
    output_file.write(accuracy_metrics)

print(f"Accuracy metrics saved to {OUTPUT_FILE}")
