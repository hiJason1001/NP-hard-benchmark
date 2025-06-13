input_file = "result/Concorde/concorde_examples.txt"
output_file = "concorde_examples_sorted.txt"

def extract_number(filename):
    import re
    match = re.search(r'example(\d+)\.hcp', filename)
    return int(match.group(1)) if match else float('inf')

with open(input_file, "r") as f:
    lines = [line for line in f if line.strip() != ""]


records = [lines[i:i+4] for i in range(0, len(lines), 4)]

records.sort(key=lambda rec: extract_number(rec[0]))

with open(output_file, "w") as f:
    for rec in records:
        f.writelines(rec)
        if not rec[-1].endswith('\n'):
            f.write('\n')