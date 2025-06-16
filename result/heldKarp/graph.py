import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ========================
# CONFIGURATION
# ========================
NAME = "heldKarp_examples"
FOLDER = "heldKarp"
TIME_UNIT = "seconds"
TITLE = f"Held-Karp Runtime vs Number of Vertices ({TIME_UNIT})"
INPUT_FILE = f"result/{FOLDER}/{NAME}.txt"
OUTPUT_FILE = f"result/{FOLDER}/{NAME}.png"
FONT = 24

# ========================

with open(INPUT_FILE, 'r') as file:
    data = file.read()

lines = data.strip().splitlines()

records = []
i = 0

while i < len(lines):
    line = lines[i].strip()

    if line.endswith(".hcp"):
        if i + 3 >= len(lines):
            break

        filename = line
        vertices_edges_line = lines[i + 1].strip()
        time_line = lines[i + 2].strip()
        result_line = lines[i + 3].strip().lower()

        vertices_match = re.search(r'Number of Vertices:\s*(\d+),\s*number of edges:\s*(\d+)', vertices_edges_line)
        if not vertices_match:
            print(f"Warning: could not parse vertices/edges at line {i+1}: {vertices_edges_line}")
            i += 4
            continue

        vertices = int(vertices_match.group(1))
        edges = int(vertices_match.group(2))

        time_match = re.search(r'(\d+) hours (\d+) minutes (\d+) seconds (\d+) milliseconds (\d+) microseconds', time_line)
        if not time_match:
            print(f"Warning: could not parse time at line {i+2}: {time_line}")
            i += 4
            continue

        hours = int(time_match.group(1))
        minutes = int(time_match.group(2))
        seconds = int(time_match.group(3))
        milliseconds = int(time_match.group(4))
        microseconds = int(time_match.group(5))

        total_seconds = (hours * 3600) + (minutes * 60) + seconds + (milliseconds / 1000) + (microseconds / 1_000_000)

        if TIME_UNIT == "seconds":
            total_time = total_seconds
        elif TIME_UNIT == "milliseconds":
            total_time = total_seconds * 1000
        elif TIME_UNIT == "microseconds":
            total_time = total_seconds * 1_000_000
        else:
            raise ValueError(f"Unsupported TIME_UNIT: {TIME_UNIT}")

        is_no = 'no' in result_line

        records.append({
            'vertices': vertices,
            'edges': edges,
            'time': total_time,
            'result': 'no' if is_no else 'yes'
        })

        i += 4
    else:
        i += 1

# ========================
# PLOT
# ========================

if not records:
    print("No records found! Please check the data format.")
else:
    df = pd.DataFrame(records)

    df_yes = df[df["result"] == "yes"]
    df_no = df[df["result"] == "no"]

    plt.figure(figsize=(10, 6))

    plt.scatter(df_yes["vertices"], df_yes["time"], color='blue', label='YES instances')

    plt.scatter(df_no["vertices"], df_no["time"], color='red', label='NO instances')

    if len(df) >= 2:
        x = df["vertices"]
        y = df["time"]

        x_fit = np.linspace(df["vertices"].min(), df["vertices"].max(), 200)
        coeffs = np.polyfit(x, np.log(y), 1)
        exp_fit = np.exp(coeffs[1]) * np.exp(coeffs[0] * x_fit)

        plt.plot(x_fit, exp_fit, color='black', label='Exponential Fit')

    plt.xlabel("Number of Vertices", fontsize=FONT)
    plt.ylabel(f"Time ({TIME_UNIT})", fontsize=FONT)
    plt.title(TITLE, fontsize=FONT)

    plt.ylim(bottom=-0.05 * df["time"].max())
    plt.grid(True)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE)
    print(f"Graph saved as {OUTPUT_FILE}")
