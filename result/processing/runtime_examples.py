import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ========================
# CONFIGURATION
# ========================
# Options: "seconds", "milliseconds", "microseconds"
TIME_UNIT = "microseconds"

INPUT_FILE = "result/examples/hybridHam.txt"
OUTPUT_FILE = "hybridHam_examples.png"

# ========================



with open(INPUT_FILE, 'r') as file:
    data = file.read()

lines = data.strip().splitlines()

records = []
i = 0

while i < len(lines):
    line = lines[i].strip()
    
    if line.endswith(".hcp"):
        if i + 2 >= len(lines):
            break
        
        vertices_edges_line = lines[i + 1].strip()
        time_line = lines[i + 2].strip()
        
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
        
        records.append({
            'vertices': vertices,
            'edges': edges,
            'time': total_time
        })
        
        i += 4
    else:
        i += 1


if not records:
    print("No records found! Please check the data format.")
else:
    df = pd.DataFrame(records)

    df_avg = df.groupby("vertices").agg({"time": "mean"}).reset_index()

    df_sorted = df.sort_values(by=["vertices", "edges"])
    df_avg_sorted = df_avg.sort_values(by="vertices")

    plt.figure(figsize=(10, 6))
    plt.scatter(df_sorted["vertices"], df_sorted["time"], color='blue', label='Original Data Points')

    plt.plot(df_avg_sorted["vertices"], df_avg_sorted["time"], color='green', linestyle='-', label='Average Line')

    plt.xlabel("Number of Vertices")
    plt.ylabel(f"Time ({TIME_UNIT})")
    plt.title(f"Time vs Number of Vertices ({TIME_UNIT})")
    plt.grid(True)

    if TIME_UNIT in ["milliseconds", "microseconds"]:
        plt.ticklabel_format(axis='y', style='plain')
        plt.ylim(bottom=0)

    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE)
    print(f"Graph saved as {OUTPUT_FILE}")

