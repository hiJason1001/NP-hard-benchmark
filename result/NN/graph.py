import re
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker
import pandas as pd
import numpy as np
import os
from scipy.special import gammaln

# ========================
# CONFIGURATION
# ========================
# NAMES = ["concorde_tsphcp", "concorde_FHCPCS", "concorde_ALL_hcp"]
NAMES = ["NN_ALL_hcp", "NN_FHCPCS", "NN_tsphcp"]
FOLDER = "result/NN"
TIME_UNIT = "seconds"
TITLE = f"NN Runtime vs Number of Vertices"
GRAPH_TYPE = "loglog"  # "loglog" or "powerfit"
DRAW_SEPARATE_LINES = True
FONT = 24
OUTPUT_FILE = f"result/NN/combined_plot.png"
# ========================

all_records = []
colors = plt.cm.tab10.colors

for idx, NAME in enumerate(NAMES):
    INPUT_FILE = f"{FOLDER}/{NAME}.txt"
    if not os.path.exists(INPUT_FILE):
        print(f"File {INPUT_FILE} not found. Skipping.")
        continue

    with open(INPUT_FILE, 'r') as file:
        data = file.read()

    lines = data.strip().splitlines()
    i = 0
    records = []

    while i < len(lines):
        if NAME == "concorde_FHCPCS" and i >= 160:
            break
        line = lines[i].strip()
        if line.endswith(".hcp") or line.endswith(".tsp"):
            if i + 3 >= len(lines):
                break

            vertices_edges_line = lines[i + 1].strip()
            time_line = lines[i + 2].strip().lower()
            result_line = lines[i + 3].strip().lower()

            vertices_match = re.search(r'Number of Vertices:\s*(\d+),\s*number of edges:\s*(\d+)', vertices_edges_line)
            if not vertices_match:
                i += 4
                continue

            vertices = int(vertices_match.group(1))
            edges = int(vertices_match.group(2))

            is_tle = "time limit exceeded" in time_line
            total_time = None

            if not is_tle:
                time_match = re.search(r'(\d+) hours (\d+) minutes (\d+) seconds (\d+) milliseconds (\d+) microseconds', time_line)
                if not time_match:
                    i += 4
                    continue

                h, m, s, ms, us = map(int, time_match.groups())
                total_seconds = h * 3600 + m * 60 + s + ms / 1000 + us / 1_000_000

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
                'time': total_time if total_time is not None else -1,  # Use -1 to flag TLEs
                'result': 'no' if is_tle else 'yes',
                'source': NAME
            })

            i += 4
        else:
            i += 1


    all_records.extend(records)

# ========================
# PLOT
# ========================
if not all_records:
    print("No records found! Please check the data format.")
else:
    df = pd.DataFrame(all_records)
    df_sorted = df.sort_values(by=["vertices", "edges"])
    plt.figure(figsize=(10, 6))

    tle_plotted = False 
    global_max_time = df[df["result"] == "yes"]["time"].max()
    tle_y = global_max_time * 1.5 if global_max_time > 0 else 1000

    if DRAW_SEPARATE_LINES:
        for idx, name in enumerate(NAMES):
            sub_df = df[df["source"] == name]
            color = colors[idx % len(colors)]

            ok_df = sub_df[sub_df["result"] == "yes"]
            tle_df = sub_df[sub_df["result"] == "no"]

            plt.scatter(ok_df["vertices"], ok_df["time"], color=color, s=5)

            if not tle_df.empty:
                plt.scatter(
                    tle_df["vertices"],
                    [tle_y] * len(tle_df),
                    color='red',
                    marker='x',
                    s=60,
                    label="Time Limit" if not tle_plotted else None
                )
                tle_plotted = True

            df_avg = ok_df.groupby("vertices").agg({"time": "mean"}).reset_index().sort_values(by="vertices")
            df_avg["time"] = df_avg["time"].replace(0, 1e-6)
            log_x = np.log(df_avg["vertices"])
            log_y = np.log(df_avg["time"])
            k, b = np.polyfit(log_x, log_y, 1)

            x_vals = np.linspace(df_avg["vertices"].min(), df_avg["vertices"].max(), 500)
            y_vals = np.exp(b) * x_vals**k
            label = f"{name} Fit ($n^{{{k:.2f}}}$)"
            # label = f"$n^{{{k:.2f}}}$"
            plt.plot(x_vals, y_vals, color=color, label=label)
    else:
        ok_df = df[df["result"] == "yes"]
        tle_df = df[df["result"] == "no"]

        plt.scatter(ok_df["vertices"], ok_df["time"], color='black', s=10, label='Solved')

        if not tle_df.empty:
            plt.scatter(
                tle_df["vertices"],
                [tle_y] * len(tle_df),
                color='red',
                marker='x',
                s=60,
                label="Time Limit"
            )

        df_avg = ok_df.groupby("vertices").agg({"time": "mean"}).reset_index()
        df_avg_sorted = df_avg.sort_values(by="vertices")
        df_avg_sorted["time"] = df_avg_sorted["time"].replace(0, 1e-6)
        log_x = np.log(df_avg_sorted["vertices"])
        log_y = np.log(df_avg_sorted["time"])
        k, b = np.polyfit(log_x, log_y, 1)

        x_vals = np.linspace(df_avg_sorted["vertices"].min(), df_avg_sorted["vertices"].max(), 500)
        y_vals = np.exp(b) * x_vals**k
        plt.plot(x_vals, y_vals, color='black', label=f'Power Law Fit ($n^{{{k:.2f}}}$)')

    if GRAPH_TYPE == "loglog":
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("Number of Vertices (log scale)", fontsize=FONT)
        plt.ylabel(f"Time ({TIME_UNIT}, log scale)", fontsize=FONT)
        ax = plt.gca()
        ax.set_ylim(top=tle_y * 1.2)
        ax.xaxis.set_major_formatter(matplotlib.ticker.LogFormatterMathtext())
        ax.yaxis.set_major_formatter(matplotlib.ticker.LogFormatterMathtext())
    else:
        plt.xlabel("Number of Vertices", fontsize=FONT)
        plt.ylabel(f"Time ({TIME_UNIT})", fontsize=FONT)
        plt.ylim(top=tle_y * 1.2)

    plt.title(TITLE, fontsize=FONT)
    plt.grid(True, which="both", ls="--")
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=FONT)
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE)
    print(f"Graph saved as {OUTPUT_FILE}")
