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

# n1 = "concorde"
n4 = "hybridHam"
n5 = "LKH"
n6 = "NN"

FILENAME = "FHCPCS"

FILES = [
    # (f"result/{n1}/{n1}_{FILENAME}.txt", f"{n1}_{FILENAME}", "power"),
    (f"result/{n4}/{n4}_{FILENAME}.txt", f"{n4}_{FILENAME}", "power"),
    (f"result/{n5}/{n5}_{FILENAME}.txt", f"{n5}_{FILENAME}", "power"),
    (f"result/{n6}/{n6}_{FILENAME}.txt", f"{n6}_{FILENAME}", "power"),
]


TIME_UNIT = "seconds"
TITLE = "All Algorithm Runtime vs Number of Vertices on FHCPCS Dataset"
GRAPH_TYPE = "loglog"  # "loglog" or "linear"
DRAW_SEPARATE_LINES = True 
FIT_TYPE = "power"  # "power" or "exponential"
FONT = 16
PLOT_TLE = True  # Set False to hide TLE points
OUTPUT_FILE = "combined_plot_power.png"
# ========================

all_records = []
colors = plt.cm.tab10.colors

for idx, (INPUT_FILE, NAME, FIT_TYPE) in enumerate(FILES):
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
                'time': total_time if total_time is not None else -1,
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
        for idx, (_, name, fit_type) in enumerate(FILES):
            sub_df = df[df["source"] == name]
            color = colors[idx % len(colors)]

            ok_df = sub_df[sub_df["result"] == "yes"]
            tle_df = sub_df[sub_df["result"] == "no"]

            plt.scatter(ok_df["vertices"], ok_df["time"], color=color, s=1)

            if PLOT_TLE and not tle_df.empty:
                plt.scatter(
                    tle_df["vertices"],
                    [tle_y] * len(tle_df),
                    color='red',
                    marker='x',
                    s=60,
                    label="Time Limit by LKH" if not tle_plotted else None
                )
                tle_plotted = True


            df_avg = ok_df.groupby("vertices").agg({"time": "mean"}).reset_index().sort_values(by="vertices")
            df_avg["time"] = df_avg["time"].replace(0, 1e-6)
            x_vals = np.linspace(df_avg["vertices"].min(), df_avg["vertices"].max(), 500)

            clean_name = name.replace("_examples", "")
            if fit_type  == "power":
                log_x = np.log(df_avg["vertices"])
                log_y = np.log(df_avg["time"])
                k, b = np.polyfit(log_x, log_y, 1)
                y_vals = np.exp(b) * x_vals ** k

                label = f"{clean_name} ($n^{{{k:.2f}}}$)"

            elif fit_type == "exponential":
                x = df_avg["vertices"].values
                y = df_avg["time"].values
                log_y = np.log(y)
                m, c = np.polyfit(x, log_y, 1)
                a = np.exp(c)
                b = np.exp(m)

                y_vals = a * (b ** x_vals)
                base = np.log2(b)
                label = f"{clean_name} ($2^{{{base:.2f}n}}$)"

            else:
                raise ValueError("FIT_TYPE must be 'power' or 'exponential'")

            plt.plot(x_vals, y_vals, color=color, label=label, linewidth=2)
    else:
        ok_df = df[df["result"] == "yes"]
        tle_df = df[df["result"] == "no"]

        plt.scatter(ok_df["vertices"], ok_df["time"], color='black', s=10, label='Solved')

        if PLOT_TLE and not tle_df.empty:
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
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(fontsize=12, ncol=2)
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE)
    print(f"Graph saved as {OUTPUT_FILE}")
