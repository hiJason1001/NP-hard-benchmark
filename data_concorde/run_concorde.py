import os
import subprocess
import time

CONCORDE_PATH = "./../../concorde/TSP/concorde"
TSP_DIR = "tsphcp_noisy"
HCP_DIR = "../data_processed/tsphcp_processed"
OUTPUT_LOG = "concorde_results.txt"
# TIME_LIMIT = 600  # second
TIME_LIMIT = 10 # second

def get_vertices_and_edges(hcp_file):
    with open(hcp_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    return int(parts[0]), int(parts[1])
    return None, None

def format_time(duration):
    ms = int((duration - int(duration)) * 1000)
    us = int((duration * 1_000_000) % 1000)
    duration = int(duration)
    h = duration // 3600
    m = (duration % 3600) // 60
    s = duration % 60
    return f"{h} hours {m} minutes {s} seconds {ms} milliseconds {us} microseconds"

def run_all_concorde():
    processed = set()
    if os.path.exists(OUTPUT_LOG):
        with open(OUTPUT_LOG, 'r') as f:
            for line in f:
                if line.strip().endswith(".tsp"):
                    processed.add(line.strip())

    tsp_files = sorted([
        os.path.join(root, file)
        for root, _, files in os.walk(TSP_DIR)
        for file in files if file.endswith(".tsp")
    ])

    with open(OUTPUT_LOG, 'a') as out:
        for tsp_path in tsp_files:
            base = os.path.basename(tsp_path)
            if base in processed:
                continue

            out.write(f"{base}\n")

            try:
                start = time.time()
                subprocess.run([CONCORDE_PATH, tsp_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=TIME_LIMIT)
                duration = time.time() - start
                out.write(f"{format_time(duration)}\nYes\n")
            except subprocess.TimeoutExpired:
                out.write("Time Limit Exceeded\nNo\n")

            os.system("rm -f *.mas *.pul *.sav *.sol *.res *.[0-9][0-9][0-9]")
            out.flush()



def run_fhcpcs_specific():
    processed_bases = set()

    if os.path.exists(OUTPUT_LOG):
        with open(OUTPUT_LOG, 'r') as f:
            for line in f:
                line = line.strip()
                if line.endswith(".hcp"):
                    base = line[:-4]
                    processed_bases.add(base)

    with open(OUTPUT_LOG, 'a') as out:
        for i in range(1001):
            base = f"graph{i}"
            if base in processed_bases:
                print(f"Skipping {base} as it is already processed.")
                continue
            
            tsp_path = os.path.join(TSP_DIR, f"{base}.tsp")
            hcp_path = os.path.join(HCP_DIR, f"{base}.hcp")

            if not os.path.isfile(tsp_path):
                print(f"{tsp_path} does not exist, skipping.")
                continue

            out.write(f"{base}.hcp\n")

            v, e = get_vertices_and_edges(hcp_path)
            if v is not None and e is not None:
                out.write(f"Number of Vertices: {v}, number of edges: {e}\n")
            else:
                out.write("Number of Vertices: ?, number of edges: ?\n")

            try:
                start = time.time()
                subprocess.run([CONCORDE_PATH, tsp_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=TIME_LIMIT)
                duration = time.time() - start
                out.write(f"{format_time(duration)}\nYes\n")
                os.system("rm -f *.mas *.pul *.sav *.sol *.res *.[0-9][0-9][0-9] *graph*")
            except subprocess.TimeoutExpired:
                out.write("Time Limit Exceeded\nNo\n")
                os.system("rm -f *.mas *.pul *.sav *.sol *.res *.[0-9][0-9][0-9] *graph*")

            out.flush()

if __name__ == "__main__":
    run_all_concorde()
    # run_fhcpcs_specific()