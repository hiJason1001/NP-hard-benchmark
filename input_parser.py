import os
import re


input_folder = "data_raw/ALL_hcp"
output_folder = "data_processed/ALL_hcp_processed"

os.makedirs(output_folder, exist_ok=True)

def parse_input_file(input_path: str, output_path: str):
    n = None
    edges = []
    reading_edges = False

    with open(input_path, 'r') as file:
        for line in file:
            stripped = line.strip()

            if not stripped:
                continue

            if not reading_edges:
                # DIMENSION: number of vertices
                match = re.search(r"DIMENSION\s*:\s*(\d+)", stripped)
                if match:
                    n = int(match.group(1))

                # If the line is two integers, assume edge section has started
                elif re.match(r"\d+\s+\d+", stripped):
                    reading_edges = True

            if reading_edges:
                if stripped == "-1":
                    break
                parts = stripped.split()
                if len(parts) == 2:
                    u, v = map(int, parts)
                    edges.append((u - 1, v - 1))  # Convert to 0-indexed

    if n is None:
        vertex_set = set()
        for u, v in edges:
            vertex_set.update([u, v])
        n = max(vertex_set) + 1

    m = len(edges)

    with open(output_path, 'w') as out:
        out.write(f"{n} {m}\n")
        for u, v in edges:
            out.write(f"{u} {v}\n")


for filename in os.listdir(input_folder):
    if filename.endswith(".hcp"):
        input_path = os.path.join(input_folder, filename)
        output_path = os.path.join(output_folder, filename)

        parse_input_file(input_path, output_path)
        print(f"Processed {filename}")