import sys
import os

def hcp_to_tsp(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    dim = 0
    edges = []
    edge_data_index = None

    # Handle variations in headers (like DIMENTION / EDGE_DATA_SELECTION)
    for i, line in enumerate(lines):
        line_clean = line.strip().upper()
        if line_clean.startswith("DIMENSION") or line_clean.startswith("DIMENTION"):
            try:
                dim = int(line.split(":")[1].strip())
            except:
                dim = int(line.split()[-1].strip())
        elif "EDGE_DATA_SECTION" in line_clean or "EDGE_DATA_SELECTION" in line_clean:
            edge_data_index = i + 1
            break

    if edge_data_index is None:
        raise ValueError(f"No valid EDGE_DATA_SECTION found in {input_file}")

    # Parse edge list
    for line in lines[edge_data_index:]:
        line = line.strip()
        if not line or line.upper() == "EOF":
            continue
        parts = line.split()
        if len(parts) != 2:
            print(f"Warning: Skipping malformed line in {input_file}: '{line}'")
            continue
        u, v = map(int, parts)
        edges.append((u - 1, v - 1))  # 0-indexed

    # Build matrix
    matrix = [[2] * dim for _ in range(dim)]
    for i in range(dim):
        matrix[i][i] = 0
    for u, v in edges:
        matrix[u][v] = matrix[v][u] = 1

    # Write TSP file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        f.write(f"NAME: converted\n")
        f.write(f"TYPE: TSP\n")
        f.write(f"DIMENSION: {dim}\n")
        f.write(f"EDGE_WEIGHT_TYPE: EXPLICIT\n")
        f.write(f"EDGE_WEIGHT_FORMAT: FULL_MATRIX\n")
        f.write(f"EDGE_WEIGHT_SECTION\n")
        for row in matrix:
            f.write(" ".join(map(str, row)) + "\n")
        f.write("EOF\n")
    


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 hcp_to_tsp.py input.hcp output.tsp")
    else:
        hcp_to_tsp(sys.argv[1], sys.argv[2])
