import sys

def hcp_to_tsp(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    dim = 0
    edges = []
    for line in lines:
        if line.startswith("DIMENSION"):
            dim = int(line.split(":")[1])
        elif line.strip() == "EDGE_DATA_SECTION":
            break

    edge_data_index = lines.index("EDGE_DATA_SECTION\n") + 1
    for line in lines[edge_data_index:]:
        line = line.strip()
        if not line or line == "EOF":
            continue
        parts = line.split()
        if len(parts) != 2:
            print(f"Warning: Skipping malformed line: '{line}'")
            continue
        u, v = map(int, parts)
        edges.append((u - 1, v - 1))  # 0-indexed

    matrix = [[2] * dim for _ in range(dim)]
    for i in range(dim):
        matrix[i][i] = 0

    for u, v in edges:
        matrix[u][v] = matrix[v][u] = 1

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
