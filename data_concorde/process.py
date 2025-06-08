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

def simple_hcp_to_tsp(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    try:
        dim, num_edges = map(int, lines[0].split())
    except Exception as e:
        raise ValueError(f"Invalid format in {input_file} (expected 'n m' on first line): {e}")

    edges = []
    for line in lines[1:]:
        try:
            u, v = map(int, line.split())
            edges.append((u, v))
        except ValueError:
            print(f"Warning: Skipping malformed edge line in {input_file}: '{line}'")

    matrix = [[2] * dim for _ in range(dim)]
    for i in range(dim):
        matrix[i][i] = 0
    for u, v in edges:
        if 0 <= u < dim and 0 <= v < dim:
            matrix[u][v] = matrix[v][u] = 1
        else:
            print(f"Warning: Skipping out-of-bounds edge ({u}, {v}) in {input_file}")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        f.write("NAME: converted\n")
        f.write("TYPE: TSP\n")
        f.write(f"DIMENSION: {dim}\n")
        f.write("EDGE_WEIGHT_TYPE: EXPLICIT\n")
        f.write("EDGE_WEIGHT_FORMAT: FULL_MATRIX\n")
        f.write("EDGE_WEIGHT_SECTION\n")
        for row in matrix:
            f.write(" ".join(map(str, row)) + "\n")
        f.write("EOF\n")


def process_all_hcp_files(input_root='data_processed/examples', output_root='data_concorde/examples'):
    for dirpath, _, filenames in os.walk(input_root):
        for filename in filenames:
            if filename.lower().endswith(".hcp"):
                input_path = os.path.join(dirpath, filename)
                relative_path = os.path.relpath(input_path, input_root)
                output_path = os.path.join(output_root, os.path.splitext(relative_path)[0] + ".tsp")
                try:
                    # hcp_to_tsp(input_path, output_path)
                    simple_hcp_to_tsp(input_path, output_path)
                    print(f"Converted: {input_path} -> {output_path}")
                except Exception as e:
                    print(f"Failed to convert {input_path}: {e}")

if __name__ == "__main__":
    process_all_hcp_files()
