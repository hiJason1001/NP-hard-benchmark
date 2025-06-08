import os
import glob

INF = 10**9

def read_graph(filepath):
    with open(filepath, 'r') as f:
        n, m = map(int, f.readline().split())
        adj = [[INF]*n for _ in range(n)]
        for i in range(n):
            adj[i][i] = 0
        for _ in range(m):
            u,v = map(int, f.readline().split())
            # u -= 1
            # v -= 1
            adj[u][v] = 1
            adj[v][u] = 1
    return n, adj

def write_tsp(filename, n, adj):
    with open(filename, 'w') as f:
        f.write(f"NAME: {os.path.basename(filename)}\n")
        f.write("TYPE: TSP\n")
        f.write(f"DIMENSION: {n}\n")
        f.write("EDGE_WEIGHT_TYPE: EXPLICIT\n")
        f.write("EDGE_WEIGHT_FORMAT: FULL_MATRIX\n")
        f.write("EDGE_WEIGHT_SECTION\n")
        for i in range(n):
            row = []
            for j in range(n):
                w = adj[i][j]
                if w == INF:
                    w = 1000000 
                row.append(str(w))
            f.write(" ".join(row) + "\n")
        f.write("EOF\n")

def main(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    input_files = glob.glob(os.path.join(input_folder, "*.hcp"))
    for infile in input_files:
        n, adj = read_graph(infile)
        base = os.path.splitext(os.path.basename(infile))[0]
        outpath = os.path.join(output_folder, base + ".tsp")
        write_tsp(outpath, n, adj)
        print(f"Converted {infile} -> {outpath}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python convert_to_tsp.py <input_folder> <output_folder>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
