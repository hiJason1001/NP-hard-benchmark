import re

FILE_WITH_RESULTS = "data_concorde/concorde_tsphcp_noisy_50.txt"
FILE_WITH_GRAPH_INFO = "data_concorde/concorde_tsphcp.txt"
OUTPUT_FILE = "data_concorde/results_with_vertices.txt"

graph_info = {}
with open(FILE_WITH_GRAPH_INFO, 'r') as f:
    lines = f.readlines()

i = 0
while i < len(lines):
    line = lines[i].strip()
    if line.endswith(".tsp") or line.endswith(".hcp"):
        filename = line.strip().lower()
        if i + 1 < len(lines):
            ve_line = lines[i + 1].strip()
            if "number of vertices" in ve_line.lower():
                graph_info[filename] = ve_line
        i += 2
    else:
        i += 1

with open(FILE_WITH_RESULTS, 'r') as fin, open(OUTPUT_FILE, 'w') as fout:
    lines = fin.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        fout.write(line + "\n")

        if line.endswith(".tsp") or line.endswith(".hcp"):
            filename = line.strip().lower().replace(".tsp", ".hcp")
            ve_line = graph_info.get(filename, "Number of Vertices: ?, number of edges: ?")
            fout.write(ve_line + "\n")
            if i + 1 < len(lines):
                fout.write(lines[i + 1])
            if i + 2 < len(lines):
                fout.write(lines[i + 2])
            i += 3
        else:
            i += 1
