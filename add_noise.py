import os
import random

def read_edge_list(filepath):
    edges = set()
    with open(filepath, 'r') as f:
        header = f.readline()
        n, m = map(int, header.strip().split())
        for line in f:
            if not line.strip():
                continue
            u, v = map(int, line.strip().split())
            if u == v:
                continue
            edges.add((min(u, v), max(u, v)))
    return edges, n

def add_random_edges(edges, n, k, seed=42):
    random.seed(seed)
    edge_set = set(edges)
    max_possible_edges = n * (n - 1) // 2
    available_edges = max_possible_edges - len(edge_set)

    if k > available_edges:
        print(f"Warning: Only {available_edges} edges can be added, but requested {k}. Adding only {available_edges}.")
        k = available_edges

    added = 0
    attempts = 0
    max_attempts = 10 * k

    while added < k and attempts < max_attempts:
        u = random.randint(0, n - 1)
        v = random.randint(0, n - 1)
        attempts += 1
        if u == v:
            continue
        edge = (min(u, v), max(u, v))
        if edge in edge_set:
            continue
        edge_set.add(edge)
        added += 1

    if added < k:
        print(f"Stopped early: only added {added} edges out of {k} requested due to lack of new unique edges.")

    return edge_set


def write_edge_list(filepath, edges):
    n = 0
    for u, v in edges:
        n = max(n, u, v)
    n += 1
    m = len(edges)

    with open(filepath, 'w') as f:
        f.write(f"{n} {m}\n")
        for u, v in sorted(edges):
            f.write(f"{u} {v}\n")

def process_directory(input_dir, output_dir, epsilon=None, k_fixed=None, seed=42):
    os.makedirs(output_dir, exist_ok=True)
    files = [f for f in os.listdir(input_dir) if f.endswith('.hcp')]

    for file in files:
        input_path = os.path.join(input_dir, file)
        edges, n = read_edge_list(input_path)
        m = len(edges)

        if epsilon is not None:
            k = int(epsilon * len(edges))
        elif k_fixed is not None:
            k = k_fixed
        else:
            raise ValueError("Either epsilon or k_fixed must be specified.")

        noisy_edges = add_random_edges(edges, n, k, seed=seed)
        output_path = os.path.join(output_dir, file)
        write_edge_list(output_path, noisy_edges)
        print(f"Processed {file}: added {k} edges → {output_path}")



# ==== PARAMETERS ====

INPUT_DIR = "data_processed/tsphcp_processed"
OUTPUT_DIR = "data_noisy/tsphcp"
EPSILON = 0.05      # Add % of possible edges
K_FIXED = None      
SEED = 99

# ==== RUN ====

process_directory(
    input_dir=INPUT_DIR,
    output_dir=OUTPUT_DIR,
    epsilon=EPSILON,
    k_fixed=K_FIXED,
    seed=SEED
)