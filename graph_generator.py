import random
import os
import networkx as nx

MAX_N = 25
NUM_GRAPHS = 100

def generate_random_graphs(x: int, output_dir: str = "data_processed/examples"):
    os.makedirs(output_dir, exist_ok=True)

    for i in range(1, x + 1):
        n = random.randint(1, MAX_N)
        max_edges = n * (n - 1) // 2

        if max_edges == 0:
            edges = []
        else:
            # Start with a connected base graph (random spanning tree)
            G = nx.generators.random_tree(n)
            edges = set((min(u, v), max(u, v)) for u, v in G.edges())

            # Add more edges randomly up to max_edges
            possible_edges = list(set((i, j) for i in range(n) for j in range(i + 1, n)) - edges)
            remaining_edges = max_edges - len(edges)
            extra_edges = random.randint(0, remaining_edges)
            additional_edges = random.sample(possible_edges, extra_edges)
            edges.update(additional_edges)

        filename = os.path.join(output_dir, f"example{i}.hcp")
        with open(filename, 'w') as f:
            f.write(f"{n} {len(edges)}\n")
            for u, v in edges:
                f.write(f"{u} {v}\n")


generate_random_graphs(NUM_GRAPHS)
