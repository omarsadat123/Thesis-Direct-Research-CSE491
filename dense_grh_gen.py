#!/usr/bin/env python3
"""
dense_grh_gen.py  — generate a .grh adjacency-list file (0-based labels)

Format:
- Line i (0 ≤ i < n) lists the neighbors of vertex i, space-separated
- No self-loops, undirected (symmetric)
- No trailing spaces on lines
- No extra newline after the *last* line

Modes:
  1) --complete                  → K_n (perfectly dense)
  2) --density P (0..1) [--seed] → Erdős–Rényi G(n, P) (symmetric)
"""

import argparse, random, sys

def write_grh(n, adj_lists, out_path):
    # Build lines first, then join with '\n' to avoid a final extra newline.
    # Ensure no trailing spaces by joining exactly the neighbor IDs.
    lines = [" ".join(map(str, nbrs)) for nbrs in adj_lists]
    data = "\n".join(lines)
    with open(out_path, "w", newline="") as f:
        f.write(data)

def gen_complete(n):
    # For K_n, neighbors of u are all vertices except u.
    # Build once with simple list arithmetic; already sorted and symmetric.
    adj = []
    for u in range(n):
        # Split to avoid an O(n) temporary if n is huge? This is fine and clear:
        nbrs = list(range(u)) + list(range(u + 1, n))
        adj.append(nbrs)
    return adj

def gen_er_symmetric(n, p, seed=None):
    # Undirected ER: generate edges in the upper triangle, then mirror.
    if seed is not None:
        random.seed(seed)
    adj = [[] for _ in range(n)]
    # Upper-triangle edge sampling
    for u in range(n):
        # v starts at u+1 to avoid self-loops and duplicates
        for v in range(u + 1, n):
            if random.random() < p:
                adj[u].append(v)
                adj[v].append(u)
    # Ensure neighbors sorted (they already are by construction)
    # but keep this in case someone edits to use another sampling strategy.
    for u in range(n):
        adj[u].sort()
    return adj

def main():
    ap = argparse.ArgumentParser(description="Generate a dense .grh file (0-based).")
    ap.add_argument("--n", type=int, required=True, help="number of vertices (>=1)")
    ap.add_argument("--out", type=str, required=True, help="output .grh path")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--complete", action="store_true", help="generate complete graph K_n")
    mode.add_argument("--density", type=float, help="ER density p in [0,1] (undirected, symmetric)")
    ap.add_argument("--seed", type=int, help="RNG seed (only used with --density)")
    args = ap.parse_args()

    if args.n < 1:
        sys.exit("n must be >= 1")

    if args.complete:
        adj = gen_complete(args.n)
    else:
        p = 1.0 if args.density is None else args.density
        if not (0.0 <= p <= 1.0):
            sys.exit("--density must be in [0,1]")
        adj = gen_er_symmetric(args.n, p, args.seed)

    # Final integrity checks (cheap)
    n = args.n
    # 0-based labels, in-range, and no self-loops:
    for u in range(n):
        for v in adj[u]:
            if not (0 <= v < n):
                sys.exit(f"Neighbor out of range at row {u}: {v}")
            if v == u:
                sys.exit(f"Self-loop at row {u}")
    # Symmetry (sampled to keep it cheap for huge n):
    SAMPLE = min(n, 100)
    step = max(1, n // SAMPLE)
    for u in range(0, n, step):
        s = set(adj[u])
        for v in s:
            if u not in set(adj[v]):
                sys.exit(f"Asymmetry: {u} in adj[{v}] is missing")
    # Sorted lines, no duplicates:
    for u in range(n):
        if adj[u] != sorted(adj[u]):
            sys.exit(f"Neighbors not sorted at row {u}")
        if len(adj[u]) != len(set(adj[u])):
            sys.exit(f"Duplicate neighbor at row {u}")

    write_grh(n, adj, args.out)

if __name__ == "__main__":
    main()
