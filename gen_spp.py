"""Generate a random directed acyclic SPP test instance and write it to
file, in the format read by `read_spp` in tol.py (first line n, second
line "s t" 1-indexed, then n lines of n integers, -1 = no arc).

Usage:
    python gen_spp.py [fname] [n] [k] [lb] [ub]

    fname - name of the output file
    n     - number of vertices
    k     - number of arcs
    lb    - minimum arc cost, inclusive
    ub    - maximum arc cost, inclusive

Vertices are indexed 0..n-1 and arcs only ever go from a lower-indexed
vertex to a higher-indexed one, so the graph is always acyclic
regardless of which k arcs happen to be chosen (same trick gen_mstp.py
uses k random edges for connectivity, just applied to direction here
instead). The source is fixed at vertex 0 and the target at vertex
n-1 (1 and n in the 1-indexed output file). Costs are drawn uniformly
at random from [lb, ub], k arcs are chosen uniformly at random without
replacement from the n*(n-1)/2 possible forward arcs, and the whole
process is retried until vertex n-1 is actually reachable from vertex
0 (so the instance actually has a source-to-target path).
"""

import random
import sys


def gen_dag(n, k, lb, ub):
    """Return an n x n cost matrix (-1 = no arc) of a random DAG with k
    random arcs, each with an independently random cost in [lb, ub].
    Arcs only ever go from a lower-indexed vertex to a higher-indexed
    one, so the result is always acyclic. Not guaranteed to have a
    path from vertex 0 to vertex n-1; see `gen_reachable_dag`."""
    graph = []
    arcs = []

    for i in range(n):
        line = []
        for j in range(n):
            if j > i:
                arcs.append([i, j])
            line.append(-1)
        graph.append(line)

    for i in range(k):
        e = random.choice(arcs)
        arcs.remove(e)
        e_cost = random.choice(range(lb, ub + 1))
        graph[e[0]][e[1]] = e_cost

    return graph


def is_reachable(graph, s, t):
    """Return True iff `t` is reachable from `s` in `graph` (an n x n
    cost matrix, -1 = no arc), via a depth-first search from `s`."""
    n = len(graph)

    stack = []
    label = []

    for i in range(n):
        label.append(0)

    stack.append(s)
    while(len(stack) != 0):
        cur_v = stack.pop()
        label[cur_v] = 1
        for i in range(n):
            if (graph[cur_v][i] != -1) and (label[i] == 0):
                stack.append(i)

    return label[t] == 1


def gen_reachable_dag(n, k, lb, ub):
    """Return an n x n cost matrix (-1 = no arc) of a random DAG with k
    random arcs, each with an independently random cost in [lb, ub].
    The result is guaranteed to have a path from vertex 0 to vertex
    n-1."""
    s, t = 0, n - 1
    while(True):
        graph = gen_dag(n, k, lb, ub)
        if not is_reachable(graph, s, t):
            continue
        break

    return graph


def write_spp(fname, graph, s, t):
    """Writes generated graph and start/end vertices to file, in the
    format read by `read_spp` from `tol.py` (`s`/`t` are written
    1-indexed)."""
    n = len(graph)

    f = open(fname, 'w')
    f.write(str(n))
    f.write('\n')
    f.write(str(s + 1) + ' ' + str(t + 1))
    f.write('\n')
    for i in range(n):
        for j in range(n):
            if graph[i][j] != -1:
                f.write(' ')
            f.write(str(graph[i][j]))
            f.write(' ')
        f.write('\n')

    return


def main():
    if len(sys.argv) < 6:
        print("Usage: python gen_spp.py [fname] [n] [k] [lb] [ub]")
        return
    else:
        fname = sys.argv[1]
        n = int(sys.argv[2])
        k = int(sys.argv[3])
        lb = int(sys.argv[4])
        ub = int(sys.argv[5])

    graph = gen_reachable_dag(n, k, lb, ub)
    write_spp(fname, graph, 0, n - 1)

    return


if __name__ == "__main__":
    main()
