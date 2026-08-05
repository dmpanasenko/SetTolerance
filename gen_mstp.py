"""Generate a random connected MSTP test instance and write it to
file, in the format read by `read_graph` in mst_tolerance.py
(first line n, then n lines of n integers, -1 = no edge).

Usage:
    python gen_mstp.py [fname] [n] [k] [lb] [ub]

    fname - name of the output file
    n     - number of vertices
    k     - number of edges
    lb    - minimum edge cost, inclusive
    ub    - maximum edge cost, inclusive

Costs are drawn uniformly at random from [lb, ub], k edges are chosen
uniformly at random without replacement from the n*(n-1)/2 possible
edges, and the whole process is retried until the resulting graph is
connected (so that it actually is a valid MSTP instance).
"""

import random
import sys


def gen_graph(n, k, lb, ub):
    """Return an n x n cost matrix (-1 = no edge) with k random edges,
    each with an independently random cost in [lb, ub]. The result is
    not guaranteed to be connected; see `is_connected`."""
    graph = []
    edges = []
    chosen_edges = []

    for i in range(n):
        line = []
        for j in range(n):
            if j > i:
                edges.append([i, j])
            line.append(-1)
        graph.append(line)

    for i in range(k):
        e = random.choice(edges)
        edges.remove(e)
        e_cost = random.choice(range(lb, ub+1))
        graph[e[0]][e[1]] = e_cost
        graph[e[1]][e[0]] = e_cost
        chosen_edges.append(e)

    return graph


def is_connected(graph):
    """Return True iff `graph` (an n x n cost matrix, -1 = no edge) is
    connected, via a depth-first search from vertex 0."""
    n = len(graph)

    stack = []
    label = []

    for i in range(n):
        label.append(0)

    start_v = 0

    stack.append(start_v)
    while(len(stack) != 0):
        cur_v = stack.pop()
        label[cur_v] = 1
        for i in range(n):
            if (graph[i][cur_v] != -1) and (label[i] == 0):
                stack.append(i)

    if sum(label) == n:
        return True
    else:
        return False


def gen_mstp(n, k, lb, ub):
    """Return an n x n cost matrix (-1 = no edge) with k random edges,
    each with an independently random cost in [lb, ub]. The result is
    guaranteed to be connected."""
    while(True):
        graph = gen_graph(n, k, lb, ub)
        if not is_connected(graph):
            continue
        break
    
    return graph


def write_graph(fname, graph):
    """Writes generated graph to file in format, 
    compatable with `read_graph` from `mst_tolerance.py`."""
    n = len(graph)

    f = open(fname, 'w')
    f.write(str(n))
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
        print("Usage: python gen_mstp.py [fname] [n] [k] [lb] [ub]")
        return
    else:
        fname = sys.argv[1]
        n = int(sys.argv[2])
        k = int(sys.argv[3])
        lb = int(sys.argv[4])
        ub = int(sys.argv[5])

    graph = gen_mstp(n, k, lb, ub)
    write_graph(fname, graph)

    return


if __name__ == "__main__":
    main()
