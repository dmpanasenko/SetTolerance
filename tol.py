"""Computation of regular set tolerances for combinatorial sum problems.

Reference implementation for:

G. Jager, D. Panasenko, "Computation of Set Tolerances with
Applications to the Minimum Spanning Tree Problem", arXiv:2605.10542.

Terminology (see the paper for full definitions):
    EUL  Exact Upper LP     - Section 3,  Theorem 7,  constraints (1)-(2)
    ELL  Exact Lower LP     - Section 4,  Theorem 10, constraint  (3)
    ILL  Include Lower LP   - Section 4,  Theorem 13, constraint  (5)
    TLL  Tolerance Lower LP - Section 4,  Theorem 13, constraints (6)-(7)

A "problem_func" callback represents the combinatorial sum problem (CSP)
being analyzed. It must have the signature
    problem_func(instance, given, banned) -> (solution, cost)
where `given` is a list of elements forced to belong to the solution,
`banned` is a list of elements forced to be excluded, and `cost` is
float('inf') if no feasible solution respects those constraints. This
module ships three implementations: `mstp` for the Minimum Spanning 
Tree Problem, `lap` for the Linear Assignment Problem and `spp` for the
Shortest Path Problem.
"""


from itertools import chain, combinations
from time import time

import numpy as np
from scipy.optimize import linear_sum_assignment as lap_solver
from scipy.optimize import linprog


# --------------------------------------------------------------------------
# Generic helpers
# --------------------------------------------------------------------------

def powerset(iterable):
    """Return every subset of `iterable`, as tuples, from size 0 up to n."""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def _solve_tolerance_lp(E, A, b, Aeq=None, beq=None, timer=None):
    """Solve max sum(alpha_i) s.t. A @ alpha <= b, Aeq @ alpha == beq, alpha >= 0.

    This is the common LP shape behind EUL/ELL/ILL/TLL; only the
    constraint rows differ between them; each of `eul`/`ell`/`ill`/`tll`
    (and their recursive variants) builds A/b/Aeq/beq and delegates here.

    Parameters:
        E: the ground-set subset the tolerance is being computed for -
            used only to label the returned per-element increases.
        A,b: inequality constraint rows/right-hand sides (may be empty).
        Aeq,beq: equality constraint rows/right-hand sides (may be empty).
        timer: optional 1-element list; timer[0] accumulates LP solve time.

    Returns:
        (increases, tolerance): `increases` maps each element of E to its
                                individual (rounded) contribution to the tolerance, `tolerance` is
                                the total (rounded) LP objective value.
                                None, float('inf')) if the LP is unbounded, meaning the set
                                tolerance is infinite.
                                (None, -1) for any other non-optimal solver status.
    """
    k = len(E)
    c = [-1] * k
    bounds = [(0, None)] * k

    if timer is not None:
        start_t = time()
    res = linprog(c, A_ub=A or None, b_ub=b or None,
                  A_eq=Aeq or None, b_eq=beq or None, bounds=bounds)
    if timer is not None:
        timer[0] += time() - start_t

    if res.status == 3:
        return None, float('inf')
    if res.status != 0:
        return None, -1

    increases = {tuple(E[i]): int(res.x[i]) for i in range(k)}
    return increases, int(-res.fun)


# --------------------------------------------------------------------------
# Instance I/O
# --------------------------------------------------------------------------

def read_mstp(fname):
    """Read an n x n adjacency/cost matrix (-1 = no edge) from a text file.

    File format: first line is n, followed by n lines of 
    n whitespace-separated integers each.
    """
    with open(fname) as f:
        n = int(f.readline())
        return [[int(u) for u in f.readline().split()] for _ in range(n)]


def print_mstp(graph):
    """Print an n x n matrix as whitespace-separated rows."""
    for row in graph:
        print(*row)
    return


def read_lap(fname):
    """Read an n x n assignment cost matrix from a text file.
    
    File format: first line is n, followed by n lines of 
    n whitespace-separated integers each."""
    with open(fname) as f:
        n = int(f.readline())
        return [[int(u) for u in f.readline().split()] for _ in range(n)]


def print_lap(matrix):
    """Print an n x n matrix as whitespace-separated rows."""
    for row in matrix:
        print(*row)
    return

def read_spp(fname):
    """Read an n x n adjacency/cost matrix (-1 = no edge) from a text file,
    and read start and end vertices.

    File format: first line is n, 
    second line is s t (whitespace-separated numbers of start and end vertices,
    numbering starts with 1)
    followed by n lines of n whitespace-separated integers each."""
    with open(fname) as f:
        n = int(f.readline())
        s, t = [int(v) - 1 for v in f.readline().split()]
        graph = [[int(u) for u in f.readline().split()] for _ in range(n)]
        return (graph, s, t)


def print_spp(graph, s, t):
    """Print start and end vertices, and print an n x n matrix 
    as whitespace-separated rows."""
    print("Start vertex:", s)
    print("End vertex:", t)
    for row in graph:
        print(*row)
    return


# --------------------------------------------------------------------------
# problem_func implementation #1: Minimum Spanning Tree Problem 
# --------------------------------------------------------------------------

def find_set(parent, v):
    """Union-find: return the representative of v's set (with path compression)."""
    if parent[v] != v:
        parent[v] = find_set(parent, parent[v])
    return parent[v]


def union_set(parent, rank, u, v):
    """Union-find: merge the sets containing u and v (union by rank)."""
    ru, rv = find_set(parent, u), find_set(parent, v)
    if ru == rv:
        return
    if rank[ru] < rank[rv]:
        parent[ru] = rv
    elif rank[ru] > rank[rv]:
        parent[rv] = ru
    else:
        parent[rv] = ru
        rank[ru] += 1
    return


def mstp(graph, given, banned):
    """MSTP problem_func: minimum spanning tree with forced-in/forced-out edges.

    Implements Kruskal's algorithm on `graph` (an n x n cost matrix with -1
    marking non-edges), restricted to spanning trees that contain every
    edge in `given` and none of the edges in `banned`.

    Trick used to enforce `given`/`banned` without a separate feasibility
    pass: `given` edges are temporarily assigned cost -1 so Kruskal's
    greedy sort picks them first, and `banned` edges are assigned a cost
    larger than any possible tree so they are only ever picked if the
    graph leaves no alternative (in which case the later feasibility
    checks below catch it and return infeasible). The -1 given-edge cost
    is added back to the running total (+1 to cancel it, then + the real
    edge cost) once the tree is complete, so the returned cost is correct.

    Returns:
        (tree, cost): `tree` is a list of [i, j] edges, `cost` is the tree
                      weight. Returns (None, float('inf')) if no spanning tree respects
                      `given`/`banned`.
    """
    n = len(graph)
    edges = [[[i, j], -1] for i, j in given]

    big_cost = 0
    for i in range(n):
        for j in range(i + 1, n):
            if graph[i][j] != -1:
                big_cost += graph[i][j]
                if [i, j] in given or [i, j] in banned:
                    continue
                edges.append([[i, j], graph[i][j]])

    for i, j in banned:
        edges.append([[i, j], big_cost + 1])

    edges.sort(key=lambda edge: edge[1])

    parent = list(range(n))
    rank = [0] * n

    tree = []
    tree_cost = 0
    for (u, v), w in edges:
        ru, rv = find_set(parent, u), find_set(parent, v)
        if ru != rv:
            union_set(parent, rank, ru, rv)
            tree.append([u, v])
            tree_cost += w
        if len(tree) == n - 1:
            break

    if len(tree) != n - 1:
        return None, float('inf')

    for i, j in banned:
        if [i, j] in tree:
            return None, float('inf')

    for i, j in given:
        if [i, j] not in tree:
            return None, float('inf')
        tree_cost += graph[i][j] + 1 

    return tree, tree_cost


# --------------------------------------------------------------------------
# problem_func implementation #2: Linear Assignment Problem
# --------------------------------------------------------------------------

def lap(matr, given, banned):
    """LAP problem_func: min-cost perfect matching with forced-in/out arcs.

    Solves the linear assignment problem on cost matrix `matr`, restricted
    to matchings that contain every (row, col) arc in `given` and none of
    the arcs in `banned`, by fixing the `given` rows/cols and running the
    Hungarian algorithm (`scipy.optimize.linear_sum_assignment`) on the
    remaining reduced matrix with `banned` arcs blocked out (cost = inf).

    Returns:
        (arcs, cost):   `arcs` is a list of (row, col) tuples, `cost` is the
                        total matching cost. Returns (None, float('inf')) if `given`
                        conflicts with itself or with `banned`, or if the reduced problem
                        is infeasible.
    """
    matrix = np.asarray(matr, dtype=float)
    n = matrix.shape[0]

    rows_used, cols_used = set(), set()
    fixed_cost = 0.0

    for i, j in given:
        if i in rows_used or j in cols_used or (i, j) in banned:
            return None, float('inf')
        rows_used.add(i)
        cols_used.add(j)
        fixed_cost += matrix[i, j]

    free_rows = [i for i in range(n) if i not in rows_used]
    free_cols = [j for j in range(n) if j not in cols_used]

    if not free_rows and free_cols:
        return list(given), fixed_cost

    reduced = matrix[np.ix_(free_rows, free_cols)].copy()

    for i, j in banned:
        if i in free_rows and j in free_cols:
            reduced[free_rows.index(i), free_cols.index(j)] = float('inf')

    try:
        r_idx, c_idx = lap_solver(reduced)
    except ValueError:
        return None, float('inf')

    arcs = list(given)
    reduced_cost = 0.0
    for r, c in zip(r_idx, c_idx):
        if reduced[r, c] == float('inf'):
            return None, float('inf')
        i, j = free_rows[r], free_cols[c]
        arcs.append((i, j))
        reduced_cost += matrix[i, j]

    return arcs, fixed_cost + reduced_cost


# --------------------------------------------------------------------------
# problem_func implementation #3: Shortest Path Problem
# --------------------------------------------------------------------------

def spp(instance, given, banned):
    """SPP problem_func: min-cost path from `s` to `t` with forced-in/out arcs.

    Solves the shortest path problem given by `instance`, which is a tuple
    (`graph`, `s`, `t`) where `graph` is an adjacency matrix of a directed
    acyclic graph (-1 = no arc, all present costs non-negative), `s` is
    the start vertex and `t` is the end vertex.

    Solutions are restricted to paths that contain every (u, v) arc in
    `given` and none of the arcs in `banned`. This is achieved by pricing
    each arc in `given` at -M-1 (M = sum of every arc's real cost in the
    graph) and excluding every arc in `banned` from the search entirely.
    The resulting instance is solved with Bellman-Ford, and the path found
    is checked to make sure it actually contains every arc from `given`
    (the discount only guarantees the path uses as many `given` arcs as
    any s-t path in the graph can - if that number is less than all of
    them, no path satisfies the constraints).

    Returns:
        (path, cost): `path` is a list of (u, v) arc tuples forming a
        path from s to t, `cost` is its total real cost. Returns
        (None, float('inf')) if no path from s to t respects
        `given`/`banned`.
    """
    graph, s, t = instance
    n = len(graph)

    given_arcs = set(map(tuple, given))
    banned_edges = set(map(tuple, banned))

    for u, v in given_arcs:
        if graph[u][v] == -1:
            return None, float('inf')

    big_M = sum(graph[i][j] for i in range(n) for j in range(n) if graph[i][j] != -1)
    discount = big_M + 1

    def weight(i, j):
        return -discount if (i, j) in given_arcs else graph[i][j]

    dist = [float('inf')] * n
    prev = [-1] * n
    dist[s] = 0

    for _ in range(n - 1):
        updated = False
        for i in range(n):
            for j in range(n):
                if graph[i][j] == -1 or (i, j) in banned_edges:
                    continue
                if dist[i] != float('inf') and dist[i] + weight(i, j) < dist[j]:
                    dist[j] = dist[i] + weight(i, j)
                    prev[j] = i
                    updated = True
        if not updated:
            break

    if dist[t] == float('inf'):
        return None, float('inf')

    path = []
    curr = t
    while prev[curr] != -1:
        path.append((prev[curr], curr))
        curr = prev[curr]
    path.reverse()

    if not given_arcs.issubset(path):
        return None, float('inf') 

    cost = sum(graph[i][j] for i, j in path)
    return path, cost


# --------------------------------------------------------------------------
# Set upper tolerance: EUL 
# --------------------------------------------------------------------------

def eul(instance, E, problem_func, S_star=None, c_star=None, timer=None):
    """Exact Upper LP: compute the regular set upper tolerance u'_P(E).

    `S_star`/`c_star` are the base instance's optimal solution/value;
    computed automatically if omitted, but pass them in explicitly when
    calling `eul` repeatedly on the same instance to avoid recomputing it
    every time.

    Returns:
        (increases, tolerance): see `_solve_tolerance_lp`; or (None, float('inf'))
                                if the base instance itself is infeasible.
    """
    k = len(E)
    if S_star is None or c_star is None:
        S_star, c_star = problem_func(instance, [], [])
    if c_star == float('inf'):
        return None, float('inf')

    E_prime_inds = [i for i in range(k) if E[i] in S_star]

    A, b, Aeq, beq = [], [], [], []
    for F in powerset(E):
        F = list(F)
        not_F = [e for e in E if e not in F]
        F_inds = [i for i in range(k) if E[i] in F]

        _, f_c = problem_func(instance, F, not_F)
        if f_c == float('inf'):
            continue

        line = [0] * k
        for i in E_prime_inds:
            line[i] += 1
        for i in F_inds:
            line[i] -= 1

        if f_c == c_star:
            Aeq.append(line)
            beq.append(0)
        else:
            A.append(line)
            b.append(f_c - c_star)

    return _solve_tolerance_lp(E, A, b, Aeq, beq, timer)


def compute_all_tolerances_eul(instance, elems, problem_func, timer=None):
    """Compute u'_P(E) (via EUL) for every non-empty subset E of `elems`."""
    start_t = time()
    S_star, c_star = problem_func(instance, [], [])
    for E in powerset(elems):
        if E:
            eul(instance, list(E), problem_func, S_star, c_star, timer)
    print(f"EUL run time (no rec): {time() - start_t:.3f}")
    if timer is not None:
        print(f"LP solving time: {timer[0]:.3f}")
    return


# --------------------------------------------------------------------------
# Set lower tolerance, approach 1: ELL 
# --------------------------------------------------------------------------

def ell(instance, E, problem_func, S_star=None, c_star=None, timer=None):
    """Exact Lower LP: compute the regular set lower tolerance l'_P(E)."""
    k = len(E)
    if S_star is None or c_star is None:
        S_star, c_star = problem_func(instance, [], [])
    if c_star == float('inf'):
        return None, float('inf')

    A, b = [], []
    for F in powerset(E):
        F = list(F)
        not_F = [e for e in E if e not in F]
        F_inds = [i for i in range(k) if E[i] in F]

        _, f_c = problem_func(instance, F, not_F)
        if f_c == float('inf'):
            continue

        line = [0] * k
        for i in F_inds:
            line[i] += 1
        A.append(line)
        b.append(f_c - c_star)

    return _solve_tolerance_lp(E, A, b, timer=timer)


def compute_all_tolerances_ell(instance, elems, problem_func, timer=None):
    """Compute l'_P(E) (via ELL) for every non-empty subset E of `elems`."""
    start_t = time()
    S_star, c_star = problem_func(instance, [], [])
    for E in powerset(elems):
        if E:
            ell(instance, list(E), problem_func, S_star, c_star, timer)
    end_t = time() - start_t
    print(f"ELL run time (no rec): {end_t:.3f}")
    if timer is not None:
        print(f"LP solving time: {timer[0]:.3f}, supporting time: {end_t - timer[0]:.3f}")
    return


# --------------------------------------------------------------------------
# Set lower tolerance, approach 2: ILL
# --------------------------------------------------------------------------

def ill(instance, E, problem_func, S_star=None, c_star=None, timer=None):
    """Include Lower LP: compute l'_P(E) using only fc(D+(F)) constraints."""
    k = len(E)
    if S_star is None or c_star is None:
        S_star, c_star = problem_func(instance, [], [])
    if c_star == float('inf'):
        return None, float('inf')

    A, b = [], []
    for F in powerset(E):
        F = list(F)
        F_inds = [i for i in range(k) if E[i] in F]

        _, f_c = problem_func(instance, F, [])
        if f_c == float('inf'):
            continue

        line = [0] * k
        for i in F_inds:
            line[i] += 1
        A.append(line)
        b.append(f_c - c_star)

    return _solve_tolerance_lp(E, A, b, timer=timer)


def compute_all_tolerances_ill(instance, elems, problem_func, timer=None):
    """Compute l'_P(E) (via ILL) for every non-empty subset E of `elems`."""
    start_t = time()
    S_star, c_star = problem_func(instance, [], [])
    for E in powerset(elems):
        if E:
            ill(instance, list(E), problem_func, S_star, c_star, timer)
    end_t = time() - start_t
    print(f"ILL run time (no rec): {end_t:.3f}")
    if timer is not None:
        print(f"LP solving time: {timer[0]:.3f}, supporting time: {end_t - timer[0]:.3f}")
    return


def ill_rec(instance, E, problem_func, c_star, fc_plus_cache, timer=None):
    """Recursive ILL: reuse fc(D+(F)) across subsets (Remark 5).

    `fc_plus_cache` maps tuple(sorted(F)) -> fc(D+(F)) for every F that has
    already been the top-level `E` of an earlier call. When
    `compute_all_tolerances_ill_rec` sweeps subsets of the ground set in
    order of increasing cardinality, every proper subset F of the current
    E was already processed at an earlier (smaller) size, so its
    fc(D+(F)) is already cached and only fc(D+(E)) itself needs a CSP solve.
    """
    k = len(E)
    if c_star == float('inf'):
        return None, float('inf')

    A, b = [], []
    for F in powerset(E):
        if not F:
            continue
        F = sorted(F)
        F_inds = [i for i in range(k) if E[i] in F]

        if len(F) < k:
            f_c = fc_plus_cache[tuple(F)]
        else:
            _, f_c = problem_func(instance, F, [])
            fc_plus_cache[tuple(F)] = f_c
        if f_c == float('inf'):
            continue

        line = [0] * k
        for i in F_inds:
            line[i] += 1
        A.append(line)
        b.append(f_c - c_star)

    return _solve_tolerance_lp(E, A, b, timer=timer)


def compute_all_tolerances_ill_rec(instance, elems, problem_func, timer=None):
    """Compute l'_P(E) (via recursive ILL) for every non-empty subset of `elems`."""
    start_t = time()
    _, c_star = problem_func(instance, [], [])
    fc_plus_cache = {(): c_star}
    for E in powerset(elems):
        if E:
            ill_rec(instance, list(E), problem_func, c_star, fc_plus_cache, timer)
    end_t = time() - start_t
    print(f"ILL run time (rec): {end_t:.3f}")
    if timer is not None:
        print(f"LP solving time: {timer[0]:.3f}, supporting time: {end_t - timer[0]:.3f}")
    return


# --------------------------------------------------------------------------
# Set lower tolerance, approach 3: TLL
# --------------------------------------------------------------------------

def tll(instance, E, problem_func, S_star=None, c_star=None, timer=None):
    """Tolerance Lower LP: compute l'_P(E) using l'_P(F) for proper subsets.

    Note this recurses into `tll` afresh for every F; 
    see `tll_rec` for the version that shares a tolerance cache 
    across subsets.
    """
    k = len(E)
    if S_star is None or c_star is None:
        S_star, c_star = problem_func(instance, [], [])
    if c_star == float('inf'):
        return None, float('inf')

    if k == 1:
        _, f_c = problem_func(instance, E, [])
        return {tuple(E[0]): f_c - c_star}, f_c - c_star

    A, b = [], []
    for F in powerset(E):
        if not F:
            continue
        F = list(F)
        F_inds = [i for i in range(k) if E[i] in F]

        if len(F) < k:
            _, f_c = tll(instance, F, problem_func, S_star, c_star)
        else:
            _, f_c = problem_func(instance, F, [])
            f_c -= c_star
        if f_c == float('inf'):
            continue

        line = [0] * k
        for i in F_inds:
            line[i] += 1
        A.append(line)
        b.append(f_c)

    return _solve_tolerance_lp(E, A, b, timer=timer)


def compute_all_tolerances_tll(instance, elems, problem_func, timer=None):
    """Compute l'_P(E) (via TLL, no caching) for every non-empty subset of `elems`."""
    start_t = time()
    S_star, c_star = problem_func(instance, [], [])
    for E in powerset(elems):
        if not E:
            continue
        tll(instance, list(E), problem_func, S_star, c_star, timer)
    end_t = time() - start_t
    print(f"TLL run time (no rec): {end_t:.3f}")
    if timer is not None:
        print(f"LP solving time: {timer[0]:.3f}, supporting time: {end_t - timer[0]:.3f}")
    return


def tll_rec(instance, E, problem_func, c_star, tolerances_cache, timer=None):
    """Recursive TLL: reuse l'_P(F) across subsets (Remark 5).

    `tolerances_cache` maps tuple(sorted(F)) -> l'_P(F) for every F that
    has already been the top-level `E` of an earlier call, populated by
    the caller (see `compute_all_tolerances_tll_rec`) after each call
    returns. As with `ill_rec`, sweeping subsets by increasing cardinality
    means every proper subset F of the current E already has a cached
    tolerance, so only fc(D+(E)) needs a CSP solve.
    """
    k = len(E)
    if c_star == float('inf'):
        return None, float('inf')

    if k == 1:
        _, f_c = problem_func(instance, E, [])
        return {tuple(E[0]): f_c - c_star}, f_c - c_star

    A, b = [], []
    for F in powerset(E):
        if not F:
            continue
        F = sorted(F)
        F_inds = [i for i in range(k) if E[i] in F]

        if len(F) < k:
            f_c = tolerances_cache[tuple(F)]
        else:
            _, f_c = problem_func(instance, F, [])
            f_c -= c_star
        if f_c == float('inf'):
            continue

        line = [0] * k
        for i in F_inds:
            line[i] += 1
        A.append(line)
        b.append(f_c)

    return _solve_tolerance_lp(E, A, b, timer=timer)


def compute_all_tolerances_tll_rec(instance, elems, problem_func, timer=None):
    """Compute l'_P(E) (via recursive TLL) for every non-empty subset of `elems`."""
    start_t = time()
    _, c_star = problem_func(instance, [], [])
    tolerances_cache = {(): c_star}
    for E in powerset(elems):
        if not E:
            continue
        _, cost = tll_rec(instance, list(E), problem_func, c_star, tolerances_cache, timer)
        tolerances_cache[tuple(sorted(E))] = cost
    end_t = time() - start_t
    print(f"TLL run time (rec): {end_t:.3f}")
    if timer is not None:
        print(f"LP solving time: {timer[0]:.3f}, supporting time: {end_t - timer[0]:.3f}")
    return
