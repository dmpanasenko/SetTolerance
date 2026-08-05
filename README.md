# Set Tolerances

Reference implementation accompanying:

> G. Jäger, D. Panasenko, "Computation of Set Tolerances with Applications
> to the Minimum Spanning Tree Problem", [arXiv:2605.10542](https://arxiv.org/abs/2605.10542).

Given a combinatorial sum problem (CSP) and a subset `E` of its ground set,
this code computes the **regular set upper/lower tolerance** of `E`: the
maximum total cost increase/decrease of the elements of `E` such that every
current optimal solution remains optimal (upper), or such that the optimal
objective value is unchanged (lower). See the paper for formal definitions.

## LP formulations

| Function              | Paper name | Computes                | Notes |
|-----------------------|------------|--------------------------|-------|
| `eul`                 | EUL        | upper tolerance          | Theorem 7, constraints (1)-(2) |
| `ell`                 | ELL        | lower tolerance          | Theorem 10, constraint (3) |
| `ill` / `ill_rec`     | ILL        | lower tolerance          | Theorem 13, constraint (4); `_rec` caches `fc(D+(F))` across subsets |
| `tll` / `tll_rec`     | TLL        | lower tolerance          | Theorem 13, constraints (5)-(6); `_rec` caches tolerances across subsets |

ELL, ILL and TLL are all provably equal (Theorem 13). 

When computing the tolerance for *every*
subset of an `m`-element ground set, ELL needs `3^m` CSP solves in total
while the recursive ILL/TLL variants need only `2^m` (Remark 5) — use
`compute_all_tolerances_ill_rec` / `_tll_rec` rather than `_ell` for that
case.

## Problem instances

Three example CSPs (`problem_func` implementations) are included:

- `mstp` — Minimum Spanning Tree Problem, 
  operating on an `n x n` cost matrix with `-1` marking non-edges.
- `lap` — Linear Assignment Problem, 
  operating on an `n x n` assignment cost matrix.
- `spp` — Shortes Path Problem from vertex `s` to vertex `t`
  in directed acyclic graph, 
  operating on `n x n` cost matrix with `-1` marking non-edges.

For provided implementations it is assumed that all costs are non-negative.

Write your own `problem_func(instance, given, banned) -> (solution, cost)`
to apply the tolerance LPs to a different CSP; `cost` must be `float('inf')`
when `given`/`banned` make the instance infeasible.

## Generating test instances

Three generator scripts create random instances in the format each `read_*`
function expects. All three follow the same generate → validity-check →
retry-until-valid → write pattern, so a generated instance is always usable
without further checking:

| Script | Usage | Guarantees |
|---|---|---|
| `gen_mstp.py` | `python gen_mstp.py [fname] [n] [k] [lb] [ub]` | `n` vertices, `k` edges chosen from the `n*(n-1)/2` possible ones, costs uniform in `[lb, ub]`; retries until the graph is **connected**. |
| `gen_lap.py` | `python gen_lap.py [fname] [n] [lb] [ub]` | `n x n` cost matrix, every cost uniform in `[lb, ub]`; no retry needed — a complete cost matrix is always feasible, so there's no `k`. |
| `gen_spp.py` | `python gen_spp.py [fname] [n] [k] [lb] [ub]` | `n` vertices, `k` arcs chosen from the `n*(n-1)/2` possible *forward* arcs (arcs only ever go from a lower-indexed vertex to a higher one, so the graph is acyclic by construction), costs uniform in `[lb, ub]`. Source is fixed at vertex `1` and target at vertex `n` (1-indexed); retries until the target is **reachable** from the source. |

```bash
python gen_mstp.py example_instances/mstp.txt 8 12 1 10
python gen_lap.py example_instances/lap.txt 5 1 10
python gen_spp.py example_instances/spp.txt 8 10 1 10
```

Note `gen_spp.py` doesn't take `s`/`t` as parameters — they're always vertex
`1` and vertex `n` in the output file. If you need a different source/target,
adjust the call to `write_spp` at the bottom of the script.

## Usage

```bash
pip install -r requirements.txt
python run_experiments.py --mstp ./example_instances/mstp.txt
python run_experiments.py --lap ./example_instances/lap.txt --methods ill_rec tll_rec
```

Input file format (`read_mstp` / `read_lap`): first line is `n`, followed
by `n` lines of `n` whitespace-separated integers (the cost matrix).\
Input file format (`read_spp`): first line is `n`, second line is `s` `t`, 
followed by `n` lines of `n` whitespace-separated integers (the cost matrix).

As a library:

```python
from tol import mstp, read_mstp, eul, compute_all_tolerances_ill_rec

graph = read_mstp("mstp2.txt")
n = len(graph)
edges = [(i, j) for i in range(n) for j in range(i + 1, n) if graph[i][j] != -1]

S_star, c_star = mstp(graph, [], [])
increases, tolerance = eul(graph, edges[:2], mstp, S_star, c_star)

compute_all_tolerances_ill_rec(graph, edges, mstp)
```

Example instances for each of the three implemented CSPs are located in the `/example_instances/` folder.

## Complexity warning

Every method here solves at least one LP per *subset* of the ground set
passed to it. This is fine for a handful of elements but grows
exponentially — see the complexity comparison in Remark 5.
