"""Benchmark script: compute all regular set tolerances for an instance, 
timing each of the LP formulations from tol.py against each other.

Usage:
    python run_experiments.py --spp ./example_instances/spp.txt
    python run_experiments.py --lap ./example_instances/lap.txt --methods eul
    python run_experiments.py --mstp ./example_instances/mstp.txt --methods ill_rec tll_rec

`--methods` accepts either individual method names (eul, ell, ill, ill_rec,
tll, tll_rec) or one of the named groups below, and any mix of the two
(duplicates are dropped, first occurrence wins):

    all-lower    ell, ill, ill_rec, tll, tll_rec    (every lower
                 tolerance method; eul computes the *upper* tolerance 
                 via a different LP, so it's left out of this group)
    ell-rec      ell, ill_rec, tll_rec              (default - plain ELL 
                 vs. both recursive lower-tolerance variants)
    ell-ill-rec  ell, ill, ill_rec, tll_rec         (ELL and ILL 
                 vs. the recursive variants)
    all          every method, including eul

    python run_experiments.py --mstp ./example_instances/mstp.txt --methods ell-rec
    python run_experiments.py --mstp ./example_instances/mstp.txt --methods all-lower eul

Warning: the ground set here has one element per edge (MSTP) / arc (SPP) /
row-column pair (LAP), and every method below solves at least one LP per
subset of that ground set (2^m of them). This is fine for graphs with a
handful of edges but becomes intractable quickly.
"""

import argparse

from tol import (
    compute_all_tolerances_eul,
    compute_all_tolerances_ell,
    compute_all_tolerances_ill,
    compute_all_tolerances_ill_rec,
    compute_all_tolerances_tll,
    compute_all_tolerances_tll_rec,
    mstp,
    lap,
    spp,
    read_mstp,
    read_lap,
    read_spp,
)

METHODS = {
    "eul": compute_all_tolerances_eul,
    "ell": compute_all_tolerances_ell,
    "ill": compute_all_tolerances_ill,
    "ill_rec": compute_all_tolerances_ill_rec,
    "tll": compute_all_tolerances_tll,
    "tll_rec": compute_all_tolerances_tll_rec,
}

PRESETS = {
    "all-lower": ["ell", "ill", "ill_rec", "tll", "tll_rec"],
    "ell-rec": ["ell", "ill_rec", "tll_rec"],
    "ell-ill-rec": ["ell", "ill", "ill_rec", "tll_rec"],
    "all": list(METHODS),
}

DEFAULT_METHODS = ["ell-rec"]


def resolve_methods(selected):
    """Expand any preset names in `selected` into their component method
    names, keep individual method names as-is, and drop duplicates while
    preserving first-occurrence order."""
    resolved = []
    for name in selected:
        for method in PRESETS.get(name, [name]):
            if method not in resolved:
                resolved.append(method)
    return resolved


def run(instance, elems, problem_func, methods):
    for name in resolve_methods(methods):
        timer = [0.0]
        METHODS[name](instance, elems, problem_func, timer)
    return


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--mstp", help="MSTP instance file (see read_mstp)")
    source.add_argument("--lap", help="LAP cost-matrix file (see read_lap)")
    source.add_argument("--spp", help="SPP instance file (see read_spp)")
    parser.add_argument(
        "--methods", nargs="+", choices=list(METHODS) + list(PRESETS),
        default=DEFAULT_METHODS,
        help="which LP formulations to run - individual method names, "
             "preset group names, or a mix of both (default: ell-rec)",
    )
    args = parser.parse_args()

    if args.mstp:
        graph = read_mstp(args.mstp)
        n = len(graph)
        edges = [(i, j) for i in range(n) for j in range(i + 1, n) if graph[i][j] != -1]
        run(graph, edges, mstp, args.methods)
    elif args.lap:
        matrix = read_lap(args.lap)
        n = len(matrix)
        arcs = [(i, j) for i in range(n) for j in range(n)]
        run(matrix, arcs, lap, args.methods)
    else:
        graph, s, t = read_spp(args.spp)
        n = len(graph)
        arcs = [(i, j) for i in range(n) for j in range(n) if graph[i][j] != -1]
        run((graph, s, t), arcs, spp, args.methods)

    return


if __name__ == "__main__":
    main()
