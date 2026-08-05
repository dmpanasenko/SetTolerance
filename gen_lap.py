"""Generate a random LAP test instance and write it to file, in the
format read by `read_lap` in tol.py (first line n, then n lines of n
integers).

Usage:
    python gen_lap.py [fname] [n] [lb] [ub]

    fname - name of the output file
    n     - size of the (square) cost matrix
    lb    - minimum cost, inclusive
    ub    - maximum cost, inclusive

Costs are drawn independently and uniformly at random from [lb, ub]
for every one of the n x n row/column pairs. Unlike an MSTP or SPP
instance, a LAP instance has no notion of a "missing" arc - every row
can be matched to every column - so there's no connectivity check or
retry loop needed; any n x n cost matrix is a valid instance.
"""

import random
import sys


def gen_lap(n, lb, ub):
    """Return an n x n cost matrix with every cost drawn independently
    and uniformly at random from [lb, ub]."""
    return [[random.choice(range(lb, ub + 1)) for _ in range(n)] for _ in range(n)]


def write_lap(fname, matrix):
    """Writes generated matrix to file, in the format read by `read_lap`
    from `tol.py`."""
    n = len(matrix)

    f = open(fname, 'w')
    f.write(str(n))
    f.write('\n')
    for i in range(n):
        for j in range(n):
            f.write(str(matrix[i][j]))
            f.write(' ')
        f.write('\n')

    return


def main():
    if len(sys.argv) < 5:
        print("Usage: python gen_lap.py [fname] [n] [lb] [ub]")
        return
    else:
        fname = sys.argv[1]
        n = int(sys.argv[2])
        lb = int(sys.argv[3])
        ub = int(sys.argv[4])

    matrix = gen_lap(n, lb, ub)
    write_lap(fname, matrix)

    return


if __name__ == "__main__":
    main()
