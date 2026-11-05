#!/usr/bin/env python3
"""
verify_formula.py — Closed-form formula verification for Canonical Transversal Numbers K(n),
Transversal Triangle R(n, k), and Rooted Microstate Volume T(n).

Computes:
  R(n, k) = sum_{lambda |- n, len=k} (c_lambda)^2 * multinomial(n-k, lambda_1 - 1, ...)
  K(n)    = sum_{k=1}^n R(n, k)
  T(n)    = n * K(n)
"""

import argparse
from math import factorial, prod
from collections import Counter


def integer_partitions(n: int, max_part: int | None = None) -> list[list[int]]:
    if max_part is None:
        max_part = n
    if n == 0:
        return [[]]
    res: list[list[int]] = []
    for p in range(min(n, max_part), 0, -1):
        for sub in integer_partitions(n - p, p):
            res.append([p] + sub)
    return res


def formula_shape(n: int, shape: list[int]) -> tuple[int, int]:
    """Returns (Transversal component R(lambda), Rooted microstates T(lambda))."""
    k = len(shape)
    counts = Counter(shape)
    c_lambda = factorial(k) // prod(factorial(c) for c in counts.values())
    multinomial = factorial(n - k) // prod(factorial(p - 1) for p in shape)
    r_val = (c_lambda ** 2) * multinomial
    return r_val, n * r_val


def compute_transversal_numbers(n: int) -> tuple[int, int, dict[int, int]]:
    """Returns (K(n), T(n), {k: R(n, k)})."""
    shapes = integer_partitions(n)
    k_tot = 0
    rooted = 0
    r_k: dict[int, int] = {k: 0 for k in range(1, n + 1)}
    for sh in shapes:
        k = len(sh)
        r_val, t_val = formula_shape(n, sh)
        k_tot += r_val
        rooted += t_val
        r_k[k] += r_val
    return k_tot, rooted, r_k


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute canonical transversal numbers K(n), transversal triangle R(n, k), and rooted T(n)."
    )
    parser.add_argument("--min-n", type=int, default=1, help="Minimum n (default: 1)")
    parser.add_argument("--max-n", type=int, default=16, help="Maximum n (default: 16)")
    parser.add_argument("--k-table", action="store_true", help="Display the Transversal Triangle R(n, k)")
    parser.add_argument("--format", choices=["table", "markdown", "raw"], default="table",
                        help="Output format (default: table)")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.k_table:
        print(f"Canonical Transversal Triangle R(n, k) for n = {args.min_n}..{args.max_n}:")
        header = ["n"] + [f"k={k}" for k in range(1, args.max_n + 1)] + ["K(n)", "T(n)"]
        if args.format == "markdown":
            print("| " + " | ".join(header) + " |")
            print("|" + "|".join(["---:"] * len(header)) + "|")
        for n in range(args.min_n, args.max_n + 1):
            k_tot, rooted, r_k = compute_transversal_numbers(n)
            row = [str(n)] + [str(r_k.get(k, "-")) for k in range(1, args.max_n + 1)] + [str(k_tot), str(rooted)]
            if args.format == "markdown":
                print("| " + " | ".join(row) + " |")
            else:
                print(f"n = {n:2d}: " + " ".join(f"{r_k.get(k, 0):8d}" for k in range(1, n + 1)) + f" | K={k_tot:12d} | T={rooted:14d}")
        return

    if args.format == "markdown":
        print("| n | Canonical K(n) | Rooted Microstates T(n) |")
        print("|---|---:|---:|")
        for n in range(args.min_n, args.max_n + 1):
            k_tot, rooted, _ = compute_transversal_numbers(n)
            print(f"| {n:2d} | {k_tot:18,d} | {rooted:20,d} |")
    elif args.format == "raw":
        k_list = [compute_transversal_numbers(n)[0] for n in range(args.min_n, args.max_n + 1)]
        t_list = [compute_transversal_numbers(n)[1] for n in range(args.min_n, args.max_n + 1)]
        print("K(n) =", k_list)
        print("T(n) =", t_list)
    else:
        print(f"Canonical Transversal Numbers K(n) and Rooted T(n) for n = {args.min_n}..{args.max_n}:")
        print("-" * 62)
        print(f"{'n':>2s} | {'Canonical K(n)':>22s} | {'Rooted Volume T(n)':>24s}")
        print("-" * 62)
        for n in range(args.min_n, args.max_n + 1):
            k_tot, rooted, _ = compute_transversal_numbers(n)
            print(f"{n:2d} | {k_tot:22,d} | {rooted:24,d}")
        print("-" * 62)


if __name__ == "__main__":
    main()
