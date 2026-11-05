#!/usr/bin/env python3
"""
verify_formula.py — Closed-form formula verification for Canonical Krell Numbers K(n),
Striding Numbers S*(n, k), and Rooted Microstate Volume T(n).

Computes:
  S*(n, k) = sum_{lambda |- n, len=k} (c_lambda)^2 * multinomial(n-k, lambda_1 - 1, ...)
  K(n)     = sum_{k=1}^n S*(n, k)
  T(n)     = n * K(n)
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
    """Returns (Striding component S*(lambda), Rooted microstates T(lambda))."""
    k = len(shape)
    counts = Counter(shape)
    c_lambda = factorial(k) // prod(factorial(c) for c in counts.values())
    multinomial = factorial(n - k) // prod(factorial(p - 1) for p in shape)
    s_star = (c_lambda ** 2) * multinomial
    return s_star, n * s_star


def compute_krell_and_striding(n: int) -> tuple[int, int, dict[int, int]]:
    """Returns (K(n), T(n), {k: S*(n, k)})."""
    shapes = integer_partitions(n)
    krell = 0
    rooted = 0
    striding_k: dict[int, int] = {k: 0 for k in range(1, n + 1)}
    for sh in shapes:
        k = len(sh)
        s_star, t_val = formula_shape(n, sh)
        krell += s_star
        rooted += t_val
        striding_k[k] += s_star
    return krell, rooted, striding_k


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute canonical Krell numbers K(n), Striding numbers S*(n, k), and rooted T(n)."
    )
    parser.add_argument("--min-n", type=int, default=1, help="Minimum n (default: 1)")
    parser.add_argument("--max-n", type=int, default=16, help="Maximum n (default: 16)")
    parser.add_argument("--k-table", action="store_true", help="Display the Striding Triangle S*(n, k)")
    parser.add_argument("--format", choices=["table", "markdown", "raw"], default="table",
                        help="Output format (default: table)")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.k_table:
        print(f"Striding Triangle S*(n, k) for n = {args.min_n}..{args.max_n}:")
        header = ["n"] + [f"k={k}" for k in range(1, args.max_n + 1)] + ["K(n)", "T(n)"]
        if args.format == "markdown":
            print("| " + " | ".join(header) + " |")
            print("|" + "|".join(["---:"] * len(header)) + "|")
        for n in range(args.min_n, args.max_n + 1):
            krell, rooted, str_k = compute_krell_and_striding(n)
            row = [str(n)] + [str(str_k.get(k, "-")) for k in range(1, args.max_n + 1)] + [str(krell), str(rooted)]
            if args.format == "markdown":
                print("| " + " | ".join(row) + " |")
            else:
                print(f"n = {n:2d}: " + " ".join(f"{str_k.get(k, 0):8d}" for k in range(1, n + 1)) + f" | K={krell:12d} | T={rooted:14d}")
        return

    if args.format == "markdown":
        print("| n | Canonical Krell K(n) | Rooted Microstates T(n) |")
        print("|---|---:|---:|")
        for n in range(args.min_n, args.max_n + 1):
            krell, rooted, _ = compute_krell_and_striding(n)
            print(f"| {n:2d} | {krell:18,d} | {rooted:20,d} |")
    elif args.format == "raw":
        k_list = [compute_krell_and_striding(n)[0] for n in range(args.min_n, args.max_n + 1)]
        t_list = [compute_krell_and_striding(n)[1] for n in range(args.min_n, args.max_n + 1)]
        print("K(n) =", k_list)
        print("T(n) =", t_list)
    else:
        print(f"Canonical Krell Numbers K(n) and Rooted T(n) for n = {args.min_n}..{args.max_n}:")
        print("-" * 62)
        print(f"{'n':>2s} | {'Canonical Krell K(n)':>22s} | {'Rooted Volume T(n)':>24s}")
        print("-" * 62)
        for n in range(args.min_n, args.max_n + 1):
            krell, rooted, _ = compute_krell_and_striding(n)
            print(f"{n:2d} | {krell:22,d} | {rooted:24,d}")
        print("-" * 62)


if __name__ == "__main__":
    main()
