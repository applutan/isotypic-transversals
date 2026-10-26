#!/usr/bin/env python3
"""
verify_formula.py — Test the closed-form formula for total rooted isotypic composition transversals.
"""

import argparse
import os
import sys
from math import factorial, prod


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


def formula_shape(n: int, shape: list[int]) -> int:
    k = len(shape)
    # Part multiplicities
    counts: dict[int, int] = {}
    for p in shape:
        counts[p] = counts.get(p, 0) + 1
    c_lambda = factorial(k) // prod(factorial(c) for c in counts.values())

    multinomial = factorial(n - k) // prod(factorial(p - 1) for p in shape)

    return n * (c_lambda ** 2) * multinomial


def formula_total(n: int) -> int:
    shapes = integer_partitions(n)
    return sum(formula_shape(n, sh) for sh in shapes)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute closed-form total rooted isotypic composition transversals T(n)."
    )
    parser.add_argument("--min-n", type=int, default=1, help="Minimum n (default: 1)")
    parser.add_argument("--max-n", type=int, default=16, help="Maximum n (default: 16)")
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"Closed-form formula T(n) for n = {args.min_n}..{args.max_n}:")
    for n in range(args.min_n, args.max_n + 1):
        tot = formula_total(n)
        print(f"n = {n:2d}: {tot:18d}")


if __name__ == "__main__":
    main()
