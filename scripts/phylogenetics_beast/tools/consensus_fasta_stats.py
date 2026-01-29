#!/usr/bin/env python3
"""
Compute length stats for a multi-FASTA file (one length per record).

Usage:
  python tools/consensus_fasta_stats.py path/to/file.fasta [path2.fasta ...]

Prints stats for each file: count, min, max, mean, median, stdev, N50, total.
"""

from __future__ import annotations

import sys
import statistics
from typing import List


def read_lengths(fp: str) -> List[int]:
    lens: List[int] = []
    with open(fp, "r") as fh:
        l = 0
        saw = False
        for line in fh:
            if line.startswith(">"):
                if saw:
                    lens.append(l)
                    l = 0
                saw = True
            else:
                l += len(line.strip())
        if saw:
            lens.append(l)
    return lens


def n50(lengths: List[int]) -> int:
    total = sum(lengths)
    half = total / 2
    for L in sorted(lengths, reverse=True):
        half -= L
        if half <= 0:
            return L
    return 0


def stats_for(lengths: List[int]):
    return {
        "count": len(lengths),
        "min": min(lengths),
        "max": max(lengths),
        "mean": statistics.fmean(lengths) if lengths else 0.0,
        "median": statistics.median(lengths) if lengths else 0,
        "stdev": statistics.stdev(lengths) if len(lengths) > 1 else 0.0,
        "N50": n50(lengths) if lengths else 0,
        "total_bp": sum(lengths),
    }


def main(argv: List[str]) -> int:
    if len(argv) < 2:
        print(__doc__.strip())
        return 2
    for path in argv[1:]:
        lens = read_lengths(path)
        if not lens:
            print(f"{path}: no records found")
            continue
        s = stats_for(lens)
        print(f"=== {path} ===")
        print(f"Records      : {s['count']}")
        print(f"Total bp     : {s['total_bp']}")
        print(f"Min / Max    : {s['min']} / {s['max']} bp")
        print(f"Mean/Median  : {s['mean']:.2f} / {s['median']:.0f} bp")
        print(f"Std-dev      : {s['stdev']:.2f} bp")
        print(f"N50          : {s['N50']} bp")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

