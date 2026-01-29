#!/usr/bin/env python3
"""
Compute length stats for exported per-region FASTA files under `regions/`.

For each `regions/region_*` directory, read all `*.fasta` files:
 - If all lengths are equal, that value is the region length.
 - Else, take the modal length as the region's representative length and note the spread.

Outputs:
 - Global stats across regions (count, min, max, mean, median, stdev, N50, total bp)
 - Top-10 longest and shortest regions
 - Count and list of regions with length spread > 0 (i.e., not perfectly equal)

This script is read-only and does not modify any files.
"""

from __future__ import annotations

import glob
import os
import statistics
from collections import Counter, defaultdict
from typing import Dict, List, Tuple


def read_fasta_len(path: str) -> int:
    """Return the length of the first record in a FASTA file (robust to wrapping)."""
    seq_len = 0
    with open(path, "r") as fh:
        header_seen = False
        for line in fh:
            if not header_seen:
                if line.startswith(">"):
                    header_seen = True
                # ignore any preamble until first '>'
                continue
            if line.startswith(">"):
                break
            seq_len += len(line.strip())
    return seq_len


def n50(lengths: List[int]) -> int:
    if not lengths:
        return 0
    total = sum(lengths)
    half = total / 2
    for L in sorted(lengths, reverse=True):
        half -= L
        if half <= 0:
            return L
    return 0


def collect_region_lengths(regions_root: str = "regions") -> Tuple[List[Tuple[str, int]], Dict[str, Dict]]:
    region_lengths: List[Tuple[str, int]] = []
    diagnostics: Dict[str, Dict] = {}

    for region_dir in sorted(glob.glob(os.path.join(regions_root, "region_*"))):
        if not os.path.isdir(region_dir):
            continue
        fasta_files = sorted(glob.glob(os.path.join(region_dir, "*.fasta")))
        if not fasta_files:
            continue

        lens = []
        for fp in fasta_files:
            try:
                L = read_fasta_len(fp)
                lens.append((os.path.basename(fp), L))
            except Exception:
                # skip unreadable files but track as 0
                lens.append((os.path.basename(fp), 0))

        only_lengths = [L for _, L in lens]
        counts = Counter(only_lengths)
        mode_len, mode_n = counts.most_common(1)[0]
        min_len, max_len = min(only_lengths), max(only_lengths)
        spread = max_len - min_len

        region_id = os.path.basename(region_dir)
        region_lengths.append((region_id, mode_len))
        diagnostics[region_id] = {
            "min": min_len,
            "max": max_len,
            "mode": mode_len,
            "spread": spread,
            "counts": dict(sorted(counts.items())),
        }

    return region_lengths, diagnostics


def main() -> None:
    region_lengths, diagnostics = collect_region_lengths()
    if not region_lengths:
        print("No regions found.")
        return

    # Global stats over region modal lengths
    ids, lens = zip(*region_lengths)
    total = sum(lens)
    stats = {
        "count": len(lens),
        "min": min(lens),
        "max": max(lens),
        "mean": statistics.fmean(lens),
        "median": statistics.median(lens),
        "stdev": statistics.stdev(lens) if len(lens) > 1 else 0.0,
        "N50": n50(list(lens)),
        "total_bp": total,
    }

    print("=== Region Lengths (by mode per region) ===")
    print(f"Regions       : {stats['count']}")
    print(f"Total bp      : {stats['total_bp']}")
    print(f"Min / Max     : {stats['min']} / {stats['max']} bp")
    print(f"Mean / Median : {stats['mean']:.2f} / {stats['median']:.0f} bp")
    print(f"Std-dev       : {stats['stdev']:.2f} bp")
    print(f"N50           : {stats['N50']} bp")
    print()

    # Top/bottom 10 regions by modal length
    sorted_regions = sorted(region_lengths, key=lambda x: x[1])
    print("Shortest regions:")
    for rid, L in sorted_regions[:10]:
        d = diagnostics[rid]
        print(f" - {rid:>10} : {L} bp  (min={d['min']}, max={d['max']}, spread={d['spread']})")
    print()
    print("Longest regions:")
    for rid, L in sorted_regions[-10:][::-1]:
        d = diagnostics[rid]
        print(f" - {rid:>10} : {L} bp  (min={d['min']}, max={d['max']}, spread={d['spread']})")
    print()

    # Regions with spread > 0
    off = [(rid, diagnostics[rid]) for rid in diagnostics if diagnostics[rid]["spread"] > 0]
    print(f"Regions with non-zero spread: {len(off)}")
    # Show a small sample for readability
    for rid, d in sorted(off, key=lambda x: x[1]["spread"], reverse=True)[:15]:
        counts_str = ", ".join(f"{k}bp×{v}" for k, v in d["counts"].items())
        print(f" - {rid}: min={d['min']}, max={d['max']}, spread={d['spread']}  [{counts_str}]")


if __name__ == "__main__":
    main()

