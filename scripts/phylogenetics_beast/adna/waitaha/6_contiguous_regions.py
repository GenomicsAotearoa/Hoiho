#!/usr/bin/env python3
"""
merge_contiguous_intervals_with_stats.py
---------------------------------------

Reads TSV (chrom  start  end) from STDIN or file
→ writes merged intervals to STDOUT
→ prints region‑length stats to STDERR

Works on Python ≥3.7 (uses only standard library).
"""

import sys, gzip, argparse, statistics
from typing import Iterable, Optional, List, TextIO, Tuple

# ────────────────────────────────────────────────────────────────── helpers
def open_in(fname: Optional[str]) -> Iterable[str]:
    if fname is None or fname == "-":
        yield from sys.stdin
    else:
        opener = gzip.open if fname.endswith((".gz", ".bgz")) else open
        with opener(fname, "rt") as fh:          # type: ignore[arg-type]
            for line in fh:
                yield line

def merged_regions(lines: Iterable[str], lenient: bool = False
                   ) -> Iterable[Tuple[str, int, int]]:
    cur_chr = cur_start = cur_end = None
    for raw in lines:
        if not raw.strip() or raw.startswith("#"):
            continue
        chrom, s_str, e_str = raw.rstrip("\n").split("\t")
        s, e = int(s_str), int(e_str)

        if cur_chr is None:
            cur_chr, cur_start, cur_end = chrom, s, e
            continue

        contiguous = chrom == cur_chr and (
            s == cur_end or (lenient and s <= cur_end)
        )
        if contiguous:
            cur_end = max(cur_end, e)
        else:
            yield cur_chr, cur_start, cur_end
            cur_chr, cur_start, cur_end = chrom, s, e

    if cur_chr is not None:
        yield cur_chr, cur_start, cur_end

def n50(lengths: List[int]) -> int:
    """Return N50 for a set of lengths (0 if empty)."""
    if not lengths:
        return 0
    half = sum(lengths) / 2
    cum = 0
    for L in sorted(lengths, reverse=True):
        cum += L
        if cum >= half:
            return L
    return 0  # shouldn’t reach

# ────────────────────────────────────────────────────────────────── main
def main() -> None:
    ap = argparse.ArgumentParser(
        description="Merge contiguous intervals and print length stats"
    )
    ap.add_argument("infile", nargs="?", help="input TSV (default: stdin)")
    ap.add_argument("--lenient", action="store_true",
                    help="merge any overlap (next_start ≤ cur_end)")
    args = ap.parse_args()

    lengths: List[int] = []
    for chrom, start, end in merged_regions(open_in(args.infile), args.lenient):
        print(f"{chrom}\t{start}\t{end}")
        lengths.append(end - start)

    # ───── stats ───────────────────────────────────────────────────────────
    if not lengths:
        print("No intervals written.", file=sys.stderr, flush=True)
        return

    stats = {
        "count": len(lengths),
        "min": min(lengths),
        "max": max(lengths),
        "mean": statistics.fmean(lengths),
        "median": statistics.median(lengths),
        "stdev": statistics.stdev(lengths) if len(lengths) > 1 else 0.0,
        "N50": n50(lengths),
        "total_bp": sum(lengths),
    }

    print(
        (
            "\n--- Region‑length statistics ---\n"
            f"Count        : {stats['count']}\n"
            f"Total bp     : {stats['total_bp']}\n"
            f"Min          : {stats['min']}\n"
            f"Max          : {stats['max']}\n"
            f"Mean         : {stats['mean']:.2f}\n"
            f"Median       : {stats['median']}\n"
            f"Std‑dev      : {stats['stdev']:.2f}\n"
            f"N50          : {stats['N50']}\n"
        ),
        file=sys.stderr,
        flush=True,
    )

if __name__ == "__main__":
    main()
