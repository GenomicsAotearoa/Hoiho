#!/usr/bin/env python3
import sys
import argparse
import gzip
import io
import re
from collections import defaultdict

def open_maybe_gzip(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r", encoding="utf-8")

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert EggNOG-mapper annotations to GOATOOLS association format."
    )
    p.add_argument("input", help="EggNOG-mapper .tsv (or .tsv.gz). Use '-' for stdin.")
    p.add_argument("-o", "--output", default="-", help="Output file (default: stdout).")
    p.add_argument(
        "--gene-extract",
        default="split",
        choices=["split", "regex", "full"],
        help="How to derive gene ID from query ID. "
             "'split' -> take text before first '.'; "
             "'regex' -> use --gene-regex; "
             "'full'  -> use full query as gene ID."
    )
    p.add_argument(
        "--gene-regex",
        default=r"^(?P<gene>[^.\s]+)",
        help="Regex with named group 'gene' to extract gene ID when --gene-extract=regex."
    )
    p.add_argument(
        "--keep-empty",
        action="store_true",
        help="Include genes with no GO terms (will produce empty right-hand side)."
    )
    return p.parse_args()

def gene_from_query(q, mode, regex):
    if mode == "full":
        return q
    if mode == "split":
        # split on first '.' (common for transcript suffixes), fall back to whole
        return re.split(r"\.", q, maxsplit=1)[0]
    # regex mode
    m = re.match(regex, q)
    return m.group("gene") if m else q

def main():
    args = parse_args()
    go_by_gene = defaultdict(dict)  # dict-as-ordered-set of GO terms per gene

    with open_maybe_gzip(args.input) as fh:
        header = None
        go_idx = None
        query_idx = None

        for line in fh:
            if not line.strip():
                continue
            # Header line starts with '#query' in emapper v2; other metadata begins with '##' or '#'
            if line.startswith("#"):
                if line.startswith("#query"):
                    header = line.lstrip("#").rstrip("\n").split("\t")
                    # Find positions of key columns
                    try:
                        query_idx = header.index("query")
                    except ValueError:
                        sys.exit("ERROR: couldn't find 'query' column in header")
                    try:
                        go_idx = header.index("GOs")
                    except ValueError:
                        sys.exit("ERROR: couldn't find 'GOs' column in header")
                continue

            # Past here: data rows
            if header is None or go_idx is None or query_idx is None:
                # If the file lacks a #query header for some reason, try best-effort by counting columns
                parts = line.rstrip("\n").split("\t")
                # Fall back: assume emapper default schema with 'GOs' as the 10th column (index 9)
                if query_idx is None:
                    query_idx = 0
                if go_idx is None:
                    go_idx = 9 if len(parts) > 9 else None
                    if go_idx is None:
                        sys.exit("ERROR: Could not determine GOs column index.")
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(query_idx, go_idx):
                continue

            query_id = parts[query_idx]
            gos_raw = parts[go_idx].strip()

            gene_id = gene_from_query(query_id, args.gene_extract, args.gene_regex)

            # Skip empty or '-' GOs unless --keep-empty is set
            if gos_raw == "-" or gos_raw == "" or gos_raw is None:
                if args.keep_empty:
                    go_by_gene.setdefault(gene_id, {})
                continue

            # GO terms are comma-separated in emapper output
            for tok in gos_raw.split(","):
                term = tok.strip()
                if not term:
                    continue
                if not term.startswith("GO:"):
                    continue
                # Use dict to preserve insertion order & uniqueness
                if term not in go_by_gene[gene_id]:
                    go_by_gene[gene_id][term] = True

    # Write output
    out = sys.stdout if args.output == "-" else open(args.output, "w", encoding="utf-8")
    with out:
        for gene in sorted(go_by_gene.keys()):
            terms = list(go_by_gene[gene].keys())
            if not terms and not args.keep_empty:
                continue
            # Semicolon-separated, no spaces
            rhs = ";".join(terms)
            out.write(f"{gene}\t{rhs}\n")

if __name__ == "__main__":
    main()
