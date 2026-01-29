#!/usr/bin/env python3
"""
For each regions/region_* directory:
1) Load all per-sample FASTA files and compare sequence lengths.
   - If all equal: skip alignment/trim.
   - If (max_len - min_len) <= DIFF_TOLERANCE: run MAFFT + trimAl, then split back and overwrite.
   - If (max_len - min_len) > DIFF_TOLERANCE: print only the non-modal-length samples and abort.

MAFFT and trimAl outputs are temporary and cleaned up after splitting.
"""

import os
import glob
import subprocess
import textwrap
import sys
from collections import Counter, defaultdict

THREADS = "32"          # MAFFT threads
WRAP    = None          # wrap FASTA lines; None = no wrap
DIFF_TOLERANCE = 64     # max allowed (max_len - min_len) before aborting

# --------------------------------------
def write_fasta(path, header, seq, wrap=WRAP):
    """(Re)write a single-sequence FASTA file."""
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        if wrap:
            fh.write(textwrap.fill(seq, wrap) + "\n")
        else:
            fh.write(seq + "\n")

# --------------------------------------
def read_single_fasta(fp):
    """Return (header_id, seq) for the first record in a FASTA file."""
    with open(fp) as fh:
        header = fh.readline().strip()
        if not header.startswith(">"):
            raise ValueError(f"{fp} does not look like FASTA")
        header_id = header[1:].split()[0]
        seq = "".join(line.strip() for line in fh if not line.startswith(">"))
    return header_id, seq

# --------------------------------------
def print_offlength_ids(lengths, mode_len):
    """Group and print IDs whose length != mode_len."""
    grouped = defaultdict(list)
    for hid, L in lengths:
        if L != mode_len:
            grouped[L].append(hid)
    if not grouped:
        print("   (No off-length samples w.r.t. mode.)")
    else:
        for L in sorted(grouped.keys()):
            ids = grouped[L]
            print(f"   Off-length {L} bp ({len(ids)}): " + ", ".join(ids))

# --------------------------------------
def process_region(region_dir: str):
    """Check lengths → optionally MAFFT → trimAl → split for one region directory."""
    print(f"▶  Processing {region_dir}")

    # 1) collect *.fasta files
    fasta_files = [f for f in glob.glob(os.path.join(region_dir, "*.fasta"))
                   if os.path.isfile(f)]
    if not fasta_files:
        print("   (no FASTA files found – skipping)")
        return

    # 2) read sequences & lengths
    id_to_path = {}
    id_to_seq  = {}
    lengths = []

    for fp in sorted(fasta_files):
        hid, seq = read_single_fasta(fp)
        id_to_path[hid] = fp
        id_to_seq[hid]  = seq
        lengths.append((hid, len(seq)))

    # 3) decide action based on length spread
    only_lengths = [L for _, L in lengths]
    min_len = min(only_lengths)
    max_len = max(only_lengths)
    spread = max_len - min_len

    counts = Counter(only_lengths)
    mode_len, mode_n = counts.most_common(1)[0]
    summary = ", ".join(f"{L}bp×{n}" for L, n in sorted(counts.items()))
    print(f"   Length summary: {summary} (mode={mode_len}bp; n={mode_n})")

    if spread == 0:
        print(f"   ✓ lengths equal ({min_len} bp) — skipping alignment")
        return
    elif spread > DIFF_TOLERANCE:
        print(f"   ✗ length spread {spread} bp (> {DIFF_TOLERANCE}) — aborting")
        print_offlength_ids(lengths, mode_len)
        sys.exit(1)
    else:
        print(f"   △ lengths differ by up to {spread} bp (≤ {DIFF_TOLERANCE}) — aligning & trimming")
        print_offlength_ids(lengths, mode_len)

    # 4) build combined FASTA (only when needed)
    combined   = os.path.join(region_dir, "combined.fa")
    aln_file   = combined + ".aln"
    trimmed    = combined + ".trim.fa"

    with open(combined, "w") as out:
        for hid in sorted(id_to_seq.keys()):
            out.write(f">{hid}\n{id_to_seq[hid]}\n")

    # 5) MAFFT
    with open(aln_file, "w") as out_aln:
        subprocess.run(
            ["mafft", "--auto", "--thread", THREADS, combined],
            stdout=out_aln,
            stderr=subprocess.DEVNULL,
            check=True,
        )

    # 6) trimAl
    subprocess.run(
        ["trimal", "-in", aln_file, "-out", trimmed, "-nogaps"],
        stderr=subprocess.DEVNULL,
        check=True,
    )

    # 7) split the trimmed alignment back out
    current_id, seq_lines = None, []
    with open(trimmed) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    write_fasta(id_to_path[current_id], current_id, "".join(seq_lines))
                current_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_id is not None:
            write_fasta(id_to_path[current_id], current_id, "".join(seq_lines))

    # 8) clean up
    for tmp in [combined, aln_file, trimmed]:
        try:
            os.remove(tmp)
        except FileNotFoundError:
            pass

    print(f"   ✓ realigned & trimmed {len(fasta_files)} samples")

# --------------------------------------
if __name__ == "__main__":
    for region in sorted(glob.glob("regions/region_*")):
        if os.path.isdir(region):
            process_region(region)
