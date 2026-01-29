#!/usr/bin/env python3
"""
Filter a BAM to remove contigs with zero mapped reads, updating header and records.

Usage:
    python filter_unused_contigs.py --input filtered.bam [--output pruned.bam]
"""
import argparse
from pathlib import Path
import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove contigs with zero mapped reads from a BAM file."
    )
    parser.add_argument(
        '--input', '-i', required=True,
        type=Path, help='Input BAM file (must be indexed)'
    )
    parser.add_argument(
        '--output', '-o', required=False,
        type=Path, help='Output pruned BAM file (default: input.pruned.bam)'
    )
    parser.add_argument(
        '--min-mapped', type=int, default=1,
        help='Minimum number of mapped reads to keep a contig (default: 1)'
    )
    return parser.parse_args()


def main():
    args = parse_args()
    in_bam = args.input
    if args.output:
        out_bam = args.output
    else:
        stem = in_bam.stem
        out_bam = in_bam.with_name(f"{stem}.pruned.bam")

    # Open input BAM
    infile = pysam.AlignmentFile(str(in_bam), 'rb')

    # Determine contigs with mapped reads
    stats = infile.get_index_statistics()
    keep_contigs = [s.contig for s in stats if s.mapped >= args.min_mapped]
    if not keep_contigs:
        print(f"No contigs with >= {args.min_mapped} mapped reads found. Exiting.")
        return

    # Build new header with only kept contigs
    header = infile.header.to_dict()
    header['SQ'] = [sq for sq in header.get('SQ', []) if sq['SN'] in keep_contigs]

    # Write pruned BAM, fixing reference IDs for reads
    outfile = pysam.AlignmentFile(str(out_bam), 'wb', header=header)
    for read in infile.fetch(until_eof=True):
        ref_name = read.reference_name
        if ref_name not in keep_contigs:
            continue
        # Update reference_id to match new header
        read.reference_id = outfile.get_tid(ref_name)
        # If paired, update next_reference_id as well
        if not read.is_unmapped and read.next_reference_name:
            read.next_reference_id = outfile.get_tid(read.next_reference_name)
        outfile.write(read)

    infile.close()
    outfile.close()

    # Index output BAM
    pysam.index(str(out_bam))

    print(f"Wrote pruned BAM to {out_bam} containing {len(keep_contigs)} contigs.")
    print(f"Index created: {out_bam}.bai")


if __name__ == '__main__':
    main()
