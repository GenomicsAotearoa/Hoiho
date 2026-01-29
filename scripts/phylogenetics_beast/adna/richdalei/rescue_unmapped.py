#!/usr/bin/env python3
import subprocess, argparse, os

def run(cmd):
    print(f"> {' '.join(cmd)}")
    subprocess.check_call(cmd)

def count_mapped(bam, mq):
    out = subprocess.check_output(
        ["samtools", "view", "-c", "-q", str(mq), bam]
    ).decode().strip()
    return int(out)

def count_reads_fastq(fq):
    # FASTQ: 4 lines per read
    lines = sum(1 for _ in open(fq))
    return lines // 4

def split_by_length(input_fastq, cuts):
    """
    Splits input_fastq into bins by sequence length.
    cuts: sorted list of cutoffs [c1, c2, ...]
    Returns list of output FASTQ filenames.
    """
    bins = [open(f"bin_{i}.fq", "w") for i in range(len(cuts)+1)]
    with open(input_fastq) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            length = len(seq.strip())
            idx = 0
            for c in cuts:
                if length > c:
                    idx += 1
                else:
                    break
            bins[idx].write(header + seq + plus + qual)
    for b in bins:
        b.close()
    return [f"bin_{i}.fq" for i in range(len(cuts)+1)]

def main():
    p = argparse.ArgumentParser(
        description="Iterative rescue remapping with length bins"
    )
    p.add_argument("--input-bam", required=True, help="Initial rescaled BAM")
    p.add_argument("--reads-fastq", required=True, help="Original FASTQ")
    p.add_argument("--ref", required=True, help="Indexed reference FASTA")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--mq-thresh", type=int, default=30,
                   help="Keep reads with MAPQ >= this")
    p.add_argument("--rescue-thresh", type=int, default=100,
                   help="Min newly rescued reads to continue")
    p.add_argument("--rescue-frac-thresh", type=float, default=0.05,
                   help="Min fraction of rescue reads mapped to continue")
    p.add_argument("--bin-cuts", type=str, default="40,100",
                   help="Comma-separated length cutoffs")
    p.add_argument("--bin-mismatch", type=str, default="0.02,0.015,0.01",
                   help="Comma-separated mismatch rates per bin")
    args = p.parse_args()

    cuts       = sorted(int(x) for x in args.bin_cuts.split(","))
    mismatches = [float(x) for x in args.bin_mismatch.split(",")]
    assert len(mismatches) == len(cuts) + 1, "bin-mismatch length must equal number of bins"

    current_bam = args.input_bam
    prev_good   = count_mapped(current_bam, args.mq_thresh)
    round_n     = 0

    while True:
        round_n += 1
        print(f"\n### Iteration {round_n} ###")

        # 1. split good vs rescue
        good_bam     = f"good_round{round_n}.bam"
        rescue_fastq = f"rescue_round{round_n}.fq"
        run(["samtools", "view", "-b", "-q", str(args.mq_thresh),
             current_bam, "-o", good_bam])
        run(["samtools", "fastq", "-f", "4", "-q", str(args.mq_thresh),
             current_bam, "-0", rescue_fastq])

        rescue_count = count_reads_fastq(rescue_fastq)
        if rescue_count == 0:
            print("No reads left to rescue. Exiting.")
            break

        # 2. split rescue reads by length
        bin_fastqs = split_by_length(rescue_fastq, cuts)

        # 3. remap each bin with its mismatch rate
        bin_bams = []
        for i, fq in enumerate(bin_fastqs):
            sai = f"round{round_n}_bin{i}.sai"
            sam = f"round{round_n}_bin{i}.sam"
            bam = f"round{round_n}_bin{i}.bam"
            mism = mismatches[i]
            run(["bwa", "aln", "-l", "1024", "-n", str(mism), "-o", "2",
                 "-t", str(args.threads), args.ref, fq, "-f", sai])
            run(["bwa", "samse", args.ref, sai, fq, "-f", sam])
            run(["samtools", "view", "-bS", sam])
            run(["samtools", "sort", "-T", f"tmp{round_n}_bin{i}", "-o", bam, "-"])
            bin_bams.append(bam)

        # 4. merge bin bams
        merged_bins = f"merged_round{round_n}_bins.bam"
        run(["samtools", "merge", merged_bins] + bin_bams)

        # 5. merge good + rescued
        merged    = f"merged_round{round_n}.bam"
        final_bam = f"final_round{round_n}.bam"
        run(["samtools", "merge", merged, good_bam, merged_bins])
        run(["samtools", "sort", "-o", final_bam, merged])
        run(["samtools", "index", final_bam])

        # 6. check progress
        current_bam  = final_bam
        good         = count_mapped(current_bam, args.mq_thresh)
        rescued      = good - prev_good
        frac_rescued = rescued / rescue_count
        print(f"Rescued {rescued} reads ({frac_rescued:.2%} of candidates)")

        if rescued < args.rescue_thresh or frac_rescued < args.rescue_frac_thresh:
            print("Plateau reached—stopping.")
            break
        prev_good = good

    print(f"\nFinal BAM: {current_bam}")
    print("Next: end-trim 2–5 bp (e.g., bam trimBam --in {0} --out trimmed.bam --trim5 2 --trim3 2)".format(current_bam))
