#!/usr/bin/env python3
"""
Fully Integrated Iterative aDNA Rescue and Rescaling Script with ASCII Plot
----------------------------------------------------------------------------

This script takes a raw aligned BAM and performs a full feedback loop for aDNA
read rescue. This version includes a "No Read Left Behind" policy, robust
pydamage handling via pysam, and an "Iterative Relaxation" remapping strategy.

Workflow:
0.  **Initial Rescale (Priming)**:
    - Runs pydamage, but only uses its output to update the original, complete
      set of mapped reads, guaranteeing no data loss.
1.  **Iterative Relaxation Loop**:
    - The script iterates through a user-provided list of BWA -n values.
    a. **Split**: Separates reads into "good" and "to_rescue" pools.
    b. **Remap**: Attempts to remap "to_rescue" reads using the -n value for the current iteration.
    c. **Account for Failures**: Explicitly captures reads that bwa failed to remap.
    d. **Merge & Sort**: Combines all read pools and crucially sorts the result.
    e. **Rescale**: Repeats the robust pydamage update process.
    f. **Repeat**: Continues for each -n value provided by the user.
2.  **Final Split**:
    - Splits the final rescued BAM into two files: one for "Anc" contigs and
      one for all other contigs (contaminant sink).
"""

import subprocess
import argparse
import os
import sys
import shutil
import tempfile
import re
from pathlib import Path
from typing import List, Dict, Any, Optional, TextIO
import pysam

# --- ANSI Color Codes for the Plot ---
C_GREEN = '\033[92m'
C_YELLOW = '\033[93m'
C_RED = '\033[91m'
C_RESET = '\033[0m'
C_BOLD = '\033[1m'

def strip_ansi_codes(text: str) -> str:
    """Removes ANSI escape codes from a string."""
    return re.sub(r'\x1b\[[0-9;]*m', '', text)

def run(cmd: List[Any], cwd: Optional[Path] = None) -> None:
    """Run a shell command, print it for logging, and exit on failure."""
    cmd_str = [str(c) for c in cmd]
    log_prefix = f"(in {cwd}) " if cwd else ""
    print(f"{log_prefix}> {' '.join(cmd_str)}", flush=True)
    try:
        subprocess.check_call(cmd_str, cwd=cwd)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(1)

def run_pydamage_in_pixi_env(
    pixi_manifest: Path, input_bam: Path, output_dir: Path, threads: int
) -> None:
    """Constructs and runs the pydamage command within its specified Pixi environment."""
    output_dir.mkdir(parents=True, exist_ok=True)
    pydamage_base_cmd = ["pydamage", "analyze", "--rescale", "--process", threads, input_bam]
    full_cmd = ["pixi", "run", "--manifest-path", pixi_manifest] + pydamage_base_cmd
    run(full_cmd, cwd=output_dir)

def robust_pydamage_rescale(
    original_mapped_bam: Path,
    final_rescaled_bam: Path,
    pixi_manifest: Path,
    threads: int,
    temp_dir: Path
) -> None:
    """
    Runs pydamage and uses its output to update the original BAM,
    guaranteeing no reads are lost, even if pydamage filters them.
    """
    pydamage_out_dir = temp_dir / "pydamage_run"
    
    run_pydamage_in_pixi_env(pixi_manifest, original_mapped_bam.resolve(), pydamage_out_dir, threads)
    leaky_rescaled_bam = pydamage_out_dir / "pydamage_results" / "pydamage_rescaled.bam"

    print("Creating dictionary of rescaled reads...")
    rescaled_reads = {}
    with pysam.AlignmentFile(leaky_rescaled_bam, "rb", threads=threads) as samfile:
        for read in samfile:
            rescaled_reads[read.query_name] = read

    print("Writing final, complete rescaled BAM to ensure no read loss...")
    with pysam.AlignmentFile(original_mapped_bam, "rb", threads=threads) as infile, \
         pysam.AlignmentFile(final_rescaled_bam, "wb", template=infile, threads=threads) as outfile:
        
        for original_read in infile:
            if original_read.query_name in rescaled_reads:
                outfile.write(rescaled_reads[original_read.query_name])
            else:
                outfile.write(original_read)
    
    print("Robust rescaling complete.")

def get_bam_stats_flagstat(bam_path: Path, mq_thresh: int) -> Dict[str, int]:
    """Gathers counts using a single, fast `samtools flagstat` call."""
    if not bam_path.exists() or bam_path.stat().st_size == 0:
        return {"total": 0, "good": 0, "low_mq": 0, "unmapped": 0}
    
    threads = os.cpu_count() or 1
    flagstat_out = subprocess.check_output(
        ["samtools", "flagstat", f"-@{threads-1}", str(bam_path)]
    ).decode()
    
    total_match = re.search(r"(\d+) \+ \d+ in total", flagstat_out)
    mapped_match = re.search(r"(\d+) \+ \d+ mapped ", flagstat_out)
    
    if not total_match or not mapped_match:
        print(f"Error: Could not parse flagstat output for {bam_path}", file=sys.stderr)
        return {"total": 0, "good": 0, "low_mq": 0, "unmapped": 0}

    total = int(total_match.group(1))
    mapped = int(mapped_match.group(1))
    unmapped = total - mapped
    
    good = int(subprocess.check_output(
        ["samtools", "view", "-c", f"-@{threads-1}", "-q", str(mq_thresh), str(bam_path)]
    ).decode().strip())
    
    low_mq = mapped - good
    
    return {"total": total, "good": good, "low_mq": low_mq, "unmapped": unmapped}

def count_reads_fastq(fastq_path: Path) -> int:
    """Count reads in a FASTQ (4 lines per read)."""
    if not fastq_path.exists() or fastq_path.stat().st_size == 0: return 0
    out = subprocess.check_output(["wc", "-l", str(fastq_path)]).decode().split()[0]
    return int(out) // 4

def get_read_names(bam_path: Path, threads: int) -> Path:
    """Extracts read names (QNAME) of MAPPED reads from a BAM file."""
    output_txt = bam_path.with_suffix(".names.txt")
    print(f"> samtools view -F 4 -@ {threads-1} {bam_path} | cut -f1 > {output_txt}")
    with open(output_txt, "w") as f_out:
        p1 = subprocess.Popen(["samtools", "view", "-F", "4", f"-@{threads-1}", str(bam_path)], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["cut", "-f1"], stdin=p1.stdout, stdout=f_out)
        p1.stdout.close()
        p2.communicate()
        if p1.wait() != 0 or p2.returncode != 0:
            print("Error: Failed to get mapped read names.", file=sys.stderr); sys.exit(1)
    return output_txt

def split_bam_by_mapq(
    input_bam: Path, good_bam: Path, rescue_bam: Path, mq_thresh: int, threads: int
) -> None:
    """
    Robustly splits an input BAM into two non-overlapping files:
    'good' (MAPQ >= mq_thresh) and 'rescue' (all other reads).
    """
    print("Robustly partitioning BAM into 'good' and 'to_rescue' pools...")
    good_read_names = set()
    cmd = ["samtools", "view", f"-@{threads-1}", "-q", str(mq_thresh), str(input_bam)]
    process = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in process.stdout.strip().split('\n'):
        if line:
            good_read_names.add(line.split('\t')[0])

    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as infile:
        with pysam.AlignmentFile(good_bam, "wb", template=infile, threads=threads) as good_outfile, \
             pysam.AlignmentFile(rescue_bam, "wb", template=infile, threads=threads) as rescue_outfile:
            
            for read in infile:
                if read.query_name in good_read_names:
                    good_outfile.write(read)
                else:
                    rescue_outfile.write(read)
    print("Partitioning complete.")

def print_progress_plot(title: str, stats: Dict[str, int], mq_thresh: int, log_fh: Optional[TextIO], bar_width: int = 50):
    """Prints a color-coded ASCII progress bar and detailed stats to console and a log file."""
    total = stats["total"]
    good = stats["good"]
    low_mq = stats["low_mq"]
    unmapped = stats["unmapped"]

    header = f"\n+- {C_BOLD}{title}{C_RESET} {'-' * (bar_width + 15 - len(title))}"
    footer = f"+{'-' * (bar_width + 16)}+"
    
    lines_to_print = [header]
    
    if total == 0:
        lines_to_print.append("| No reads to plot.")
    else:
        good_chars = int((good / total) * bar_width)
        low_mq_chars = int((low_mq / total) * bar_width)
        unmapped_chars = bar_width - good_chars - low_mq_chars

        bar = (f"{C_GREEN}{'¦' * good_chars}"
               f"{C_YELLOW}{'¦' * low_mq_chars}"
               f"{C_RED}{'¦' * unmapped_chars}{C_RESET}")

        lines_to_print.extend([
            f"| {bar} |",
            f"| {' ' * (bar_width + 2)} |",
            f"| {C_GREEN}¦ Good (MAPQ>={mq_thresh}){C_RESET:<48} : {good:>12_} ({good/total:.2%})",
            f"| {C_YELLOW}¦ Low-MQ Mapped{C_RESET:<48} : {low_mq:>12_} ({low_mq/total:.2%})",
            f"| {C_RED}¦ Unmapped{C_RESET:<48} : {unmapped:>12_} ({unmapped/total:.2%})",
            f"| {'-' * (bar_width + 2)} |",
            f"| {C_BOLD}Total Reads{C_RESET:<48} : {total:>12_}"
        ])
    
    lines_to_print.append(footer)

    for line in lines_to_print:
        print(line, flush=True)
        if log_fh:
            log_fh.write(strip_ansi_codes(line) + "\n")
    if log_fh:
        log_fh.flush()

def print_stats_block(title: str, stats: Dict[str, Any], log_fh: Optional[TextIO]):
    """Prints a formatted block of non-plot statistics to console and a log file."""
    header = f"--- {title} ---"
    footer = "-" * (len(title) + 8)
    
    print(header, flush=True)
    if log_fh: log_fh.write(header + "\n")

    for key, value in stats.items():
        line = f"  {key:<35}: {value:_}" if isinstance(value, int) else f"  {key:<35}: {value}"
        print(line, flush=True)
        if log_fh: log_fh.write(line + "\n")

    print(footer, flush=True)
    if log_fh: log_fh.write(footer + "\n"); log_fh.flush()

def split_final_bam(final_bam: Path, output_dir: Path, base_name: str, threads: int):
    """Splits the final BAM into 'Anc' and 'contam' files based on contig names."""
    print("\n>>> Splitting final BAM into target (Anc) and contaminant files <<<")

    anc_bam    = output_dir / f"{base_name}.rescued.anc.bam"
    contam_bam = output_dir / f"{base_name}.rescued.contam.bam"

    # Gather all contig names from the BAM header
    with pysam.AlignmentFile(final_bam, "rb", threads=threads) as sam:
        all_contigs = list(sam.references)
    anc_contigs = [c for c in all_contigs if c.startswith("Anc")]
    contam_contigs = [c for c in all_contigs if c not in anc_contigs]

    if not anc_contigs:
        print(f"{C_YELLOW}Warning: no 'Anc' contigs found in {final_bam}{C_RESET}")
        return final_bam, None

    # 1) extract Anc contigs
    cmd = ["samtools", "view", "-b", f"-@{threads-1}", str(final_bam)] + anc_contigs + ["-o", str(anc_bam)]
    run(cmd)
    run(["samtools", "index", f"-@{threads-1}", str(anc_bam)])

    # 2) extract everything _except_ Anc contigs
    if contam_contigs:
        cmd = ["samtools", "view", "-b", f"-@{threads-1}", str(final_bam)] + contam_contigs + ["-o", str(contam_bam)]
        run(cmd)
        run(["samtools", "index", f"-@{threads-1}", str(contam_bam)])
    else:
        contam_bam = None

    return anc_bam, contam_bam

def main():
    parser = argparse.ArgumentParser(
        description="Fully integrated iterative aDNA rescue and rescaling script with ASCII plot.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--input-bam", required=True, type=Path, help="Initial raw BAM file from alignment (mapped + unmapped).")
    parser.add_argument("--ref", required=True, type=Path, help="Indexed reference FASTA (for BWA).")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Directory to save all output files.")
    parser.add_argument("--pydamage-pixi-manifest", required=True, type=Path, help="Path to the pixi.toml file for the pydamage environment.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for BWA, samtools, and pydamage.")
    parser.add_argument("--mq-thresh", type=int, default=30, help="MAPQ threshold for 'good' reads.")
    # *** FIX: Removed obsolete arguments ***
    parser.add_argument("--n-values", type=float, nargs='+', default=[0.01, 0.02, 0.03], help="List of BWA aln -n mismatch values to use in successive iterations. Default: 0.01 0.02 0.03")
    parser.add_argument("--log-file", type=Path, help="Optional file to write all statistics and progress plots to.")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    base_name = args.input_bam.stem
    
    log_file_handle = open(args.log_file, 'w') if args.log_file else None

    try:
        initial_stats = get_bam_stats_flagstat(args.input_bam, args.mq_thresh)
        print_progress_plot("Initial Input BAM State", initial_stats, args.mq_thresh, log_file_handle)

        initial_rescaled_bam = args.output_dir / f"{base_name}.iter0.bam"
        current_bam = initial_rescaled_bam

        if not initial_rescaled_bam.exists():
            print("\n>>> Performing Initial Rescaling (Priming Step) <<<")
            with tempfile.TemporaryDirectory(prefix="rescue_prime_", dir=args.output_dir) as temp_dir_str:
                temp_dir = Path(temp_dir_str)
                mapped_only = temp_dir / "mapped_only.bam"
                unmapped_only = temp_dir / "unmapped_only.bam"
                robustly_rescaled_mapped = temp_dir / "robustly_rescaled.bam"

                run(["samtools", "view", "-b", f"-@{args.threads-1}", "-F", "4", args.input_bam, "-o", mapped_only])
                run(["samtools", "index", f"-@{args.threads-1}", mapped_only])
                run(["samtools", "view", "-b", f"-@{args.threads-1}", "-f", "4", args.input_bam, "-o", unmapped_only])
                
                robust_pydamage_rescale(mapped_only, robustly_rescaled_mapped, args.pydamage_pixi_manifest, args.threads, temp_dir)
                
                run(["samtools", "merge", "-f", f"-@{args.threads-1}", initial_rescaled_bam, robustly_rescaled_mapped, unmapped_only])
                run(["samtools", "index", f"-@{args.threads-1}", initial_rescaled_bam])
        
        iter0_stats = get_bam_stats_flagstat(initial_rescaled_bam, args.mq_thresh)
        print_progress_plot("Post-Priming State (iter0.bam)", iter0_stats, args.mq_thresh, log_file_handle)
        prev_iter_good_reads = iter0_stats["good"]

        for round_n, n_val in enumerate(args.n_values, 1):
            print(f"\n>>> Starting Iteration {round_n} with -n {n_val} <<<")

            final_bam_iter = args.output_dir / f"{base_name}.iter{round_n}.bam"

            if final_bam_iter.exists() and final_bam_iter.with_suffix(".bam.bai").exists():
                print(f"--- Found existing output for iteration {round_n}. Resuming from next iteration. ---")
                current_bam = final_bam_iter
                prev_iter_good_reads = get_bam_stats_flagstat(current_bam, args.mq_thresh)["good"]
                continue

            with tempfile.TemporaryDirectory(prefix=f"rescue_iter{round_n}_", dir=args.output_dir) as temp_dir_str:
                temp_dir = Path(temp_dir_str)
                
                good_bam_iter = temp_dir / "good.bam"
                to_rescue_bam = temp_dir / "to_rescue.bam"
                
                split_bam_by_mapq(current_bam, good_bam_iter, to_rescue_bam, args.mq_thresh, args.threads)

                rescue_fastq = temp_dir / "rescue.fq"
                run(["samtools", "fastq", f"-@{args.threads}", to_rescue_bam, "-0", rescue_fastq])

                num_candidates = count_reads_fastq(rescue_fastq)
                print_stats_block(f"Iteration {round_n}: Rescue Candidates", {"Reads to rescue": num_candidates}, log_file_handle)
                if num_candidates == 0:
                    print("No reads left to rescue. Finalizing."); break

                sai_file = temp_dir / "rescue.sai"
                bwa_samse_output = temp_dir / "bwa_samse_output.bam"
                
                run(["bwa", "aln", "-l", "1024", "-n", str(n_val), "-o", "2", "-t", str(args.threads), args.ref, rescue_fastq, "-f", sai_file])
                
                samse_cmd = ["bwa", "samse", str(args.ref), str(sai_file), str(rescue_fastq)]
                view_cmd = ["samtools", "view", "-bS", f"-@{args.threads}", "-", "-o", str(bwa_samse_output)]
                print(f"> {' '.join(samse_cmd)} | {' '.join(view_cmd)}", flush=True)
                p1 = subprocess.Popen(samse_cmd, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(view_cmd, stdin=p1.stdout)
                p1.stdout.close()
                p2.communicate()
                if p1.wait() != 0 or p2.returncode != 0:
                    print("Error: BWA or samtools view in the alignment pipe failed.", file=sys.stderr); sys.exit(1)
                
                newly_aligned_bam = temp_dir / "newly_aligned.bam"
                run(["samtools", "view", "-b", f"-@{args.threads-1}", "-F", "4", bwa_samse_output, "-o", newly_aligned_bam])
                
                remapped_names = get_read_names(newly_aligned_bam, args.threads)
                failed_to_remap_bam = temp_dir / "failed_to_remap.bam"
                run(["samtools", "view", "-b", f"-@{args.threads-1}", "-N", remapped_names, str(to_rescue_bam), "-o", str(failed_to_remap_bam)])
                
                newly_rescued_stats = get_bam_stats_flagstat(newly_aligned_bam, args.mq_thresh)
                
                print_stats_block(f"Iteration {round_n}: Remapping Results", {
                    "Candidates re-aligned (any MAPQ)": newly_rescued_stats["total"],
                    f"Rescued with MAPQ>={args.mq_thresh}": newly_rescued_stats["good"],
                }, log_file_handle)

                # *** FIX: Early stopping logic removed ***

                merged_unsorted_bam = temp_dir / "merged.unsorted.bam"
                merged_unrescaled_bam = temp_dir / "merged.unrescaled.bam"
                run(["samtools", "merge", "-f", f"-@{args.threads-1}", merged_unsorted_bam, good_bam_iter, newly_aligned_bam, failed_to_remap_bam])
                run(["samtools", "sort", f"-@{args.threads-1}", "-T", temp_dir / "sort_tmp", merged_unsorted_bam, "-o", merged_unrescaled_bam])

                mapped_for_rescale = temp_dir / "mapped_for_rescale.bam"
                unmapped_final = temp_dir / "unmapped_final.bam"
                robustly_rescaled_mapped_iter = temp_dir / "robustly_rescaled_iter.bam"
                
                run(["samtools", "view", "-b", f"-@{args.threads-1}", "-F", "4", merged_unrescaled_bam, "-o", mapped_for_rescale])
                run(["samtools", "index", f"-@{args.threads-1}", mapped_for_rescale])
                run(["samtools", "view", "-b", f"-@{args.threads-1}", "-f", "4", merged_unrescaled_bam, "-o", unmapped_final])
                
                robust_pydamage_rescale(mapped_for_rescale, robustly_rescaled_mapped_iter, args.pydamage_pixi_manifest, args.threads, temp_dir)
                
                final_bam_unsorted = temp_dir / "final.unsorted.bam"
                run(["samtools", "merge", "-f", f"-@{args.threads-1}", final_bam_unsorted, robustly_rescaled_mapped_iter, unmapped_final])
                run(["samtools", "sort", f"-@{args.threads-1}", "-T", temp_dir / "final_sort_tmp", final_bam_unsorted, "-o", final_bam_iter])
                run(["samtools", "index", f"-@{args.threads-1}", final_bam_iter])

                current_bam = final_bam_iter
                current_stats = get_bam_stats_flagstat(current_bam, args.mq_thresh)
                current_good_reads = current_stats["good"]
                increase = current_good_reads - prev_iter_good_reads
                percent_increase = (increase / prev_iter_good_reads * 100) if prev_iter_good_reads > 0 else float('inf')

                print_progress_plot(f"State After Iteration {round_n}", current_stats, args.mq_thresh, log_file_handle)
                print_stats_block(f"Iteration {round_n} Summary", {
                    "Absolute increase in good reads": increase,
                    "Percentage increase": f"{percent_increase:.2f}%"
                }, log_file_handle)
                prev_iter_good_reads = current_good_reads
        
        print(f"\n>>> Rescue process finished. <<<")
        
        anc_bam, contam_bam = split_final_bam(current_bam, args.output_dir, base_name, args.threads)
        
        print("\n--- Final Output Files ---")
        if anc_bam:
            print(f"Target reads: {C_GREEN}{anc_bam}{C_RESET}")
        if contam_bam:
            print(f"Contaminant reads: {C_YELLOW}{contam_bam}{C_RESET}")

        print("\n--- Recommended Next Steps ---")
        if anc_bam:
            print(f"1. (CRUCIAL) Remove duplicates from target BAM: samtools markdup -r -S {anc_bam} ...")
            print(f"2. (Optional) Hard-trim ends: bam trimBam --in <dedup.bam> --out ... --trim5 2 --trim3 2")
            print("3. Proceed to variant calling with stringent filters.")
        else:
            print("Could not generate final split BAMs. Please check logs for errors.")

    finally:
        if log_file_handle:
            log_file_handle.close()

if __name__ == "__main__":
    main()
