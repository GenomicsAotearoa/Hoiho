#!/usr/bin/env python3
"""
Liftover ancient VCF and re-genotype against a modern target VCF.

Outputs
-------
1. <out-prefix>.bcf                 – Lifted ancient sample, genotyped only at sites present in the modern target VCF.
2. <out-prefix>.tsv                 – One row per input SNP with liftover status, reason, and the new modern GT if successful.
3. <out-prefix>.unlifted_examples.txt – ≤10 example unmapped sites per contig.
4. <out-prefix>.contig_stats.tsv    – Per-contig summary table.
Diagnostics table also printed to stderr.
"""

import argparse, gzip, csv, sys, os, pysam, subprocess
from collections import defaultdict, Counter
from intervaltree import IntervalTree

_complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")

class PSLrow:
    __slots__ = ("strand","qName","tName","blockCount","blockSizes","qStarts","tStarts")
    def __init__(self,f):
        self.strand = f[8]
        self.qName  = f[9]
        self.tName  = f[13]
        self.blockCount = int(f[17])
        self.blockSizes = [int(x) for x in f[18].rstrip(',').split(',')]
        self.qStarts    = [int(x) for x in f[19].rstrip(',').split(',')]
        self.tStarts    = [int(x) for x in f[20].rstrip(',').split(',')]
    def map_pos(self,qpos0:int):
        for size,q0,t0 in zip(self.blockSizes,self.qStarts,self.tStarts):
            if q0 <= qpos0 < q0+size:
                off=qpos0-q0
                if self.strand=='+':
                    return self.tName, t0+off, '+'
                else:
                    return self.tName, t0+(size-1-off), '-'
        return None

# ---------- PSL helpers ----------
def load_psl(psl_path):
    d=defaultdict(list)
    op=gzip.open if psl_path.endswith((".gz",".bgz")) else open
    with op(psl_path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith(("psLayout","match")): continue
            flds=ln.rstrip().split('\t')
            if len(flds)<21: continue
            d[flds[9]].append(PSLrow(flds))
    return d

def build_interval_trees(psl_dict):
    trees={}
    for chrom,rows in psl_dict.items():
        T=IntervalTree()
        for r in rows:
            for sz,q0,_ in zip(r.blockSizes,r.qStarts,r.tStarts):
                T[q0:q0+sz]=r
        trees[chrom]=T
    return trees

# ---------- FASTA helper ----------
def preload_seq(fa_path,chrom,cache):
    if chrom in cache: return cache[chrom]
    with pysam.FastaFile(fa_path) as fa:
        seq=fa.fetch(chrom).upper()
    cache[chrom]=seq
    return seq

# ---------- VCF helper ----------
def load_target_vcf(vcf_path):
    target_sites = {}
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            target_sites[(rec.chrom, rec.pos)] = (rec.ref, rec.alts)
    return target_sites

# ---------- header ----------
def make_header(mod_fa,sample_names):
    hdr=pysam.VariantHeader()
    with pysam.FastaFile(mod_fa) as fa:
        for c in fa.references:
            hdr.contigs.add(c,length=fa.get_reference_length(c))
    hdr.add_meta("source","anc2modern_liftover_diagnostics")
    hdr.add_meta('INFO',items=[('ID','AA'),('Number','1'),('Type','String'), ('Description','Ancestral allele')])
    hdr.add_meta('INFO',items=[('ID','OriginalStrand'),('Number','1'),('Type','String'), ('Description','Original strand orientation')])
    hdr.add_meta('FORMAT',items=[('ID','GT'),('Number','1'),('Type','String'), ('Description','Genotype')])
    for s in sample_names:
        hdr.samples.add(s)
    return hdr

# ---------- per‑contig liftover ----------
# MODIFIED: Handles ambiguous mappings, moves non-SNP filter, accepts contig_order
def liftover_chrom(ch, vcf_in, tree, mod_fa, target_sites, tmp_bcf, tmp_tsv, tmp_ex, samples, hdr, stats, contig_order):
    tsv = open(tmp_tsv, "w", newline='')
    w = csv.writer(tsv, delimiter='\t')
    w.writerow(["lifted", "reason", "modern_GT", "modern_chrom", "modern_pos",
                "liftedoverState", "aDNA_chrom", "aDNA_pos", "aDNA_ref", "aDNA_GT_str", "strand"])
    
    examples = []
    lifted_records = []

    for rec in vcf_in.fetch(ch):
        stats[ch]['total'] += 1
        gt = rec.samples[0].get("GT")
        adna_gt_str = "/".join('.' if i is None else str(i) for i in (gt or []))

        if gt is None:
            # Skip sites with no genotype call
            continue

        # CORRECTED: Moved non-SNP check to the top of the loop for efficiency.
        if len(rec.ref) != 1 or any(len(a) != 1 for a in rec.alts or []):
            stats[ch]['non_snp'] += 1
            w.writerow(["false", "non_snp", "", "NA", "", "", 
                        rec.chrom, rec.pos, rec.ref, adna_gt_str, ""])
            continue
        
        alleles = [rec.ref] + list(rec.alts or [])
        lifted, reason = False, "no_PSL"
        tchr, tpos0, strand, modern_gt_str = "", None, "", ""
        lifted_over_state = ""

        if tree:
            hits = tree[rec.pos - 1]
            if hits:
                # CRITICAL FIX: Check for and reject ambiguous mappings to multiple locations.
                if len(hits) > 1:
                    reason = "ambiguous_psl_map"
                    stats[ch]['ambiguous_psl_map'] += 1
                else:
                    row = next(iter(hits)).data
                    tgt = row.map_pos(rec.pos - 1)
                    if tgt:
                        tchr, tpos0, strand = tgt
                        lifted = True
                        
                        lifted_bases_set = set()
                        original_lifted_bases = []
                        for i in gt:
                            if i is not None:
                                base = alleles[i].upper()
                                if strand == '-': base = base.translate(_complement)
                                lifted_bases_set.add(base)
                                original_lifted_bases.append(base)
                            else:
                                original_lifted_bases.append(None)
                        lifted_over_state = ",".join(sorted(lifted_bases_set))

                        target_site_key = (tchr, tpos0 + 1)
                        if target_site_key in target_sites:
                            target_ref, target_alts = target_sites[target_site_key]
                            target_alleles = [target_ref] + list(target_alts or [])
                            
                            if lifted_bases_set.issubset(set(target_alleles)):
                                reason = "lifted_and_genotyped"
                                stats[ch]['lifted_and_genotyped'] += 1
                                
                                allele_map = {allele: i for i, allele in enumerate(target_alleles)}
                                new_gt_indices = [allele_map.get(b) for b in original_lifted_bases]
                                new_gt = tuple(new_gt_indices)
                                modern_gt_str = "/".join('.' if i is None else str(i) for i in new_gt)

                                new = hdr.new_record(contig=tchr, start=tpos0, alleles=(target_ref, *(target_alts or [])))
                                new.info['AA'], new.info['OriginalStrand'] = rec.ref if strand == '+' else rec.ref.translate(_complement), strand
                                new.samples[0]['GT'] = new_gt
                                lifted_records.append(new)
                            else:
                                reason = "allele_mismatch"
                                stats[ch]['allele_mismatch'] += 1
                        else:
                            reason = "not_in_target_vcf"
                            stats[ch]['not_in_target_vcf'] += 1
                    else:
                        reason = "gap"
            else:
                reason = "no_PSL"

        if not lifted:
            if reason == "no_PSL": stats[ch]['no_psl'] += 1
            elif reason == "gap": stats[ch]['gap'] += 1

        w.writerow([str(lifted).lower(), reason, modern_gt_str, tchr or "NA", 
                    (tpos0 + 1) if tpos0 is not None else "", lifted_over_state, 
                    rec.chrom, rec.pos, rec.ref, adna_gt_str, strand or ""])

        if not lifted and len(examples) < 10:
            examples.append(f"{rec.chrom}:{rec.pos}:{rec.ref}/" + ",".join(rec.alts or []))

    tsv.close()
    with open(tmp_ex, "w") as xf: xf.write("\n".join(examples))

    if lifted_records:
        # CORRECTED: Use the pre-calculated contig_order dictionary passed as an argument.
        lifted_records.sort(key=lambda r: (contig_order.get(r.contig, -1), r.start))
        
        with pysam.VariantFile(tmp_bcf, "w", header=hdr) as vcf_out:
            for rec in lifted_records:
                vcf_out.write(rec)
        
        try:
            subprocess.run(['bcftools', 'index', tmp_bcf], check=True, capture_output=True)
        except FileNotFoundError:
            print("Error: 'bcftools' command not found. Please ensure it is installed and in your system's PATH.", file=sys.stderr)
            sys.exit(1)
        except subprocess.CalledProcessError as e:
            print(f"Error indexing {tmp_bcf} with bcftools. Stderr:\n{e.stderr.decode()}", file=sys.stderr)
            sys.exit(1)

# ---------- MAIN ----------
def main():
    p = argparse.ArgumentParser(description="Liftover ancient VCF and re-genotype against a modern target VCF.")
    p.add_argument("--anc-vcf", required=True, help="Ancient VCF/BCF file to liftover.")
    p.add_argument("--psl", required=True, help="PSL file for coordinate mapping (ancient-to-modern).")
    p.add_argument("--mod-fasta", required=True, help="Modern reference FASTA file.")
    p.add_argument("--target-vcf", required=True, help="Modern VCF file with target sites and alleles.")
    p.add_argument("--out-prefix", default="aDNA_on_mod", help="Prefix for output files.")
    args = p.parse_args()

    psl = load_psl(args.psl)
    trees = build_interval_trees(psl)

    print("Loading target VCF sites...", file=sys.stderr)
    target_sites = load_target_vcf(args.target_vcf)
    print(f"Loaded {len(target_sites)} sites from {args.target_vcf}", file=sys.stderr)

    vcf = pysam.VariantFile(args.anc_vcf)
    samples = list(vcf.header.samples)
    hdr = make_header(args.mod_fasta, samples)

    # CORRECTED: Create contig_order once and add new stat category.
    contig_order = {contig: i for i, contig in enumerate(hdr.contigs)}
    stat_keys = ['total', 'lifted_and_genotyped', 'no_psl', 'gap', 'non_snp', 
                 'not_in_target_vcf', 'allele_mismatch', 'ambiguous_psl_map']
    contig_stats = defaultdict(lambda: Counter({k: 0 for k in stat_keys}))

    vcfs, tsvs, exs = [], [], []
    for ch in vcf.header.contigs:
        tree = trees.get(ch)
        vtmp = f"{args.out_prefix}.{ch}.bcf"
        ttmp = f"{args.out_prefix}.{ch}.tsv"
        etmp = f"{args.out_prefix}.{ch}.unlifted.txt"
        # CORRECTED: Pass contig_order to the function.
        liftover_chrom(ch, vcf, tree, args.mod_fasta, target_sites,
                       vtmp, ttmp, etmp, samples, hdr, contig_stats, contig_order)
        vcfs.append(vtmp); tsvs.append(ttmp); exs.append(etmp)

    # --- Merging and Cleanup ---
    merged_vcf = f"{args.out_prefix}.bcf"
    valid_vcfs = [f for f in vcfs if os.path.exists(f) and os.path.getsize(f) > 0 and os.path.exists(f + ".csi")]

    if valid_vcfs:
        command = ['bcftools', 'concat', '-a', '-o', merged_vcf, '-O', 'b'] + valid_vcfs
        print(f"Running command: {' '.join(command)}", file=sys.stderr)
        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
        except FileNotFoundError:
            print("Error: 'bcftools' command not found. Please ensure it's in your PATH.", file=sys.stderr)
            sys.exit(1)
        except subprocess.CalledProcessError as e:
            print(f"Error running bcftools concat:\nStderr: {e.stderr}", file=sys.stderr)
            sys.exit(1)
    else:
        print("Warning: No variants were successfully lifted. Creating an empty BCF.", file=sys.stderr)
        with pysam.VariantFile(merged_vcf, "wb", header=hdr) as outv: pass

    merged_tsv = f"{args.out_prefix}.tsv"
    with open(merged_tsv, "w") as out:
        first = True
        for f in tsvs:
            if not os.path.exists(f): continue
            with open(f) as fin:
                for i, line in enumerate(fin):
                    if i == 0 and not first: continue
                    first = False
                    out.write(line)

    merged_ex = f"{args.out_prefix}.unlifted_examples.txt"
    with open(merged_ex, "w") as out:
        for f in exs:
            if not os.path.exists(f) or os.path.getsize(f) == 0: continue
            out.write(f"# {os.path.basename(f)}\n")
            out.write(open(f).read() + "\n")
    
    for f_list in [vcfs, tsvs, exs]:
        for f in f_list:
            if os.path.exists(f): os.remove(f)
            if os.path.exists(f + ".csi"): os.remove(f + ".csi")

    # --- Print & save contig stats ---
    # MODIFIED: Added 'ambiguous' to stats output
    stat_file = f"{args.out_prefix}.contig_stats.tsv"
    with open(stat_file, "w") as sf:
        sf.write("contig\ttotal\tgenotyped\tno_psl\tgap\tnot_in_target\tallele_mismatch\tambiguous_map\tnon_snp\tlift_rate(%)\n")
        print("\nContig diagnostics", file=sys.stderr)
        print("%-15s %10s %10s %8s %8s %10s %10s %10s %8s %8s" %
              ("contig", "total", "genotyped", "no_PSL", "gap", "no_target", "mismatch", "ambiguous", "nonSNP", "rate"),
              file=sys.stderr)
        
        totals = Counter()
        for c, st in sorted(contig_stats.items(), key=lambda kv: kv[1]['total'], reverse=True):
            totals.update(st)
            rate = 100 * st['lifted_and_genotyped'] / st['total'] if st['total'] else 0
            print("%-15s %10d %10d %8d %8d %10d %10d %10d %8d %7.2f%%" %
                  (c, st['total'], st['lifted_and_genotyped'], st['no_psl'],
                   st['gap'], st['not_in_target_vcf'], st['allele_mismatch'], 
                   st['ambiguous_psl_map'], st['non_snp'], rate), file=sys.stderr)
            sf.write(f"{c}\t{st['total']}\t{st['lifted_and_genotyped']}\t{st['no_psl']}\t"
                     f"{st['gap']}\t{st['not_in_target_vcf']}\t{st['allele_mismatch']}\t"
                     f"{st['ambiguous_psl_map']}\t{st['non_snp']}\t{rate:.2f}\n")

    print(f"\nGenerated:\n • {merged_vcf}\n • {merged_tsv}\n • {merged_ex}\n • {stat_file}\n", file=sys.stderr)

if __name__ == "__main__":
    main()