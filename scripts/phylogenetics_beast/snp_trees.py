import random
import polars as pl
from cyvcf2 import VCF
import sys

# =============================================================================
# Phase 1: Configuration
# =============================================================================
print("Phase 1: Configuring script...")

# --- Input File Paths (Update these to match your system) ---
VCF_PATH = "merged.a9.filtered.qual99_fmissing0.2.maf0.05.biallelic.bcf"
WAITAHA_PATH = 'adna/waitaha/waitaha_filtered.parquet'
RICHDALEI_PATH = 'adna/richdalei/richdalei_filtered.parquet'
CACTUS_PATH = 'adna/cactus_samples.parquet'
HALSTATS_PATH = "seabird_alignment_halstats"

# --- Output Parameters ---
NUM_SNPS = 1000
NUM_REPLICATES = 10
OUTPUT_PREFIX = "random_snp_analysis_n1000"
MIN_CALL_RATE = 0.80 # 80% minimum call rate for Cactus sites

# --- Helper function for IUPAC codes ---
def get_iupac_code(ref, alt):
    """Returns the IUPAC ambiguity code for a heterozygous site."""
    bases = tuple(sorted((ref.upper(), alt.upper())))
    codes = {
        ('A', 'C'): 'M', ('A', 'G'): 'R', ('A', 'T'): 'W',
        ('C', 'G'): 'S', ('C', 'T'): 'Y', ('G', 'T'): 'K',
    }
    # Handle cases with 'N' or other non-standard bases gracefully
    if not ref or not alt: return 'N'
    if ref == 'N': return alt
    if alt == 'N': return ref
    return codes.get(bases, 'N') # Default to N if pair is unexpected

# =============================================================================
# Phase 2: Data Loading and Sample List Preparation
# =============================================================================
print("Phase 2: Loading data and preparing sample lists...")

# --- Load DataFrames lazily to save memory ---
try:
    waitaha_sites_lazy = pl.scan_parquet(WAITAHA_PATH)
    richdalei_sites_lazy = pl.scan_parquet(RICHDALEI_PATH)
    cactus_sites_lazy = pl.scan_parquet(CACTUS_PATH)
except FileNotFoundError as e:
    print(f"ERROR: Could not find an input parquet file: {e}. Please check paths.")
    sys.exit(1)


# --- Define Sample Lists ---
# VCF Samples (Modern a9 penguins)
try:
    vcf = VCF(VCF_PATH)
    a9_samples = vcf.samples
except Exception as e:
    print(f"ERROR: Could not open VCF file '{VCF_PATH}': {e}")
    sys.exit(1)

# aDNA Samples
adna_samples = ["waitaha", "richdalei"]

# Seabird Samples from halstats and Cactus
try:
    # CORRECTED LINE: Added separator='\t'
    seabirds_df = pl.read_csv("seabird_alignment_halstats", skip_lines=4)['GenomeName'].to_list()
    seabirds_df = [s for s in seabirds_df if s is not None]
    seabirds_df = [s for s in seabirds_df if not s.startswith("Anc")]
    # Filter out c90 (it's subantarctic islands and we have them from the regular pop)
    seabirds_df = [s for s in seabirds_df if not s.startswith("c90")]
    seabirds_df = [s for s in seabirds_df if not s.startswith("a9")] # We have a9 in the SNPs as well
    # Remove Megadyptes_antipodes
    seabirds_df = [s for s in seabirds_df if not s.startswith("Megadyptesantipodes")]
    # Filter out samples as in the notebook
    # all_seabirds = [s for s in all_seabirds if not s.startswith(("c90", "a9", "Megadyptesantipodes"))]
    all_seabirds = seabirds_df
except Exception as e:
    print(f"WARNING: Could not read '{HALSTATS_PATH}' correctly: {e}")
    print("Falling back to using column names from cactus parquet for seabird list.")
    all_seabirds = [col.replace('_GT', '') for col in cactus_sites_lazy.collect_schema().names() if col.endswith('_GT')]


# Penguin Samples (subset of seabirds)
penguin_prefixes = ("Aptenodytes", "Spheniscus", "Pygoscelis", "Eudyptula", "Eudyptes")
penguins = [sp for sp in all_seabirds if sp.startswith(penguin_prefixes)]

# --- Define the sets of samples to analyze ---
cactus_penguins = penguins
cactus_seabirds = all_seabirds

default_outgroups = ["Eudyptesmoseleyi_genomic", "Spheniscushumboldti_genomic", "Eudyptesfilholi_genomic"]

TARGET_SAMPLE_SETS = {
    "penguins": a9_samples + adna_samples + cactus_penguins,
    "all_seabirds": a9_samples + adna_samples + cactus_seabirds,
    "vcf_and_adna_only": a9_samples + adna_samples,
    "just_a_few_outgroups": a9_samples + adna_samples + default_outgroups
}
# Ensure uniqueness
for name, sample_list in TARGET_SAMPLE_SETS.items():
    TARGET_SAMPLE_SETS[name] = sorted(list(set(sample_list)))


# =============================================================================
# Phase 3: Generate the Master SNP List
# =============================================================================
print("Phase 3: Generating the master list of high-quality intersecting SNPs...")

# --- Step A: aDNA Intersection ---
print("  - Finding intersection of aDNA sites...")
adna_intersection_lazy = waitaha_sites_lazy.join(
    richdalei_sites_lazy,
    on=["modern_chrom", "modern_pos"],
    how="inner"
).select("modern_chrom", "modern_pos")

# --- Step B: VCF Sites ---
print("  - Extracting all sites from VCF...")
vcf_records = list(vcf) # Read all records into memory
vcf_sites_df = pl.DataFrame({
    "modern_chrom": [r.CHROM for r in vcf_records],
    "modern_pos": [r.POS for r in vcf_records],
    "VCF_ref": [r.REF for r in vcf_records],
    "VCF_alt": [r.ALT[0] if r.ALT else None for r in vcf_records],
}).drop_nulls() # Ensure we only have sites with defined REF and ALT

# --- Step C: Cactus Filter ---
print(f"  - Filtering Cactus sites with >= {MIN_CALL_RATE*100}% call rate...")
cactus_gt_cols = [col for col in cactus_sites_lazy.collect_schema().names() if col.endswith('_GT')]
num_cactus_samples = len(cactus_gt_cols)
min_non_null = int(num_cactus_samples * MIN_CALL_RATE)

filtered_cactus_lazy = cactus_sites_lazy.with_columns(
    pl.sum_horizontal([pl.col(c).is_not_null() for c in cactus_gt_cols]).alias("num_non_null")
).filter(pl.col("num_non_null") >= min_non_null).select("modern_chrom", "modern_pos")

# --- Step D: Final Intersection ---
print("  - Joining all data sources to create the final SNP pool...")
master_snp_list_df = adna_intersection_lazy.join(
    pl.LazyFrame(vcf_sites_df), on=["modern_chrom", "modern_pos"], how="inner"
).join(
    filtered_cactus_lazy, on=["modern_chrom", "modern_pos"], how="inner"
).collect() # Collect the final list into memory for fast sampling

if len(master_snp_list_df) < NUM_SNPS:
    print(f"FATAL ERROR: The final number of intersecting SNPs is {len(master_snp_list_df)}, which is less than the requested {NUM_SNPS}.")
    sys.exit(1)

print(f"  -> Success! Master list created with {len(master_snp_list_df)} eligible SNPs.")

# Collect the source dataframes for faster lookups inside the loop
waitaha_sites_df = waitaha_sites_lazy.collect()
richdalei_sites_df = richdalei_sites_lazy.collect()
cactus_sites_df = cactus_sites_lazy.collect()


# =============================================================================
# Phase 4: Main Replicate Generation Loop
# =============================================================================
print("\nPhase 4: Starting replicate generation...")

for set_name, target_samples in TARGET_SAMPLE_SETS.items():
    if not target_samples:
        print(f"\n--- Skipping empty sample set: '{set_name}' ---")
        continue
    print(f"\n--- Processing sample set: '{set_name}' ({len(target_samples)} samples) ---")

    for rep_index in range(1, NUM_REPLICATES + 1):
        print(f"  Generating replicate {rep_index} of {NUM_REPLICATES}...")

        # --- Step 5a: Sample SNPs ---
        chosen_sites_df = master_snp_list_df.sample(n=NUM_SNPS, with_replacement=False, shuffle=True)

        # --- Step 6a: Initialize sequence holders ---
        nexus_sequences = {sample: [] for sample in target_samples}
        fasta_sequences = {sample: [] for sample in target_samples}

        # --- Step 6b: Iterate through chosen SNPs and build sequences ---
        for site in chosen_sites_df.iter_rows(named=True):
            chrom, pos = site['modern_chrom'], site['modern_pos']
            ref_allele, alt_allele = site['VCF_ref'], site['VCF_alt']

            # Create fast lookup dictionaries for genotypes at this site
            try:
                vcf_gts = {s: g for s, g in zip(a9_samples, next(vcf(f"{chrom}:{pos}")).genotypes)}
            except StopIteration:
                print(f"Warning: Could not find site {chrom}:{pos} in VCF during replicate generation. Skipping site.")
                continue # Skip this site for all samples if it's missing from VCF

            waitaha_gt = waitaha_sites_df.filter((pl.col("modern_chrom") == chrom) & (pl.col("modern_pos") == pos))['waitaha_GT'][0]
            richdalei_gt = richdalei_sites_df.filter((pl.col("modern_chrom") == chrom) & (pl.col("modern_pos") == pos))['richdalei_GT'][0]
            cactus_row = cactus_sites_df.filter((pl.col("modern_chrom") == chrom) & (pl.col("modern_pos") == pos))

            for sample in target_samples:
                nexus_code = '-'
                fasta_code = '-'

                # Get genotype from the correct source
                if sample in adna_samples:
                    gt_str = waitaha_gt if sample == "waitaha" else richdalei_gt
                    if gt_str == "0/0": nexus_code = '0'
                    elif gt_str in ["0/1", "1/0"]: nexus_code = '1'
                    elif gt_str == "1/1": nexus_code = '2'
                elif sample in a9_samples:
                    g1, g2, _ = vcf_gts[sample]
                    if g1 != -1 and g2 != -1: # Check for missing
                        if g1 == 0 and g2 == 0: nexus_code = '0'
                        elif g1 == 1 and g2 == 1: nexus_code = '2'
                        else: nexus_code = '1' # Het
                else: # Must be a cactus sample
                    gt_col_name = sample + "_GT" if not sample.endswith("_GT") else sample
                    if gt_col_name in cactus_row.columns:
                        gt_str = cactus_row[gt_col_name][0]
                        if gt_str == "0/0": nexus_code = '0'
                        elif gt_str == "1/1": nexus_code = '2'
                    else:
                        print(f"Warning: Sample '{sample}' not found in Cactus data for site {chrom}:{pos}. Skipping sample.")
                        continue

                # Convert nexus code to fasta code
                if nexus_code == '0': fasta_code = ref_allele
                elif nexus_code == '1': fasta_code = get_iupac_code(ref_allele, alt_allele)
                elif nexus_code == '2': fasta_code = alt_allele

                nexus_sequences[sample].append(nexus_code)
                fasta_sequences[sample].append(fasta_code)

        # --- Step 7: Write output files for the replicate ---
        out_nexus_path = f"{OUTPUT_PREFIX}_{set_name}_rep{rep_index}.nex"
        out_fasta_path = f"{OUTPUT_PREFIX}_{set_name}_rep{rep_index}.fasta"

        # Write NEXUS
        with open(out_nexus_path, "w") as f:
            f.write("#NEXUS\n")
            f.write(f"[SNP matrix for sample set '{set_name}', replicate {rep_index}]\n")
            f.write("[0=homREF, 1=het, 2=homALT, -=missing]\n\n")
            f.write("Begin data;\n")
            f.write(f"\tDimensions ntax={len(target_samples)} nchar={NUM_SNPS};\n")
            f.write('\tFormat datatype=integerdata symbols="012" gap=-;\n')
            f.write("\tMatrix\n")
            for sample, seq_list in nexus_sequences.items():
                f.write(f"{sample}\t{''.join(seq_list)}\n")
            f.write("\t;\nEnd;\n")

        # Write FASTA
        with open(out_fasta_path, "w") as f:
            for sample, seq_list in fasta_sequences.items():
                f.write(f">{sample}\n{''.join(seq_list)}\n")

        print(f"    -> Wrote NEXUS: {out_nexus_path}")
        print(f"    -> Wrote FASTA: {out_fasta_path}")

print("\nAll replicates for all sample sets are done!")