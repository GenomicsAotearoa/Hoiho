# Example Bash script: run_raxml.sh
# Make sure raxml-ng is installed and in your PATH.

# Model specification can be changed (e.g., GTR+G, GTR+G+I, etc.).
# Adjust --threads and --bs-trees (bootstrap replicates) to suit your system and analysis preferences.
# --seed sets the random seed for reproducibility.
# --all performs the tree search + bootstrapping + consensus in one step.

for fasta in random_snps_rep*.fasta
do
    prefix=$(basename "$fasta" .fasta)  # e.g. "random_snps_rep1"
    pixi run raxml-ng \
        --all \
        --msa "$fasta" \
        --model GTR+G \
        --bs-trees 100 \
        --threads 4 \
        --seed 12345 \
        --prefix "$prefix"
done

