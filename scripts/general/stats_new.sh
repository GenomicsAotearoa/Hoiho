INPUT_FILE="merged.a9.filtered.qual99_fmissing0.2.maf0.05.biallelic.bcf"
UNFILTERED_FILE="merged.a9.unfiltered.bcf"

bcftools index "$INPUT_FILE"
bcftools stats -s - "$INPUT_FILE" > bcftools_stats
grep "^AF" bcftools_stats > af
grep "^PSC" bcftools_stats > psc

bcftools roh "$INPUT_FILE" > roh

bcftools +check-sparsity --n-markers 100 "$UNFILTERED_FILE" > sparsity
# bcftools +guess-ploidy -vr "$CHROM" "$INPUT_FILE" > guess_ploidy
