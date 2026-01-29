bcftools index merged.filtered.qual20_fmissing0.2.2alleles.snpsonly.pp6_pp19_removed.filtered_qual99.maf0.05.bcf
bcftools stats -s - merged.filtered.qual20_fmissing0.2.2alleles.snpsonly.pp6_pp19_removed.filtered_qual99.bcf > bcftools_stats
grep "^AF" bcftools_stats > af
grep "^PSC" bcftools_stats > psc

bcftools roh merged.filtered.qual20_fmissing0.2.2alleles.snpsonly.pp6_pp19_removed.filtered_qual99.bcf > roh

bcftools +check-sparsity --n-markers 100 merged.unfiltered.bcf > sparsity
bcftools +guess-ploidy -vr VULA01006153.1 merged.filtered.qual20_fmissing0.2.2alleles.snpsonly.pp6_pp19_removed.filtered_qual99.maf0.05.bcf > guess_ploidy
