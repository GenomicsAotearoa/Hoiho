bcftools query -i 'N_MISSING==0' -f '%CHROM\t%POS\n' ../merged.a9.filtered.qual20_fmissing0.2.2alleles.snpsonly.pp6.19.removed.bcf > snp_positions.txt
shuf -n 20000 snp_positions.txt > subset_positions.txt
bcftools view -R subset_positions.txt -Ob -o subset.bcf ../merged.a9.filtered.qual20_fmissing0.2.2alleles.snpsonly.pp6.19.removed.bcf
./ugnix/scripts/bcf2ba3 subset.bcf > hoiho.ba3
./ugnix/scripts/poptrans ../pop_assignments_comprehensive.txt hoiho.ba3 > hoiho_pops_comprehensive.ba3


