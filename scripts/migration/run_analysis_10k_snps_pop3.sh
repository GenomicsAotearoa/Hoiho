bcftools query -i 'N_MISSING==0' -f '%CHROM\t%POS\n' ../merged.a9.filtered.qual20_fmissing0.2.2alleles.snpsonly.pp6.19.removed.bcf > snp_positions.txt
shuf -n 10000 snp_positions.txt > subset_positions.txt
bcftools view -R subset_positions.txt -Ob -o subset.bcf ../merged.a9.filtered.qual20_fmissing0.2.2alleles.snpsonly.pp6.19.removed.bcf
./ugnix/scripts/bcf2ba3 subset.bcf > hoiho.ba3
./ugnix/scripts/poptrans pop_assignments.txt hoiho.ba3 > hoiho_pops.ba3
echo "./BA3-3.0.5.7/BA3SNP -v -g --trace -o quicktest -u -s 42 ./hoiho_pops.ba3 -b 100000 -i 2000000 -a 0.22 -f 0.004 -m 0.05"
