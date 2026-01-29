grep "chromosome Z" GCA_027474245.1_bSphHub1.pri_genomic.fna
grep "chromosome W" GCA_027474245.1_bSphHub1.pri_genomic.fna
samtools faidx GCA_027474245.1_bSphHub1.pri_genomic.fna "CM049923.1" > chrZ.fa
samtools faidx GCA_027474245.1_bSphHub1.pri_genomic.fna CM049922.1 > chrW.fa
minimap2 -t 48 -x asm20 -o ../chrZ.paf  ../GCA_010078485.1_BGI_Mant.V1_genomic.fna chrZ.fa
minimap2 -t 48 -x asm20 -o ../chrW.paf  ../GCA_010078485.1_BGI_Mant.V1_genomic.fna chrW.fa
cut -f 6 chrZ.paf | sort | uniq > chrZ_possible_chrs
cut -f 6 chrW.paf | sort | uniq > chrW_possible_chrs
cat chrW_possible_chrs chrZ_possible_chrs | sort | uniq > possible_sex_chrs
# Check python notebook to get the non_sex_contigs file
gatk3 -T SelectVariants -R GCA_010078485.1_BGI_Mant.V1_genomic.fna -V merged.filtered.qual20_fmissing0.2.2alleles.snpsonly.pp6_pp19_removed.filtered_qual99.maf0.05.vcf -L non_sex_contigs.list -o no_sex_chrs.vcf

bcftools view -Ob -o no_sex_chrs.bcf no_sex_chrs.vcf
plink2 --bcf no_sex_chrs.bcf --pca --allow-extra-chr --vcf-half-call missing --out no_sex_chrs
