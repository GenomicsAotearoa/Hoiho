../REPrise-rs/target/release/REPrise-rs -input a9_genome_masked.fa -output a9 -dist 1 -pa 16
RepeatMasker -e rmblast -q -lib clstr_out.reprof -xsmall -gff ../a9_genome_masked.fa -pa 16
nice RepeatMasker -e rmblast -q -lib clstr_out.reprof -xsmall -gff ../c90.fasta -pa 16

