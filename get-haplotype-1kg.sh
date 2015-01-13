#!/bin/bash
TABIX=/usr/local/genetics/tabix/tabix
VCFTOOLS=/usr/local/genetics/vcftools-0.1.7/cpp/vcftools
PLINK=/usr/local/genetics/plink/plink
CHAOS=/well/lindgren/alexd/chaos/xchaos
HAPLOVIEW=~/bin/Haploview.jar

#$TABIX -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr${1}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ${1}:${2}-${3} > temp.$5.vcf
$TABIX -fh /well/lindgren/alexd/1000genomes/ALL.chr${1}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz ${1}:${2}-${3} > temp.$5.vcf

grep -E 'SNP|#' temp.$5.vcf > temp.$5.snp.vcf

$VCFTOOLS --vcf temp.$5.snp.vcf --keep EUR.sample.txt --plink --out temp.$5

$PLINK --ped temp.$5.ped --map temp.$5.map --recodeHV --out temp.$5 --noweb

java -jar $HAPLOVIEW -nogui -pedfile temp.$5.ped -info temp.$5.info -out temp.$5 -blockoutput

$CHAOS motif vcf2fasta --i temp.$5.snp.vcf --o temp.$5.fasta -lf 1 -rf 1
#cut -f2 temp.$5.map > temp.$5.snp
#grep -f temp.$5.snp chr${1}.fasta | sed 's/    /\n/' > temp.$5.fasta
#grep -A1 -f temp.2.snp chr1.fasta | sed '/^--$/d' > temp.$5.fasta

./find-methyl-alex.pl temp.$5.fasta > temp.$5.methyl.list

./score-blocks-pos2.pl temp.$5.info temp.$5.methyl.list temp.$5.GABRIELblocks ${4} 1

## Calculate Dosage

$PLINK --ped temp.$5.ped --map temp.$5.map  --recodeA --out temp.$5 --noweb

R --vanilla "--args mfile=\"temp.$5.methyl.list\" rfile=\"temp.$5.raw\" sfile=\"hap.markers/$5.$4.hap.markers\"  cg=\"${6}\"" < EUR.score.R
