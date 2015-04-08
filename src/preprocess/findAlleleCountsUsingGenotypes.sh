sample=$1
output=$2

#grep chr22 $dbsnp | awk '{print $1,$2}' > $tmp_dir/__tmp_$$.bed
#cat $dbsnp | awk '{print $1, $2}' > $tmp_dir/__tmp_$$.bed

dbsnp=/filer/aryan/dbSnp.red
ref=/filer/aryan/cheo/data/hg19.fa

############## G1 !!
grep -Ff <(grep chr1 $dbsnp | awk '{print $2}')  $sample/chr1/trio.genotype.vcf | awk '{if(substr($10,1,1)==substr($10,3,1) && substr($10,3,1)==1 && substr($12,1,1)!=substr($12,3,1)) print $1, $2}' > /tmp/__snps_$$
#grep -Ff <(grep chr1 $dbsnp | awk '{print $2}')  $sample/chr1/trio.genotype.vcf | awk '{if(substr($10,1,1)==substr($10,3,1) && substr($12,1,1)!=substr($12,3,1)) print $1, $2}' > /tmp/__snps_$$

samtools mpileup -l /tmp/__snps_$$ -f $ref -uDSIB  $sample/chr1/__plasma.part.bam  | bcftools view  -g - | pypy genotyped_vcfParser.py > $output

