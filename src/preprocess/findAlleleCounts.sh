sample=$1
chr=$2
output=$3

#grep chr22 $dbsnp | awk '{print $1,$2}' > $tmp_dir/__tmp_$$.bed
#cat $dbsnp | awk '{print $1, $2}' > $tmp_dir/__tmp_$$.bed

dbsnp=/filer/aryan/fsda_files/feb/dbsnp/dbsnp_filtered_chr1.ref
ref=/filer/aryan/cheo/data/hg19.fa

samtools mpileup -l <(awk '{print $1, $2}' $dbsnp)  -f $ref -t DPR,DP -uIB  $sample/$chr/__plasma.part.bam  -Q15 -q10 | bcftools call - -Am    | pypy vcfParser.py $dbsnp > $output

