dbsnp=/filer/aryan/dbSnp.red
ref=/filer/aryan/cheo/data/hg19.fa
sample=$1

samtools mpileup -l $dbsnp -f $ref -DSIB  $sample/chr1/__plasma.part.bam  -Q20 -q10
