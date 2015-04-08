sample=$1
snps=$2
output=$3

home="/home/aryan/fsda/feb/"
tmp_dir=/tmp/

cat /filer/aryan/NIPT_frag/I1/chr1/fetul_frags | ../tools/create_hist.sh | awk '{for (i=0; i<NF; i++) sum+=$i; for (i=0; i<NF; i++) printf "%f ", $i/sum; printf "\n"}' > $output

exit

cat $snps | awk '{sum=$3+$4; l=$3; h=$4; allele=$5; if (h < l) {t=l; l=h; h=t; allele=$6} if (l/sum < 0.1) print $1, $2, allele}' > $tmp_dir/__snps_$$

while read line
do
	chr=`echo $line | awk '{print $1}'`
	pos=`echo $line | awk '{print $2}'`
	allele=`echo $line | awk '{print $3}'`

	$readWrapper view $chr $pos $pos > $tmp_dir/__tmp_reads_$$ 
	pypy read_extract.py $tmp_dir/__tmp_reads_$$ $pos $allele  >> $tmp_dir/__frags_$$

done < $tmp_dir/__snps_$$

cat $tmp_dir/__frags_$$ | $home/src/tools/create_hist.sh  > $output

