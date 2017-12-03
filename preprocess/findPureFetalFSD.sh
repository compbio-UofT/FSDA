tmp_dir=$TMPDIR
if [ "$tmp_dir" == "" ]; then
	tmp_dir="/tmp/"
fi
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#source $DIR/../config_file.sh

sample=$1
prep=$2

#for (( i=1; i<23; i++ ))
for file in $prep/allele_counts/*.allele_counts
do
	filename=`basename $file`
	filename=${filename%.allele_counts}
	cat $file | awk '{sum=$3+$4; l=$3; h=$4; allele=$5; if (h < l) {t=l; l=h; h=t; allele=$6} if (sum>10 && l/sum < 0.1 && l>0) print $1, $2, allele}' > $tmp_dir/__snps_$filename"_"$$

	while read line
	do
		chr=`echo $line | awk '{print $1}'`
		pos=`echo $line | awk '{print $2}'`
		allele=`echo $line | awk '{print $3}'`

		samtools view $sample/$filename".bam"  $chr:$pos-$pos > $tmp_dir/__tmp_reads_$filename"_"$$ 
		python read_extract.py $tmp_dir/__tmp_reads_$filename"_"$$ $pos $allele  >> $tmp_dir/__frags_$filename"_"$$

	done < $tmp_dir/__snps_$filename"_"$$

	rm $tmp_dir/__tmp_reads_$filename"_"$$ 
	rm $tmp_dir/__snps_$filename"_"$$ 
done

wait

for file in $tmp_dir/__frags_*_$$
do
	cat $file >> $tmp_dir/__frags_$$
	rm $file
done

cat $tmp_dir/__frags_$$ | $DIR/../tools/create_hist.sh 
rm $tmp_dir/__frags_$$

