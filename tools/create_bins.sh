tmp_dir=$TMPDIR
if [ "$tmp_dir" == "" ]; then
	tmp_dir="/tmp/"
fi
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#source $DIR/../config_file.sh

centromeres=$DIR/../files/centromeres

sample=$1
output=$2

bin_size=1000000

#for (( i=20; i<23; i++ ))
for (( i=1; i<23; i++ ))
do
	chr=chr$i
	last=`grep -w $chr $centromeres | tail -n 1 | awk  '{print $3}'`
	for (( begin=0; begin<=last; begin+=bin_size ))
	do
		end=$(( begin+bin_size ))
		cent=`grep -w $chr $centromeres | awk -v begin=$begin -v end=$end 'BEGIN {bad=0} {if ($2 < end && $3 > begin) bad=1 } END {print bad}'`
		if (( cent==1 ))
		then
			continue
		fi
		#>&2 echo $chr $begin $end
		echo $chr $begin $end

		for file in $sample/*.bam
		do
			samtools view $file $chr:$begin-$end | awk '{print $9}' >> $tmp_dir/__frag_sizes_$$
		done
		cat $tmp_dir/__frag_sizes_$$ | $DIR/create_hist.sh 
		rm  $tmp_dir/__frag_sizes_$$
	done
done  > $tmp_dir/__bins_$$

python $DIR/create_bins.py $tmp_dir/__bins_$$ $output 
rm $tmp_dir/__bins_$$

