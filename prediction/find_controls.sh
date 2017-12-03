tmp_dir=$TMPDIR
if [ "$tmp_dir" == "" ]; then
	tmp_dir="/tmp/"
fi
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#source $DIR/../config_file.sh

prep=$1

sample=$2

ref=$3

chrTest=$4
beginTest=$5
endTest=$6

ks_th=$7

for file in $ref/*.bam
do
	samtools view $file $chrTest:$beginTest-$endTest | awk '{print $9}'  >> $tmp_dir/__ref_frags_sizes_$$ 
done

cat $tmp_dir/__ref_frags_sizes_$$ | $DIR/../tools/create_hist.sh > $tmp_dir/__ref_dist_$$
rm $tmp_dir/__ref_frags_sizes_$$

python $DIR/find_controls.py $prep/ref_bins.pickle $prep/sample_bins.pickle $tmp_dir/__ref_dist_$$ $ks_th "$chrTest $beginTest $endTest" 
rm $tmp_dir/__ref_dist_$$

