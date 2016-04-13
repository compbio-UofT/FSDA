tmp_dir=$TMPDIR
if [ "$tmp_dir" == "" ]; then
	tmp_dir="/tmp/"
fi
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/../config_file.sh

prep=$1

sample=$2

chrTest=$3
beginTest=$4
endTest=$5

ks_th=$6

samtools view $config_ref_sample/$chrTest/__plasma.part.bam $chrTest:$beginTest-$endTest | awk '{print $9}' | $DIR/../tools/create_hist.sh > $tmp_dir/__ref_dist_$$

#echo $ks_th "$chrTest $beginTest $endTest" 

python $DIR/find_controls.py $prep/ref_bins.pickle $prep/sample_bins.pickle $tmp_dir/__ref_dist_$$ $ks_th "$chrTest $beginTest $endTest" 
rm $tmp_dir/__ref_dist_$$

