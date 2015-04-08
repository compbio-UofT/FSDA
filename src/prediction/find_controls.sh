#TMPDIR=/tmp/
tmp_dir=$TMPDIR #"/filer/aryan/tmp/"
#tmp_dir="/filer/aryan/tmp/"

ref_adrs=$1
prep_adrs=$2

sample_test=$3
sample_ref=$4

chrTest=$5
beginTest=$6
endTest=$7

ks_th=$8

samtools view $sample_ref/$chrTest/__plasma.part.bam $chrTest:$beginTest-$endTest | awk '{print $9}' | ../tools/create_hist.sh > $tmp_dir/__ref_dist_$$

python find_controls.py $ref_adrs/bins.pickle $prep_adrs/bins.pickle $tmp_dir/__ref_dist_$$ $ks_th "$chrTest $beginTest $endTest" 
rm $tmp_dir/__ref_dist_$$

