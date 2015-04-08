home="/home/aryan/fsda/feb/"

sample=$1
ref_adrs=$2
tmp_dir=/tmp/

../tools/create_bins.sh $sample $ref_adrs/bins.pickle

exit

samtools view $sample/chr1/__plasma.part.bam chr1:10000000-110000000 | awk '{print $9}' | ../tools/create_hist.sh > $tmp_dir/__frags_$$

pypy mixture.py $tmp_dir/__frags_$$ 25 > $ref_adrs/mog_params

rm $tmp_dir/__frags_$$

