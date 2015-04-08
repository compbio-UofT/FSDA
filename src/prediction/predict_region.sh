home="/home/aryan/fsda/feb/"
tmp_dir=$TMPDIR
#tmp_dir="/filer/aryan/tmp/"

################# Parameteres #################

chrTest=$1
beginTest=$2
endTest=$3

ks_th=$4
ne_th=$5

admixture=$6
cov=$7
cnv=$8

sample_test=/dupa-filer/laci/I1 ### TODO
sample_ref=/dupa-filer/laci/G1 ### TODO

prep_adrs=/filer/aryan/fsda_files/feb/prep/
prep_adrs=/filer/aryan/fsda_files/feb/newFiles/I1/prep/
ref_adrs=/filer/aryan/fsda_files/feb/newFiles/G1/ref/

################# Extract Reads #################

samtools view $sample_test/$chrTest/__plasma.part.bam $chrTest:$beginTest-$endTest | awk '{print $9}' | $home/src/tools/create_hist.sh  > $tmp_dir/__test_dist_$$

################# Similarity Network #################

$home/src/prediction/find_controls.sh $ref_adrs $prep_adrs $sample_test $sample_ref $chrTest $beginTest $endTest $ks_th > $tmp_dir/__ctrl_dists_$$
neigh=`cat $tmp_dir/__ctrl_dists_$$ | wc -l`
if (( neigh<ne_th ))
then
	echo " --" $neigh
	exit
fi

#####################   Simulation   ######################

pure_fetal_hist=pure_fetal_hist_one_line
sim_r=$admixture
sim_cov=$cov
cnv=$cnv
#sim_r=0.13
#sim_cov=40
#cnv=nor

python down_admixture.py $tmp_dir/__test_dist_$$ $pure_fetal_hist 0.13 $sim_r --c_real 78 --c_goal $sim_cov --cnv $cnv >> $tmp_dir/__tmp_$$
mv $tmp_dir/__tmp_$$ $tmp_dir/__test_dist_$$
python down_admixture.py $tmp_dir/__ctrl_dists_$$ $pure_fetal_hist 0.13 $sim_r --c_real 78 --c_goal $sim_cov >> $tmp_dir/__tmp_$$
mv $tmp_dir/__tmp_$$ $tmp_dir/__ctrl_dists_$$

##################### Predict ######################

python $home/src/prediction/monteCarlo.py $prep_adrs/fetalDist.prp $sim_r $tmp_dir/__ctrl_dists_$$ $tmp_dir/__test_dist_$$ $ref_adrs/mog_params #| awk '{if ($1!="PREDICTION:") print $1, $2}'  
rm $tmp_dir/__ctrl_dists_$$ $tmp_dir/__test_dist_$$

