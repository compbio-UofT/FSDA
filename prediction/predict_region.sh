tmp_dir=$TMPDIR
if [ "$tmp_dir" == "" ]; then
	tmp_dir="/tmp/"
fi
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/../config_file.sh


################# Parameteres #################

sample=$1
prep=$2

chr=$3
begin=$4
end=$5

ks_th=$6
ne_th=$7

real_fetal_fraction=$8
#admixture=$9
#cnv=${10}

#######################################################################################################################
####################################################### Process #######################################################
#######################################################################################################################

################# Extract test distribution #################


for file in $sample"/"*".bam"
do
	samtools view $file $chr:$begin-$end | awk '{print $9}' >> $tmp_dir/__test_frags_sizes_$$ 
done

cat $tmp_dir/__test_frags_sizes_$$ | $DIR/../tools/create_hist.sh > $tmp_dir/__test_dist_$$ 
rm $tmp_dir/__test_frags_sizes_$$

################# Extract control distributions #################

$DIR/find_controls.sh $prep $sample $chr $begin $end $ks_th > $tmp_dir/__ctrl_dists_$$

neigh=`cat $tmp_dir/__ctrl_dists_$$ | wc -l`
if (( neigh<ne_th ))
then
	echo "Not enough suitable controls for the targeted region." $neigh
	rm $tmp_dir/__ctrl_dists_$$ $tmp_dir/__test_dist_$$
	exit
fi

#####################   Simulation   ######################

sim_r=$real_fetal_fraction

#if [ "$cnv" != "" ]; then
#	pure_fetal_hist=$DIR/pure_fetal_hist_one_line
#	sim_r=$admixture
#	sim_cov=78
#	cnv=$cnv
#	real_r=0.13
#	real_cov=78
#
#	python $DIR/down_admixture.py $tmp_dir/__test_dist_$$ $pure_fetal_hist $real_fetal_fraction $sim_r --c_real $real_cov --c_goal $sim_cov --cnv $cnv >> $tmp_dir/__tmp_$$
#	mv $tmp_dir/__tmp_$$ $tmp_dir/__test_dist_$$
#	python $DIR/down_admixture.py $tmp_dir/__ctrl_dists_$$ $pure_fetal_hist $real_fetal_fraction $sim_r --c_real $real_cov --c_goal $sim_cov >> $tmp_dir/__tmp_$$
#	mv $tmp_dir/__tmp_$$ $tmp_dir/__ctrl_dists_$$
#fi


##################### Prediction ######################

python $DIR/prob_method.py $prep/pure_fetal_hist $sim_r $tmp_dir/__ctrl_dists_$$ $tmp_dir/__test_dist_$$ 
rm $tmp_dir/__ctrl_dists_$$ $tmp_dir/__test_dist_$$

