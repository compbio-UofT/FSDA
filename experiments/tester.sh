tmp_dir=$TMPDIR
if [ "$tmp_dir" == "" ]; then
	tmp_dir="/tmp/"
fi
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/../config_file.sh


sample=$1
prep=$2
chr=$3
batch_begin=$4
batch_size=$5
bin_size=$6
ks_th=$7
ne_th=$8
mog=$9
sum=${10}
admixture=${11}
cnv=${12}
code=${13}

centromeres=$DIR/../files/centromeres

res=$DIR/../../res_final/$code

rm -f $res

for (( i=0; i<$batch_size; i++ ))
do
	begin=$(( $batch_begin+$i*$bin_size ))
	end=$(( $begin+$bin_size ))
	cent=`grep -w $chr $centromeres | awk -v begin=$begin -v end=$end 'BEGIN {bad=0} {if ($2 < end && $3 > begin) bad=1 } END {print bad}'`
	if (( cent==1 ))
	then
		continue
	fi
	chr_end=`grep -w $chr $centromeres | tail -n 1 | awk '{print $3}'`
	if (( end>chr_end ))
	then
		continue
	fi

	echo $chr $begin $end" "  >> $res
	$DIR/../prediction/predict_region.sh $sample $prep $chr $begin $end $ks_th $ne_th $mog $sum $admixture $cnv >> $res 
done

