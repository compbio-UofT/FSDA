tmp_dir=$TMPDIR
if [ "$tmp_dir" == "" ]; then
	tmp_dir="/tmp/"
fi
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/../config_file.sh

admixture=$1
cov=$2
cnv=$3
bin_size=$4
ks_th=0.003
ne_th=1
mog=nomog
sum=sum

#sample=/dupa-filer/laci/I1
#prep=/filer/aryan/fsdaFINAL/prep/


if [ $# -eq 0 ]
then
	echo "No argument was given"
	echo "./test_all_chrs.sh admixture bin_size cov cnv"
#	exit
fi

for bin_size in 1000000
do
	batch_size=100000000/$bin_size  
	echo $bin_size
	echo $batch_size
	for admixture in  0.07
	do
		echo $admixture
		for cov in 20
		do
			echo $cov
			for cnv in nor dup del
			do
				echo $cnv
				sample=/filer/aryan/fsdaFINAL/downsampled_bams/$cov"x/"
				prep=/filer/aryan/fsdaFINAL/prep$cov"x/"
				for (( c=1; c<23; c++ ))
				do
					end=`grep -w chr$c $DIR/../files/centromeres | tail -n 1 | awk '{print $3}'`
					for (( i=0; i<=$end; i+=$bin_size*$batch_size ))
					do
						qsub -o /dev/null -e /dev/null  -cwd -q all.q -V -R y -l h_vmem=20G -l h_rt=15:00:00 -S /bin/bash tester_wrapper.sh $sample $prep chr$c $i $batch_size $bin_size $ks_th $ne_th $mog $sum $admixture $cnv  res_$admixture"_"$cov"_"$cnv"_"$bin_size"_"$ks_th"_"$ne_th"_"$mog"_"$sum"_"$c"_"$i

					done
				done
			done
		done
	done
done

