admixture=$1
cov=$2
cnv=$3
bin_size=$4
ks_th=$5
ne_th=$6
batch_size=10


if [ $# -eq 0 ]
then
	echo "No argument was given"
	echo "./test_all_chrs.sh admixture bin_size cov cnv"
#	exit
fi

for bin_size in 1000000
do
	echo $bin_size
	for admixture in 0.07
	do
		echo $admixture
		for cov in 20
		do
			echo $cov
			for ks_th in 0.002 0.003 0.004
			do
				echo $ks_th
				for ne_th in 1 10 20
				do
					echo $ne_th
					for cnv in del dup nor
					do
						echo $cnv
						for (( c=1; c<23; c++ ))
						do
							end=`grep -w chr$c /dupa-filer/laci/centromeres | tail -n 1 | awk '{print $3}'`
							for (( i=0; i<=$end; i+=$bin_size*$batch_size ))
							do
								current=`qstat | wc -l`
								current=$(( current-375 ))
								while (( current > 800 ))
								do
									sleep 5
									current=`qstat | wc -l`
									current=$(( current-375 ))
								done
								qsub -o /dev/null -e /dev/null -cwd -q all.q -V -R y -l h_vmem=20G -l h_rt=15:00:00 -S /bin/bash tester.sh chr$c $i $batch_size $bin_size $ks_th $ne_th $admixture $cov $cnv  __res_$admixture"_"$cov"_"$cnv"_"$bin_size"_"$ks_th"_"$ne_th"_"$c"_"$i

								#exit
								#qsub -o /dev/null -e __err_$admixture"_"$cov"_"$cnv"_"$bin_size"_"$ks_th"_"$ne_th"_"$c"_"$i -cwd -q all.q -V -R y -l h_vmem=2G -l h_rt=15:00:00 -S /bin/bash tester.sh chr$c $i $batch_size $bin_size $ks_th $ne_th $admixture $cov $cnv  __res_$admixture"_"$cov"_"$cnv"_"$bin_size"_"$ks_th"_"$ne_th"_"$c"_"$i
							done
						done
					done
				done
			done
		done
	done
done
wait

