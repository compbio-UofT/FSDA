home="/home/aryan/fsda/feb/"
chr=$1
batch_begin=$2
batch_size=$3
bin_size=$4
ks_th=$5
ne_th=$6
admixture=$7
cov=$8
cnv=$9

code=${10}

cd /home/aryan/fsda/feb/src/tester
centromeres=/dupa-filer/laci/centromeres


tmp_dir="$home/tmp-dir"
ref_dir="$home/gen-ref"
pre_dir="$home/pre-dir"

#prefix=res_$cnv"_"$cov"x"
#prefix=res_strict_$cnv"_"$cov"x"
prefix=res_strict_noMoGs_$cnv"_"$cov"x"
#prefix=res_noMoGs_$cnv"_"$cov"x"

res=$home/res/__$code #$prefix"_"$bin_size"_"$admixture"_"$chr"_"$batch_begin
rm $res
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
	cd ../prediction/
	./predict_region.sh $chr $begin $end $ks_th $ne_th $admixture $cov $cnv >> $res 
done

