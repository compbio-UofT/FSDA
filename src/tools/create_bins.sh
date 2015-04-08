home="/home/aryan/fsda/feb/"
tmp_dir=/tmp/
centromeres=/dupa-filer/laci/centromeres
sample=$1
output=$2

bin_size=1000000
chrs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
#chrs=chr1

for chr in $chrs
do
	last=`grep -w $chr /dupa-filer/laci/centromeres | tail -n 1 | awk  '{print $3}'`
	for (( begin=0; begin<=last; begin+=bin_size ))
	do
		end=$(( begin+bin_size ))
		cent=`grep -w $chr $centromeres | awk -v begin=$begin -v end=$end 'BEGIN {bad=0} {if ($2 < end && $3 > begin) bad=1 } END {print bad}'`
		if (( cent==1 ))
		then
			continue
		fi
		echo $chr $begin $end
		samtools view  $sample/$chr/__plasma.part.bam $chr:$begin-$end | awk '{print $9}' | $home/src/tools/create_hist.sh 
	done
done  > $tmp_dir/__bins_$$

python $home/src/tools/create_bins.py $tmp_dir/__bins_$$ $output #$ref_adrs/bins.ref

