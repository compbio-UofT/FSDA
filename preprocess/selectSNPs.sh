tmp_dir=$TMPDIR
if [ "$tmp_dir" == "" ]; then
	tmp_dir="/tmp/"
fi
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/../config_file.sh

centromeres=$DIR/../files/centromeres
population_1000gp=$DIR/../files/european_samples

mkdir -p $config_ref/snplist

chr=$1

ref_panel=$config_1000gp/ALL.$chr"."phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf

pypy $DIR/selectSNPs.py  $ref_panel $population_1000gp | awk '{if ($6<0.99 && $6>0.01) printf "%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5}' > $tmp_dir/__dbsnp_$$ 

while read line
do
	chr=`echo $line | awk '{print $1}'`
	q=`echo $line | awk '{print $2}'`
	p=`echo $line | awk '{print $3}'`
	cat $tmp_dir/__dbsnp_$$ | awk -v q=$q -v p=$p -v chr=$chr '{if ($1!=chr || ($2>p || $2<q)) print $0}' > $tmp_dir/__dbsnp_$$_tmp
	mv $tmp_dir/__dbsnp_$$_tmp $tmp_dir/__dbsnp_$$
done < <(grep -w "$chr" $centromeres)

mv $tmp_dir/__dbsnp_$$ $config_ref/snplist/snplist.$chr".ref"

