chr=$1
centromeres=/dupa-filer/laci/centromeres
ref_panel=/filer/aryan/1000gp/ALL.$chr"."phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
population_1000gp=/filer/aryan/1000gp/european_samples
tmp_dir=/tmp/
output=/filer/aryan/fsda_files/feb/dbsnp/dbsnp_$chr"."ref
output_filtered=/filer/aryan/fsda_files/feb/dbsnp/dbsnp_filtered_$chr"."ref

#pypy selectSNPs.py  $ref_panel $population_1000gp > /filer/aryan/fsda_files/feb/dbsnp/dbsnp_$chr1"."ref

cat $output | awk '{if ($6<0.99 && $6>0.01) print $0}' > $tmp_dir/__dbsnp_$$

while read line
do
        chr=`echo $line | awk '{print $1}'`
        q=`echo $line | awk '{print $2}'`
        p=`echo $line | awk '{print $3}'`
        cat $tmp_dir/__dbsnp_$$ | awk -v q=$q -v p=$p -v chr=$chr '{if ($1!=chr || ($2>p || $2<q)) print $0}' > $tmp_dir/__dbsnp_$$_tmp
        mv $tmp_dir/__dbsnp_$$_tmp $tmp_dir/__dbsnp_$$
done < <(grep -w "$chr" $centromeres)

mv $tmp_dir/__dbsnp_$$  $output_filtered

