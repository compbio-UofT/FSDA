dbsnp=/dupa-filer/laci/dbSnp.vcf
centromeres=/dupa-filer/laci/centromeres
output=/filer/aryan/dbSnp.red
tmp_dir=/tmp/

cat $dbsnp | awk '{if (substr($1,0,1)!="#" && length($4)==1 && length($5)==1) printf "chr%d\t%d\t%c\t%c\n", $1, $2, $4, $5}' > $tmp_dir/__dbsnp_$$

while read line
do
	chr=`echo $line | awk '{print $1}'`
	q=`echo $line | awk '{print $2}'`
	p=`echo $line | awk '{print $3}'`
	cat $tmp_dir/__dbsnp_$$ | awk -v q=$q -v p=$p -v chr=$chr '{if ($1!=chr || ($2>p || $2<q)) print $0}' > $tmp_dir/__dbsnp_$$_tmp
	mv $tmp_dir/__dbsnp_$$_tmp $tmp_dir/__dbsnp_$$
done < $centromeres

cat $tmp_dir/__dbsnp_$$ | awk '{print $1, $2}' > $output
rm $tmp_dir/__dbsnp_$$

