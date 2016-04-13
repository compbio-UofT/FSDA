tmp_dir=$TMPDIR
if [ "$tmp_dir" == "" ]; then
	tmp_dir="/tmp/"
fi
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#source $DIR/../config_file.sh


config_hg19=/filer/aryan/cheo/data/hg19.fa
snplist=/dupa-filer/aryan/last_ref_fsda/dbSNP_snplist


sample=$1
ref=$2
prep=$3

mkdir -p $prep/allele_counts/

########################################################################################################################################
########################################################## Find allele counts ##########################################################
########################################################################################################################################


for file in $sample/*.bam
do
	filename=`basename $file`
	filename=${filename%.bam}
	samtools mpileup -l <(awk '{print $1, $2}' $snplist) -f $config_hg19 -t DPR,DP -uIB  $file  -Q10 -q10 | bcftools call - -Am | pypy $DIR/vcfParser.py $snplist > $prep/allele_counts/$filename".allele_counts" #&
done


########################################################################################################################################
####################################################### Find fetal distribution ########################################################
########################################################################################################################################

$DIR/findPureFetalFSD.sh $sample $prep > $prep/pure_fetal_hist

########################################################################################################################################
####################################################### Create distribution bins #######################################################
########################################################################################################################################

$DIR/../tools/create_bins.sh $sample $prep/sample_bins.pickle

$DIR/../tools/create_bins.sh $ref $prep/ref_bins.pickle
