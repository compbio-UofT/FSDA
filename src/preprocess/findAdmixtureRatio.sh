tmp_dir="/tmp/"
snps=$1
output=$2

#Rscript findAdmixtureRatio.r $snps $tmp_dir/__admixture_$$
python findAdmixtureRatio.py $snps $output

