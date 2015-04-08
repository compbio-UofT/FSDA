home="/home/aryan/fsda/feb/"

sample=$1
#readWrapper=/filer/aryan/fsda_files/feb/prep/readWrapper.sh

prep_adrs=$2
#prep_adrs=/filer/aryan/fsda_files/feb/prep/

echo "Find SNPs"
#./findAlleleCountsUsingGenotypes.sh $sample $prep_adrs/snplist.prp
#./findAlleleCounts.sh $sample $prep_adrs/snplist.prp

echo "Find Admixture ratio" ## needs fix
#./findAdmixtureRatio.sh $prep_adrs/snplist.prp $prep_adrs/admixtureRatio.prp

echo "Find fetal FSD" ## needs fix
#./findPureFetalFSD.sh $sample $prep_adrs/snplist.prp $prep_adrs/fetalDist.prp

echo "Create bins"
$home/src/tools/create_bins.sh $sample $prep_adrs/bins.pickle

echo "Find MoGs" ## needs fix
#./findMoG.sh $prep_adrs/bins.prp $prep_adrs/mog.prp

