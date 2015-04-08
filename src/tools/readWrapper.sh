chr=$1
begin=$2
end=$3

if [ "$4" = "ref" ]
then
	addr=/dupa-filer/laci/G1/
else
	addr=/dupa-filer/laci/I1/
fi

samtools view $addr/$chr/__plasma.part.bam $chr:$begin-$end

