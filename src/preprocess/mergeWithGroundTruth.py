# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random
import sys

#grep -w -f <(cat /filer/aryan/fsda_files/feb/dbsnp/dbsnp_filtered_chr1.ref | awk '{print $2}') /dupa-filer/laci/I1/chr1/trio.phase.vcf | awk '{split($10,a,":"); printf "chr%s\t%s\t%d\t%d\t", $1, $2, substr(a[1],1,1), substr(a[1],3,1); split($12,a,":"); printf "%d\t%d\n", substr(a[1],1,1), substr(a[1],3,1); }' > groundTruth_I1

def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('genotypes', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('snplist', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()

    dbsnp={}
    for line in open(args.genotypes):
        if (line[0]=="#"):
            continue
        tokens=line.strip().split("\t")
        dbsnp["chr"+tokens[0]+" "+tokens[1]]=str(tokens[9][0])+"\t"+str(tokens[9][2])+"\t"+str(tokens[11][0])+"\t"+str(tokens[11][2])

    for line in open(args.snplist):
        info=line.strip().split("\t")
        key=info[0]+" "+info[1]
        if (key not in dbsnp):
            print line.strip()+"\t2\t2\t2\t2"
        else:
            print line.strip()+"\t"+dbsnp[key]

if __name__ == '__main__':
    main()

