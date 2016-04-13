# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random
import sys

def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('dbsnp', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()

    dbsnp={}
    for line in open(args.dbsnp):
        tokens=line.strip().split("\t")
        dbsnp[tokens[0]+" "+tokens[1]]=(tokens[3],tokens[4]) #,tokens[5])


    for line in sys.stdin:
        if (line[0]=="#"):
            continue
        info=line.strip().split("\t")
        params=info[7].strip().split(";")
        for p in params:
            if p.startswith("DP4="):
                count=map(int,p[4:].strip().split(","))
                cA=count[0]+count[1]
                cB=count[2]+count[3]
                #if (min(cA,cB)>0): # and info[4]!="."): 
                dbsnpinfo=dbsnp[info[0]+" "+info[1]]
                if (info[3]==dbsnpinfo[0] and (info[4]==dbsnpinfo[1] or info[4]==".")):
                    print info[0]+"\t"+info[1]+"\t"+str(cA)+"\t"+str(cB)+"\t"+info[3]+"\t"+info[4]#+"\t"+dbsnpinfo[2]
                break

if __name__ == '__main__':
    main()

