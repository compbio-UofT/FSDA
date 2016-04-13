import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('snps', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('population', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()

    pop=[]
    for sample in open(args.population):
        pop.append(sample.strip().split("\t")[0])
    pop=set(pop)

    for line in open(args.snps):
        if (line[1]=="#"):
            continue
        elif (line[0]=="#"):
            header=line.strip().split("\t")
            pop_inds=[]
            for i in range(9,len(header)):
                if header[i] in pop:
                    pop_inds.append(i)
            print (len(pop_inds))
            continue
        token=line.strip().split("\t")
        if len(token[3])!=1 or len(token[4])!=1:
            continue
        refs=0
        for i in pop_inds:
            if token[i][0]=='0':
                refs+=1
            if token[i][2]=='0':
                refs+=1
        print "chr"+token[0]+"\t"+token[1]+"\t"+token[2]+"\t"+token[3]+"\t"+token[4]+"\t"+str(float(refs)/len(pop_inds)/2)

if __name__ == '__main__':
    main()

