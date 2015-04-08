import argparse
import pickle
import matplotlib.pyplot as plt
import numpy

def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('snps', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    #parser.add_argument('out', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()
   
    snps=[]
    for snp in open(args.snps):
        tmp=snp.strip().split()
        snps+=[(int(tmp[2]),int(tmp[3]))]
    ratios=numpy.array([float(x[0])/(x[0]+x[1]) for x in snps])
    print ratios
    plt.hist(ratios,bins=40)
    plt.title("Allele ratio Histogram")
    plt.xlabel("ratio")
    plt.ylabel("Frequency")
    plt.show()
    '''
    print snps
    print len(snps)
    '''
    exit()

    r=map(float, inpfile.readline().strip().split(" "))
    p=map(float, inpfile.readline().strip().split(" "))
    admixture={}
    for i in range(len(r)):
        admixture[r[i]]=p[i]/sum(p)
    '''
    for x in sorted(admixture):
        print x, admixture[x]
    '''

    pickle.dump(admixture, open(args.out , "wb" ) )

if __name__ == '__main__':
    main()

