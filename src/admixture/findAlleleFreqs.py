import argparse

def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('snps', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('population', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()

    pop=[]
    for sample in open(args.population):
        pop.append(sample.strip())
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
            continue
        token=line.strip().split("\t")
        refs=0
        for i in pop_inds:
            if token[i][0]=='0':
                refs+=1
            if token[i][2]=='0':
                refs+=1
        print float(refs)/len(pop_inds)/2
        

    exit()
   
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

