# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random
import sys
import numpy as np
import pickle
from scipy.special import beta,betaln
#import matplotlib.pyplot as plt

def logsumexp(v):
    tmp=np.array(np.copy(v))
    tmp_mx=np.max(tmp)
    tmp-=tmp_mx
    tmp=np.exp(tmp)
    tmp_mx+=np.log(np.sum(tmp))
    return tmp_mx

def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('snps', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('out', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()


    '''
    '''

    r=[0.11,  0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15]
    p=[0.0625,  0.0625, 0.0625, 0.0625, 0.5, 0.0625, 0.0625, 0.0625, 0.0625]

    admixture={}

    for i in range(len(r)):
        admixture[r[i]]=p[i]#/sum(p)

    pickle.dump(admixture, open(args.out , "wb" ) )
    exit()

    '''
    '''

    snps=[]
    for line in open(args.snps):
        tokens=line.strip().split("\t")
        if ((int(tokens[2])+ int(tokens[3]))==40) and (float(tokens[2])/(int(tokens[2])+ int(tokens[3]))<0.2): ####################################### FILTER SNPS
        #if ((int(tokens[2])+ int(tokens[3]))>40) and (int(tokens[2])*int(tokens[3])>0) and (float(tokens[2])/(int(tokens[2])+ int(tokens[3]))<0.2): ####################################### FILTER SNPS
            snps+=[[int(tokens[2]), int(tokens[3]), float(tokens[6])]]
    print(len(snps))
    '''
    plt.hist([float(s[0])/(s[0]+s[1]) for s in snps], 20)
    plt.show()
    '''
    gtypes=[("00","00"),("00","01"),("11","11"),("11","01"),("01","00"),("01","01"),("01","11")]

    gpriors=np.zeros(shape=(len(snps),len(gtypes)))
    lref=np.log(np.array(snps)[:,2])
    lalt=np.log(1-np.array(snps)[:,2])
    for j,g in enumerate(gtypes):
        if g[0]=="00":
            gpriors[:,j]+=2*lref
            if g[1]=="00":
                gpriors[:,j]+=lref
            if g[1]=="01":
                gpriors[:,j]+=lalt

        if g[0]=="11":
            gpriors[:,j]+=2*lalt
            if g[1]=="11":
                gpriors[:,j]+=lalt
            if g[1]=="01":
                gpriors[:,j]+=lref

        if g[0]=="01":
            gpriors[:,j]+=lalt+lref
            if g[1]=="00":
                gpriors[:,j]+=lref
            if g[1]=="11":
                gpriors[:,j]+=lalt
            if g[1]=="01":
                gpriors[:,j]+=0

    last=0
    for r in np.arange(0.10,0.2,0.005):
        rllh=[]
        for eps in np.arange(0.00,0.005,0.001):
            for rho in np.arange(0.005,0.02,0.005):
                llh=np.copy(gpriors) 
                #llh=np.zeros(shape=gpriors.shape)
                for j,g in enumerate(gtypes):
                    #eps=0.001
                    if g[0]=="00":
                        if g[1]=="00":
                            a=1-eps
                            b=0+eps
                        if g[1]=="01":
                            a=1-r/2
                            b=r/2
                    if g[0]=="11":
                        if g[1]=="11":
                            a=0+eps
                            b=1-eps
                        if g[1]=="01":
                            a=r/2
                            b=1-r/2
                    if g[0]=="01":
                        if g[1]=="00":
                            a=1.0/2+r/2
                            b=(1.0-r)/2
                        if g[1]=="11":
                            a=(1.0-r)/2
                            b=1.0/2+r/2
                        if g[1]=="01":
                            a=1.0/2
                            b=1.0/2
                    a*=(1-rho)/rho
                    b*=(1-rho)/rho
                    llh[:,j]-=betaln(a, b)
                    for i,s in enumerate(snps):
                        llh[i,j]+=betaln(s[0]+a,s[1]+b)

                llh_maxpr=np.amax(llh, axis=1)
                llh-=np.matrix(llh_maxpr).T
                llh=np.exp(llh)
                llh_maxpr+=np.log(np.sum(llh,axis=1))
                rllh+=[np.sum(llh_maxpr)]
        tmp=logsumexp(rllh)
        print r, tmp-last
        last=tmp
    exit()

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
                    print info[0]+"\t"+info[1]+"\t"+str(cA)+"\t"+str(cB)+"\t"+info[3]+"\t"+info[4]+"\t"+dbsnpinfo[2]
                break

if __name__ == '__main__':
    main()

