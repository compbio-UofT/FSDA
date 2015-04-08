# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random
import sys
import numpy as np
import pickle
from scipy.special import beta,betaln
import matplotlib.pyplot as plt

def logsumexp(v):
    tmp=np.array(np.copy(v))
    tmp_mx=np.max(tmp)
    tmp-=tmp_mx
    tmp=np.exp(tmp)
    tmp_mx+=np.log(np.sum(tmp))
    return tmp_mx

def logbetabin(ca,cb,r, rho, gen):
    qmap={"hom-het-alt":r/2, "hom-het-ref":1-r/2, "het-het":1.0/2, "het-hom-alt":1.0/2-r/2, "het-hom-ref":1.0/2+r/2}

    a=qmap[gen]
    b=1-a

    a*=(1-rho)/rho
    b*=(1-rho)/rho

    loc=betaln(ca+a,cb+b)
    loc-=betaln(a, b)

    return loc



def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('snps', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('out', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()


    snps=np.zeros((200,200))
    lsnps=[]
    for line in open(args.snps):
        tokens=line.strip().split("\t")
#        if ((int(tokens[2])+ int(tokens[3]))==40) and (float(tokens[2])/(int(tokens[2])+ int(tokens[3]))<0.2): ####################################### FILTER SNPS
#        if ((int(tokens[2])+ int(tokens[3]))>40) and (int(tokens[2])*int(tokens[3])>0) and (float(tokens[2])/(int(tokens[2])+ int(tokens[3]))<0.2): ####################################### FILTER SNPS
        if (int(tokens[2])+ int(tokens[3]))>40 and (int(tokens[2])*int(tokens[3])>0) and (tokens[7]!=tokens[8] or tokens[9]!=tokens[10]):
            snps[int(tokens[2])+int(tokens[3])][int(tokens[2])]+=1 #[[int(tokens[2]), int(tokens[3]), float(tokens[6])]]
            lsnps+=[[int(tokens[2]), int(tokens[3]), float(tokens[6])]]

    '''
    plt.hist([float(s[0])/(s[0]+s[1]) for s in lsnps], 50)
    plt.show()
    '''

    r=0.13
    pi={"hom-het-alt":0.15, "hom-het-ref":0.15, "het-het":0.35, "het-hom-alt":0.13, "het-hom-ref":0.22}
    rho=0.001

    mini=200
    minj=200
    maxi=0
    maxj=0
    for i in range(200):
        for j in range(200):
            if j>i:
                break
            if snps[i][j]>0:
                mini=min(i,mini)
                minj=min(j,minj)
                maxi=max(i,maxi)
                maxj=max(j,maxj)

    for step in range(100):
        resp=np.zeros((200,200,len(pi)))
        for i in range(mini,maxi):
            for j in range(minj,maxj):
                if j>i:
                    break
                for k,gen in enumerate(pi):
                    resp[i][j][k]=np.log(pi[gen])+logbetabin(j, i-j, r, rho, gen)
                resp[i][j]=np.exp(np.array(resp[i][j])-logsumexp(np.array(resp[i][j])))
        total=0
        newpi={"hom-het-alt":0., "hom-het-ref":0., "het-het":0., "het-hom-alt":0., "het-hom-ref":0.}
        for i in range(mini,maxi):
            for j in range(minj,maxj):
                if j>i:
                    break
                total+=snps[i][j]
                for k,gen in enumerate(pi):
                    newpi[gen]+=resp[i][j][k]*snps[i][j]
        for gen in newpi:
            newpi[gen]/=total
        pi["het-het"]=(newpi["het-het"]+newpi["het-hom-alt"]+newpi["het-hom-ref"])/2
        pi["het-hom-alt"]=(newpi["het-het"]+newpi["het-hom-alt"]+newpi["het-hom-ref"])/2   *   newpi["het-hom-alt"]/(newpi["het-hom-alt"]+newpi["het-hom-ref"])
        pi["het-hom-ref"]=(newpi["het-het"]+newpi["het-hom-alt"]+newpi["het-hom-ref"])/2   *   newpi["het-hom-ref"]/(newpi["het-hom-alt"]+newpi["het-hom-ref"])
        pi["hom-het-alt"]=newpi["hom-het-alt"]
        pi["hom-het-ref"]=newpi["hom-het-ref"]
        newr=0
        newrho=0
        totalr=0
        for k,gen in enumerate(pi):
            for i in range(mini,maxi):
                alpha=0
                beta=0
                mu1=0
                mu2=0
                tot=0
                for j in range(minj,maxj):
                    if j>i:
                        break
                    if snps[i][j]==0:
                        continue
                    mu1 += resp[i][j][k] *  j  * snps[i][j]
                    mu2 += resp[i][j][k] * j*j * snps[i][j]
                    tot += resp[i][j][k]
                #print tot
                if tot==0:
                    continue
                mu1/=tot
                mu2/=tot
                alpha = (tot*mu1 - mu2) / ( tot*(mu2/mu1- mu1 - 1.0) + mu1 ) 
                beta = ( tot - mu1 )*( tot - mu2/mu1 ) / ( tot*(mu2/mu1- mu1 - 1.0) + mu1 ) 
                if alpha+beta==0:
                    print alpha+beta, mu1, mu2, tot, i
                q=alpha/(alpha+beta)
                rmap={"hom-het-alt":2*q, "hom-het-ref":2-2*q, "het-hom-alt":1-2*q, "het-hom-ref":2*q-1}
                if gen is not "het-het":
                    newr+=tot*rmap[gen] 
                newrho+=1.0/(1+alpha+beta) * tot
                totalr+=tot
        r=newr/totalr
        newrho=newrho/totalr

    exit()
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

