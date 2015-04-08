# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random
import sys
import math
from random import randint
from math import log 

def em(hist, mean, var, ratio, num):
    n=len(mean)
    for iteration in range(num):
        newMean=[0]*n
        newVar=[0]*n
        newRatio=[0]*n
        ct=[0]*n
        for l in range(len(hist)):
            p=[dist(l,mean[i],var[i])*ratio[i] for i in range(n)]
            if (sum(p)==0):
                p=[1.0]*n
            p=[x/sum(p) for x in p]

            newMean=[newMean[i]+hist[l]*l*p[i] for i in range(n)]
            ct=[ct[i]+hist[l]*p[i] for i in range(n)]

        newMean=[newMean[i]/ct[i] for i in range(n)]
        newRatio=[x/sum(hist) for x in ct]

        for l in range(len(hist)):
            p=[dist(l,newMean[i],var[i])*ratio[i] for i in range(n)]
            if (sum(p)==0):
                p=[1.0]*n
            p=[x/sum(p) for x in p]

            newVar=[newVar[i]+hist[l]*(l-newMean[i])*(l-newMean[i])*p[i] for i in range(n)]
            for i in range(len(newVar)):
                if newVar[i]<0.001:
                    newVar[i]=0.001

        newVar=[newVar[i]/ct[i] for i in range(n)]
#        newVar=[sum(newVar)/sum(ct)]*n

        var=newVar
        mean=newMean
        ratio=newRatio

    likelihood=0
    for l in range(len(hist)):
        tmp=0
        for i in range(n):
            tmp+=dist(l,mean[i],var[i])*ratio[i]
        likelihood+=hist[l]*math.log(tmp)

    return likelihood, mean, var, ratio

def normpdf(x, mean, var):
    pi = 3.1415926
    denom = (2*pi*var)**.5
    try:
        num = math.exp(-(float(x)-float(mean))**2/(2*var))
    except:
        sys.stderr.write("Error: ")
        sys.stderr.write("var: "+str(var)+"\n")
        
    return num/denom

def merge(mean, var, ratio):
    m=(mean[0]*ratio[0] + mean[1]*ratio[1])/sum(ratio)
    v=( (var[0]+mean[0]**2)*ratio[0] + (var[1]+mean[1]**2)*ratio[1] )/(sum(ratio)) - m**2
    if (v==0):
        sys.stderr.write("holy $#!\n")
        sys.stderr.write(mean)
        sys.stderr.write("\n")
        sys.stderr.write(var)
        sys.stderr.write("\n")
        sys.stderr.write(ratio)
        sys.stderr.write("\n")
    return m, v, sum(ratio)

def smalize(hist, mean, var, ratio, num):
    for pp in range(num):
        bestM=mean
        bestV=var
        bestR=ratio
        bestL=-2**31
        sys.stderr.write("p: "+str(pp)+"\n")
        for i in range(len(mean)):
            for j in range(i,len(mean)):
                m=mean[:]
                v=var[:]
                r=ratio[:]
#                if abs(m[i]-m[j]) > abs(

                m[i], v[i], r[i] = merge([m[i], m[j]], [v[i], v[j]], [r[i], r[j]])
                m.pop(j)
                v.pop(j)
                r.pop(j)
                if (i==j):
                    r=[x/sum(r) for x in r]

                ll, m, v, r=em(hist, m, v, r,10)
                if ll>bestL:
                    bestL=ll
                    bestM=m[:]
                    bestV=v[:]
                    bestR=r[:]
        mean=bestM[:]
        var=bestV[:]
        ratio=bestR[:]
        ll, mean, var, ratio=em(hist, mean, var, ratio, 500)
    return ll, mean, var, ratio

def dist(x, mean, var):

    return normpdf(x, mean, var)
   # return foldednormpdf(x, mean, var)

def main():

    parser = argparse.ArgumentParser(description='salam')
    parser.add_argument('size_hist', type=str, nargs=1, help='paths to the file containing the histogram of fragment sizes')
    parser.add_argument('gaussians_count', type=int, nargs=1, help='the number of gaussians')
    args = parser.parse_args()

    hist_file = open(args.size_hist[0], "r" )
#    hist=[int(x) for x in hist_file]
    hist=map(float, hist_file.readline().strip().split(" "))

    n=args.gaussians_count[0]+10

    n-=4
    mean=[80+x*150/n for x in range(n)]
    mean.append(350)
    mean.append(340)
    mean.append(370)
    mean.append(320)
    n+=4
    n=len(mean)

    batch=2

    var=[30]*n
    ratio=[1.0/n]*n

    '''
    ll, mean, var, ratio=em(hist, mean, var, ratio, 500)
    for pp in range(1):
        bestM=mean
        bestV=var
        bestR=ratio
        bestL=ll
        for tries in range(100):#30**batch):
            m=mean[:]
            v=var[:]
            r=ratio[:]
            for b in range(batch):
                m.append(randint(120,180))
                v.append(30)
            for b in range(batch):
                r.append(1.0/len(m))
            r=[x/sum(r) for x in r]

            ll, m, v, r=em(hist, m, v, r,30)
            if ll>bestL:
                bestL=ll
                bestM=m[:]
                bestV=v[:]
                bestR=r[:]
        mean=bestM[:]
        var=bestV[:]
        ratio=bestR[:]
        ll, mean, var, ratio=em(hist, mean, var, ratio, 500)
        
    '''
    ll, mean, var, ratio=em(hist, mean, var, ratio, 2000)
    ll, mean, var, ratio=smalize(hist, mean, var, ratio, 10)
        #print ll

    for x in mean:
        print x,
    print ""

    for x in var:
        print x,
    print ""

    for x in ratio:
        print x,
    print ""


if __name__ == '__main__':
    main()

