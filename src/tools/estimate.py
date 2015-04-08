# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random
import sys
import math
from random import randint
from math import log 

def em(hist, mean, var, ratio):

    for it in range(20):
        n=len(mean)
        ct=[0]*n
        for l in range(len(hist)):
            p=[dist(l,mean[i],var[i])*ratio[i] for i in range(n)]
            if (sum(p)==0):
                p=[1.0]*n
            p=[x/sum(p) for x in p]

            ct=[ct[i]+hist[l]*p[i] for i in range(n)]

        ratio=[x/sum(hist) for x in ct]

    return mean, var, ratio

def normpdf(x, mean, var):
    pi = 3.1415926
    denom = (2*pi*var)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

def dist(x, mean, var):

    return normpdf(x, mean, var)
   # return foldednormpdf(x, mean, var)

def main():

    parser = argparse.ArgumentParser(description='salam')
    parser.add_argument('size_hist', type=str, nargs=1, help='paths to the file containing the histogram of fragment sizes')
    parser.add_argument('stats', type=str, nargs=1, help='paths to the file containing the histogram of fragment sizes')
    args = parser.parse_args()

    hist_file = open(args.size_hist[0], "r" )
    hist=[int(x) for x in hist_file]

    stats_file = open(args.stats[0], "r" )
    mean=map(float,stats_file.readline().strip().split(" "))
    var=map(float,stats_file.readline().strip().split(" "))

    #mean=[100+x*150/n for x in range(n)]
#    mean=[195, 160, 180, 170, 130, 140, 150, 160, 165,170, 180, 350]
    n=len(mean)
    '''
    mean.append(350)
    mean.append(350)
    n+=2
    '''
    ratio=[1.0/n]*n

    mean, var, ratio=em(hist, mean, var, ratio)
    '''
    ll, mean, var, ratio=em(hist, mean, var, ratio)
    for pp in range(10):
        bestM=mean
        bestV=var
        bestR=ratio
        bestL=ll
        for tries in range(50):
            m=mean[:]
            v=var[:]
            r=ratio[:]
            m.append(randint(70,180))
            v.append(10)
            r.append(1.0/len(m))
            r=[x/sum(r) for x in r]

            ll, m, v, r=em(hist, m, v, r)
            if ll>bestL:
                bestL=ll
                bestM=m[:]
                bestV=v[:]
                bestR=r[:]
        mean=bestM
        var=bestV
        ratio=bestR
        
    '''
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

