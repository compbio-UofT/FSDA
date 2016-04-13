#######################################################################################################################
## Cell Free DNA Fragment Size Distribution Analysis (FSDA) for non-invasive prenatal CNV prediction in fetal genome ##
## Aryan Arbabi - October 2014                                                                                       ##
## Computational Biology Lab, Department of Computer Science, University of Toronto                                  ## 
## arbabi@cs.toronto.edu                                                                                             ##
#######################################################################################################################

import argparse
import numpy.random
import sys
import math
import scipy.stats
import operator
import random

maxSize=500
   
#######################################################################################################################
#################################################### Main Function ####################################################
#######################################################################################################################

def main():
    ############## Parsing Arguments ##############
    parser = argparse.ArgumentParser(description='salam')
    parser.add_argument('hist', type=str, help='paths to the fetal dirichlet parameters')
    parser.add_argument('fetal_hist', type=str, help='paths to the fetal dirichlet parameters')
    parser.add_argument('r_real', type=float, help='paths to the fetal dirichlet parameters')
    parser.add_argument('r_goal', type=float, help='paths to the plasma dirichlet parameters')
    parser.add_argument('--c_real', type=float, help='paths to the plasma dirichlet parameters')
    parser.add_argument('--c_goal', type=float, help='paths to the plasma dirichlet parameters')
    parser.add_argument('--cnv', type=str, help='paths to the plasma dirichlet parameters')
    args = parser.parse_args()
    r_real=args.r_real
    r_goal=args.r_goal
    if args.c_real==None:
        c_real=1
        c_goal=1
    else:
        c_real=args.c_real
        c_goal=args.c_goal

    ############## Input files ##############
    hist=[]
    with open(args.hist, 'r') as f:
        for line in f:
            hist += [ numpy.array(map(float,line.strip().split(" "))) +0.00000001]
    with open(args.fetal_hist, 'r') as f:
        fetal_hist = numpy.array(map(float,f.readline().strip().split(" ")))
    fetal_dist=fetal_hist/sum(fetal_hist)
    plasma_dist=[ h/float(sum(h)) for h in hist ]

    fetal_remove_rate=1-(r_goal)*(1-r_real)/((r_real)*(1-r_goal))
    down_cov_rate=c_goal/(c_real-fetal_remove_rate*r_real*c_real)

    if (args.cnv=="del"):
        r_del=r_goal/(2-r_goal)
        fetal_remove_rate=1-(r_del)*(1-r_real)/((r_real)*(1-r_del))
#    remove_rate=fetal_remove_rate * r_real*fetal_dist/plasma_dist
    '''
    print "-->"
    print fetal_dist
    print fetal_remove_rate
    print r_goal
    print r_real
    '''
    remove_rate=[ fetal_remove_rate * r_real*fetal_dist/pd for pd in plasma_dist]
    #print remove_rate
    remove_rate=[ map( lambda x : fetal_remove_rate if x!=x or x>1 else x , rr) for rr in remove_rate]
    #print remove_rate
    #exit()

    if (down_cov_rate>1):
        down_cov_rate=1

    if (args.cnv=="dup"):
        preserve_rate=(1-fetal_remove_rate*r_real)*down_cov_rate
        dup=[numpy.random.multinomial(sum(h)*preserve_rate*r_goal/2, fetal_dist) for h in hist]

    for h_in, h in enumerate(hist):
        if h_in>0:
            sys.stdout.write("\n")
            #print ""
        for i in range(len(h)):
            if (h[i]==0):
                sys.stdout.write("0 ")
                #print 0,
                continue
            additional=0
            if (args.cnv=="dup"):
                additional=dup[h_in][i]
            sys.stdout.write(str(numpy.random.binomial(h[i],(1-remove_rate[h_in][i])*down_cov_rate)+additional) + str(" "))
            #print numpy.random.binomial(h[i],(1-remove_rate[h_in][i])*down_cov_rate)+additional,
    sys.stdout.flush()
    sys.stdout.close()

if __name__ == '__main__':
    main()

