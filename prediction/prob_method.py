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
import pickle
#import matplotlib.pyplot as plt

maxSize=500

#######################################################################################################################
############################################### Copy number prediction ################################################
#######################################################################################################################

def predict(Dctrl, Dfetal, Htest, r):
    ans=[0,0,0]
    for cc in range(1,4):
        cf=(cc-2.0)*r/(2.0+(cc-2.0)*r)
        tmp=(1-cf)*Dctrl+cf*Dfetal
        tmp=tmp[0,90:250]
        if numpy.min(tmp)<=0:
            tmp-=numpy.min(tmp)-0.00000001
        tmp/=numpy.sum(tmp)
        ans[cc-1]=numpy.log(tmp) * Htest[0,90:250].T
        #ans[cc-1]=numpy.log(tmp) * Htest.T
    ans=[(x - max(ans)) for x in ans]
    ans=numpy.exp(ans)
    ans/=numpy.sum(ans)
    return ans
   
#######################################################################################################################
#################################################### Main Function ####################################################
#######################################################################################################################

def main():
    ############## Parsing Arguments ##############
    parser = argparse.ArgumentParser(description='salam')
    parser.add_argument('fetal_dist', type=str, help='paths to the fetal dirichlet parameters')
    parser.add_argument('admixture', type=str, help='paths to the plasma dirichlet parameters')
    parser.add_argument('ctrl_hists', type=str,  help='paths to the plasma dirichlet parameters')
    parser.add_argument('test_hist', type=str, help='paths to the fetal dirichlet parameters')
    args = parser.parse_args()

    fetal_dist=numpy.matrix(map(float,open(args.fetal_dist).readline().strip().split(" ")))
    r=float(args.admixture) #pickle.load(open(args.admixture,"rb"))
    Hctrls=[]
    for c in open(args.ctrl_hists):
        Hctrls+=[ numpy.matrix([map(float,c.strip().split(" "))]) ]

    Htest=numpy.matrix(map(float,open(args.test_hist).readline().strip().split(" ")))

    Dfetal=fetal_dist/numpy.sum(fetal_dist)
    ctrl_sum=numpy.sum(Hctrls,axis=0)

	Hctrls=[ctrl_sum]

    names=["Monosomy", "Normal", "Trisomy"]
    anscount=[0,0,0]
    
    for h in Hctrls:

        Dctrl=h/numpy.sum(h)

        ans=predict(Dctrl, Dfetal, Htest, r)

        ans*=100
        max_index=numpy.argmax(ans)
        ans=["%0.2f" % x for x in ans]

        anscount[max_index]+=1

    best_ans=numpy.argmax(anscount)
    print "PREDICTION:\t",
    print names[best_ans]+" "+str(anscount[best_ans])


if __name__ == '__main__':
    main()

