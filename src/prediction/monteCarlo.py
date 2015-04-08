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
############################################### Printing weight vector ################################################
#######################################################################################################################

def printW(f,W):
    for r in range(W.shape[0]):
        for c in range(W.shape[1]):
            f.write(str(W[r,c])+' ')
        f.write('\n')

#######################################################################################################################
############################################ Dirichlet log likelihood PDF #############################################
#######################################################################################################################

def dirichlet_llhpdf(x, alpha):
    def gammaln(n):
        """Compute logarithm of Euler's gamma function for discrete values."""
        if n < 1:
            return float('inf')
        if n < 3:
            return 0.0
        c = [76.18009172947146, -86.50532032941677, \
             24.01409824083091, -1.231739572450155, \
             0.001208650973866179, -0.5395239384953 * 0.00001]
        x, y = float(n), float(n)
        tm = x + 5.5
        tm -= (x + 0.5) * math.log(tm)
        se = 1.0000000000000190015
        for j in range(6):
            y += 1.0
            se += c[j] / y
        return -tm + math.log(2.5066282746310005 * se / x)
    return (numpy.log(x)*(numpy.matrix(alpha)-1).T-(numpy.sum(map(gammaln,alpha))-gammaln(numpy.sum(alpha))))[0,0]

#######################################################################################################################
###################################################### Log sum ########################################################
#######################################################################################################################
def logsum(v):
    return numpy.max(v) + numpy.log( numpy.sum(numpy.exp(v-numpy.max(v))))

#######################################################################################################################
############################################# Fragment counts likelihood ##############################################
#######################################################################################################################

def prob(r, cc, Wpl, Wf, G, hist_test):
    cf=(cc-2.0)*r/(2.0+(cc-2.0)*r)
    Wtest= cf*Wf + (1-cf)*Wpl 
    #return numpy.log((Wtest*G)/numpy.sum((Wtest*G)))*hist_test.T
    #return numpy.log((Wtest*G)[:,90:250]/numpy.sum((Wtest*G)[:,90:250]))*hist_test[0,90:250].T
    #return numpy.log((Wtest*G)[:,90:250]/numpy.sum((Wtest*G)[:,90:250]))*hist_test[0,90:250].T
    tmp=((Wtest*G)[:,90:250])
    tmp/=tmp.sum(axis=1)
    return numpy.log(tmp)*hist_test[0,90:250].T

#######################################################################################################################
############################################### Copy number prediction ################################################
#######################################################################################################################

def findAns(MCtarget, Wpl, Wf, G, hist_test, hist_ctrl, r):
    M=numpy.zeros([Wpl.shape[0],1])

    ############## MonteCarlo sample likelihood ##############
    for i,x in enumerate(Wpl):
        p=numpy.array([dirichlet_llhpdf(x,MCtarget[j]) for j in range(len(MCtarget))])
        M[i,0]=logsum(p)

    ############## Ctrl likelihood ##############
#    print (Wpl*G)[:,90:250].shape
#    M+=numpy.log((Wpl*G)/numpy.sum((Wpl*G))) * hist_ctrl.T
    tmp=((Wpl*G)[:,90:250])
    tmp/=tmp.sum(axis=1)
    M+=numpy.log(tmp)*hist_ctrl[0,90:250].T

#    M+=numpy.log((Wpl*G)[:,90:250]/numpy.sum((Wpl*G)[:,90:250],axis=1)[:,numpy.newaxis]) * hist_ctrl[0,90:250].T
    
    ############## Test likelihoods ##############
    ans=[logsum(M+ prob(r, cc+1, Wpl, Wf, G, hist_test)) for cc in range(3)]

    ############## Normalizing probabilities ##############
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
    parser.add_argument('gstats', type=str, help='paths to the fetal dirichlet parameters')
    parser.add_argument("--dump", help="increase output verbosity", type=str)
    args = parser.parse_args()

    fetal_dist=numpy.matrix(map(float,open(args.fetal_dist).readline().strip().split(" ")))
    r=float(args.admixture) #pickle.load(open(args.admixture,"rb"))
    Hctrls=[]
    for c in open(args.ctrl_hists):
        Hctrls+=[ numpy.matrix([map(float,c.strip().split(" "))]) ]

    Htest=numpy.matrix(map(float,open(args.test_hist).readline().strip().split(" ")))

    ############## Creating G Matrix ##############
    gstats_file = open(args.gstats, "r")
    gstats={}
    gstats["m"]=map(float,gstats_file.readline().strip().split(" "))
    gstats["v"]=map(float,gstats_file.readline().strip().split(" "))
    Wctrls=numpy.matrix(map(float,gstats_file.readline().strip().split(" ")))

    G = numpy.zeros(shape=(len(gstats["m"]),maxSize))
    for i in range(len(gstats["m"])):
        for s in range(maxSize):
            G[i,s]=scipy.stats.norm(gstats["m"][i],numpy.sqrt(gstats["v"][i])).pdf(s)
    G=numpy.matrix(G)

    '''
    r=0
    for x in admixture:
        if (r not in admixture) or (admixture[r] < admixture[x]):
            r=x
    '''

    Dfetal=fetal_dist/numpy.sum(fetal_dist)
    Dtest=Htest/numpy.sum(Htest)
    ctrl_sum=numpy.sum(Hctrls,axis=0)

    ########################## MOG ##############################
    for it in range(200):
        normalizedG=numpy.multiply(G,Wctrls.T)
        normalizedG=normalizedG / numpy.sum(normalizedG,axis=0)
        Wctrls=(normalizedG*ctrl_sum.T).T/numpy.sum(ctrl_sum)


    '''
    for i in range(500):
        #print (Wctrls*G)[0,i]
        print (ctrl_sum)[0,i]
    exit()
    '''

    ########################## MOG ##############################

    Dctrl=ctrl_sum/numpy.sum(ctrl_sum)
    Dctrl=Wctrls*G

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
    

    ############## Main output ##############
    ans*=100
    max_index=numpy.argmax(ans)
    ans=["%0.2f" % x for x in ans]

    names=["Monosomy", "Normal", "Trisomy"]
    print "PREDICTION:\t",
    print names[max_index]+" "+ans[max_index]

    '''
    for i in range(3):
        print names[i]+" "+ans[i]
    '''


    exit()

    ############## Input files ##############
    mog_fetal_file = open(args.mog_fetal[0], "r")
    mog_ctrl_file = open(args.mog_ctrl[0], "r")
    mog_test_file = open(args.mog_test[0], "r")
    hist_ctrl_file = open(args.hist_ctrl[0], "r")
    hist_test_file = open(args.hist_test[0], "r")
    gstats_file = open(args.gstats[0], "r")

    ############## Creating G Matrix ##############
    gstats={}
    gstats["m"]=map(float,gstats_file.readline().strip().split(" "))
    gstats["v"]=map(float,gstats_file.readline().strip().split(" "))

    G = numpy.zeros(shape=(len(gstats["m"]),maxSize))
    for i in range(len(gstats["m"])):
        for s in range(maxSize):
            G[i,s]=scipy.stats.norm(gstats["m"][i],numpy.sqrt(gstats["v"][i])).pdf(s)
    G=numpy.matrix(G)
    
    ############## Reading Mixture of Gaussian Distributions for test, ctrl and pure fetal fragments ##############
    Wf=numpy.array(map(float,mog_fetal_file.readline().strip().split(" ")))
    Wctrl=numpy.array(map(float,mog_ctrl_file.readline().strip().split(" "))) 
    Wtest=numpy.array(map(float,mog_test_file.readline().strip().split(" ")))

    ############## Reading fragment size counts for test and ctrl regions ##############
    Hctrl=numpy.matrix(map(float,hist_ctrl_file.readline().strip().split(" "))) 
    Htest=numpy.matrix(map(float,hist_test_file.readline().strip().split(" "))) 

    ############## Reading fetal admixture ratio from arguments ##############
    r=args.admixture

    ############## Creating monte carlo sample matrix ##############
    cof_alpha=10000.0

    MCtarget=[((2+(cc-2)*r)/2.0 * Wtest - (cc-2)*r/2.0 * Wf)*cof_alpha for cc in range(1,4)]
    MCtarget+=[((2+(cc-2)*r)/2.0 * Wctrl - (cc-2)*r/2.0 * Wf)*cof_alpha for cc in range(1,4)]
    MCtarget+=[((2+(cc-2)*r)/4.0 * Wtest - (cc-2)*r/4.0 * Wf + Wctrl/2.0)*cof_alpha for cc in range(1,4)]
    #MCtarget+=[Wctrl*cof_alpha]
    '''
    print (MCtarget[0]*G)[0,1:10]
    print numpy.asarray((MCtarget[0]*G)[0,1:10]).flatten()
    plt.plot(numpy.asarray(MCtarget[0]*G/cof_alpha).flatten())
    plt.plot(numpy.asarray(MCtarget[1]*G/cof_alpha).flatten())
    plt.plot(numpy.asarray(MCtarget[2]*G/cof_alpha).flatten())
    plt.plot(numpy.asarray(Hctrl/numpy.sum(Hctrl)).flatten())
    plt.plot(numpy.asarray(Htest/numpy.sum(Htest)).flatten())
    plt.show()
    '''
    #Wpl=numpy.vstack([numpy.random.dirichlet(p,1000) for p in MCtarget])

    ############## Predicting the copy number and printing ##############
    Dctrl=Hctrl/numpy.sum(Hctrl)
    Dctrl=Dctrl/numpy.sum(Dctrl)
    Dctrl=Wctrl*G
    ans=[0,0,0]
    '''
    for i in range(Dctrl.shape[1]):
        print Dctrl[0,i]
    exit()
    '''
    for cc in range(1,4):
        cf=(cc-2.0)*r/(2.0+(cc-2.0)*r)
        tmp=(1-cf)*Dctrl+cf*Wf*G
        #printW(open("__tmp_cc"+str(cc),"w"),tmp)
        tmp=tmp[0,90:250]
        if numpy.min(tmp)<=0:
            tmp-=numpy.min(tmp)-0.00000001
        tmp/=numpy.sum(tmp)
        ans[cc-1]=numpy.log(tmp) * Htest[0,90:250].T
        #ans[cc-1]=numpy.log(tmp) * Htest.T
    ans=[(x - max(ans)) for x in ans]
    ans=numpy.exp(ans)
    ans/=numpy.sum(ans)
    #ans=findAns(MCtarget, Wpl, Wf, G, Htest, Hctrl, r)
    

    ############## Main output ##############
    ans*=100
    max_index=numpy.argmax(ans)
    ans=["%0.2f" % x for x in ans]

    names=["Monosomy", "Normal", "Trisomy"]
    print "PREDICTION:\t",
    print names[max_index]+" "+ans[max_index]

    for i in range(3):
        print names[i]+" "+ans[i]


    ############## Additional output ##############
    if args.dump!=None:
        f = open(args.dump,'w')
        printW(f,Wf*G)
        f.close()

if __name__ == '__main__':
    main()

