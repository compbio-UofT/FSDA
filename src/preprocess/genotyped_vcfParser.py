# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random
import sys

def main():


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
                if (min(cA,cB)>=0): # and info[4]!="."): 
                    print info[0], info[1], cA, cB, info[3], info[4]
                break

if __name__ == '__main__':
    main()

