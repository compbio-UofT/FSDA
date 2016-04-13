import argparse
import pickle
import sys

def dis(v,u):
    ks=0
    for i in range(len(v)):
        ks=max(ks,abs(v[i]-u[i]))
    return ks

def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('ref', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('test', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('dist', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('ks_th', type=float, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('reg', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()

    ref=pickle.load(open(args.ref , "rb" ))
    test=pickle.load(open(args.test , "rb" ) )
    dist=map(float, open(args.dist).readline().strip().split(" "))
    if (sum(dist)<1000):
        exit()
    dist=[x/sum(dist) for x in dist]
    for i in range(1,len(dist)):
        dist[i]+=dist[i-1]
    ctrls=[]
    cdf={}
    for d in ref:
        cdfd=[x/sum(ref[d]) for x in ref[d]]
        for i in range(1,len(cdfd)):
            cdfd[i]+=cdfd[i-1]
        cdf[d]=cdfd
        

    thresholds=[0]
    while thresholds[-1]<0.005:
        thresholds.append(thresholds[-1]+0.0005)
    counts={}
    for t in thresholds:
        counts[t]=[{'min_NN':1,'max_NN':4,'count':0},{'min_NN':5,'max_NN':9,'count':0},{'min_NN':10,'max_NN':100000,'count':0}]
    chert_ct=0
    for d1 in cdf:
        ks_list=[]
        for d2 in cdf:
            if d1!=d2:
                ks_list.append( dis(cdf[d1],cdf[d2]))
        NN={}
        for t in thresholds:
            NN[t]=0
        ks_list=sorted(ks_list)
        for x in ks_list:
            if x > thresholds[-1]:
                break
            for t in thresholds:
                if x<=t:
                    NN[t]+=1
        for t in thresholds:
            for bin_NN in counts[t]:
                if NN[t]<=bin_NN['max_NN'] and NN[t]>=bin_NN['min_NN']:
                    bin_NN['count']+=1
    #    if chert_ct==30:
    #        break
        chert_ct+=1
    for t in thresholds:
        sys.stdout.write(str(t)+"\t")
        for bin_NN in counts[t]:
            sys.stdout.write(str(bin_NN['count'])+"\t")
        print ""
    exit()







    for d in ref:
        ref_d=[x/sum(ref[d]) for x in ref[d]]
        for i in range(1,len(ref_d)):
            ref_d[i]+=ref_d[i-1]
        ctrls.append((d,dis(dist,ref_d)))
    ctrls=sorted(ctrls, key=lambda tup : tup[1])

    tot=0
    reg_split=args.reg.strip().split(" ")
    for i in range(len(ref)):
#        if (tot==20):
#            break
        if (ctrls[i][1]>args.ks_th):
            break
        ct_split=ctrls[i][0].strip().split(" ")
        if (ct_split[0]==reg_split[0]) and (min(abs(int(ct_split[2])-int(reg_split[1])),abs(int(ct_split[1])-int(reg_split[2])))<10000000):
            continue
        tot+=1
        for v in test[ctrls[i][0]]:
            sys.stdout.write(str(v)+" ")
        sys.stdout.write("\n")
    sys.stdout.flush()
    sys.stdout.close()


if __name__ == '__main__':
    main()

