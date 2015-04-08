import argparse
import pickle
import sys

def dis(v,u):
    ks=0
    for i in range(len(v)):
        ks=max(ks,abs(v[i]-u[i]))
    return ks

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
    for d in ref:
        ref_d=[x/sum(ref[d]) for x in ref[d]]
        for i in range(1,len(ref_d)):
            ref_d[i]+=ref_d[i-1]
        ctrls.append((d,dis(dist,ref_d)))
    ctrls=sorted(ctrls, key=lambda tup : tup[1])

    tot=0
    reg_split=args.reg.strip().split(" ")
    for i in range(len(ref)):
        if (tot==20):
            break
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

