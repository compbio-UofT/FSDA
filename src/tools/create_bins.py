import argparse
import pickle

def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('inp', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('out', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()

    name=True
    dists={}
    lastName=""
    for line in open(args.inp):
        if name:
            lastName=line.strip()
        else:
            dists[lastName]=map(float,line.strip().split(" "))
        name=not name
    pickle.dump(dists , open( args.out, "wb" ) )

if __name__ == '__main__':
    main()

