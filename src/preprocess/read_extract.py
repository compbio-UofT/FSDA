# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random

ALL_PROPERLY_ALIGNED = 0x2
UNMAPPED = 0x4
sourceIndex={'IM1':9,'IP1':10,'IC1':11}

def mapping_parser(m):
    '''
    Parse a read in SAM format, return a dictionary with filled in fields of interest.
    '''
    if isinstance(m, str):
        m = m.strip().split('\t')
        d = {}
        d['flag'] = int(m[1])   # flags
        d['chr'] = m[2]         # chr
        d['pos'] = int(m[3])    # pos
        d['mapq'] = int(m[4])   # mapping quality
        d['cigar'] = m[5]       # cigar string
        d['length'] = m[8]       # fragment length
        d['seq'] = m[9]         # sequence
        d['qual'] = m[10]       # sequencing quality

    return d

def parse_cigar_string(s):
    '''
    Parse given CIGAR string to a list of operators.
    '''
    if s[0] == '*': return []
    
    res = []
    crt_len = ''
    i = 0
    while i < len(s):
        if str.isdigit(s[i]):
            crt_len += s[i]
        else:
            res.append([s[i], int(crt_len)])
            crt_len = ''
        i += 1  
    return res

def upper_bound(arr, ch, pos, q, p):

    if (p-q==1):
        if arr[q]['chr']==ch and arr[q]['pos']==pos:
            return q
        else:
            return p
    
    mid=(q+p)/2;

    if arr[mid]['chr']>ch:
        return upper_bound(arr, ch, pos, q, mid);
    elif arr[mid]['chr']<ch:
        return upper_bound(arr, ch, pos, mid, p);

    if arr[mid]['pos']>pos:
        return upper_bound(arr, ch, pos, q, mid);
    else:
        return upper_bound(arr, ch, pos, mid, p);
 
def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('read_file', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('pos', type=int, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    parser.add_argument('allele', type=str, help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()

    reads_file=open(args.read_file, "r")
    pos=args.pos
    allele=args.allele

    # For each read count the snips in the read and check if it belongs to the target hapoltype:
    for read in reads_file:
        parsed_read= mapping_parser(read)

        # If the read is not aligned properly, ignore it
        #if parsed_read['flag'] & ALL_PROPERLY_ALIGNED == 0: continue
        if parsed_read['flag'] & UNMAPPED != 0: continue 
        if parsed_read['cigar'] == '*': continue

        pos_qr = 0
        pos_db = parsed_read['pos']
        op = parse_cigar_string(parsed_read['cigar'])

        # Count the snips for each haplotype using the cigar string
        for o in op:
            if o[0] == 'H': continue
            elif o[0] in 'SI': pos_qr += o[1]
            elif o[0] in 'ND': pos_db += o[1]
            elif o[0] in 'M=X':
                for i in range(o[1]):
                    if pos_db==pos: #and ord(parsed_read['qual'][pos_qr]) >= 33+20:
                        if allele==parsed_read['seq'][pos_qr]:
                            print parsed_read['length']
                            pos_db += 1
                            pos_qr += 1
                            break
                    pos_db += 1
                    pos_qr += 1

if __name__ == '__main__':
    main()

