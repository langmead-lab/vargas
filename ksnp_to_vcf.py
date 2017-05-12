#!/bin/python
'''
Ravi Gaddipati

Consumes a KSNP file (HISAT format) and outputs a VCF file.

Usage:
python3.5 ksnp_to_vcf.py <snps.ksnp> <output> <minpos> <maxpos> <chr>

snps.ksp : HISAT formatted SNP file
output : output file
minpos,maxpos : Only include SNPs within this range. All if both are 0
chr : chromosome
'''

import pprint as pp
import sys


def load_ksnps(filename):
    '''
    Map position to a tuple of REF and ALTs
    POS : (REF, [(ALT, AF), ...])
    '''
    alts = {}
    ids = {}
    with open(filename, 'r') as f:
        for line in f:
            s = line.split()
            if len(s) != 5:
                print("Invalid number of fields (8) in line:")
                print(line)
                exit(1)
            pos = int(s[3]) + 1 # ksnp is 0 based, vcf is 1 based
            ref = 'N'
            alt = s[4]
            af = 0.0
            if pos in alts:
                alts[pos][1].append(tuple((alt, af)))
                ids[pos].append(s[0])
            else:
                alts[pos] = tuple((ref, [tuple((alt,af))]))
                ids[pos] = [s[0]]
    return alts,ids

def to_vcf_record(ksnps, ids, pos, num_samples, chrom):
    line = chrom + '\t' + str(pos) + '\t'
    for x in ids[pos]:
        line += x + ';'
    line = line[0:-1] + '\t' + str(ksnps[pos][0]) + '\t'

    af_str = ""
    gt_str = "0|0\t"
    for i in range(len(ksnps[pos][1])):
        line += ksnps[pos][1][i][0]
        af_str += str(ksnps[pos][1][i][1])
        gt_str += str(i+1) + '|' + str(i+1)
        if i != len(ksnps[pos][1]) - 1:
            line += ','
            af_str += ','
            gt_str += '\t'
    for i in range(len(ksnps[pos][1]), num_haplotypes, 1):
        gt_str += '\t0|0'

    line += '\t' + '100' + '\t' + 'PASS' + '\t' + "AF=" + af_str + '\tGT\t' + gt_str

    return line 

def main():
    '''
    python ksnp_to_vcf.py snps.ksnp sortfile 1,2,4,8,16 output
    '''

    if (len(sys.argv) != 6):
        print("\npython ksnp_to_vcf.py <snps.ksnp> <output> <minpos> <maxpos> <chr:chrlen>\n\n"
              "snps.ksp : HISAT formatted SNP file\n"
              "output : output file\n"
              "minpos,maxpos : Only include SNPs within this range. All if both are 0\n"
              "chr : Chromosome\n")
        exit(1)

    ksnp_file = sys.argv[1]


    num_k = None
    minpos = None
    maxpos = None

    ksnps,ids = load_ksnps(ksnp_file)
    n_percent = float(len(ksnps)) / float(100)

    contiglens = []

    prefix = sys.argv[2]
    minpos = int(sys.argv[3])
    maxpos = int(sys.argv[4])
    sp = sys.argv[5].strip().split(',')
    for x in sp:
        contiglens.append(tuple(x.split(':')))

    max_alts = 0
    for k in ksnps:
        if len(ksnps[k][1]) > max_alts:
            max_alts = len(ksnps[k][1])

    
    vcf_header = "##fileformat=VCFv4\n##KSNPFILE=" + ksnp_file 
    for chpair in contiglens:
        vcf_header += "\n##contig=<ID=" + chpair[0] + ">\n"
    vcf_header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    vcf_header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
    vcf_header += "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">\n"
    vcf_header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"
    for i in range(max_alts):
        vcf_header += str(i)
        if i != max_alts - 1:
            vcf_header += '\t'
    vcf_header += '\n'


    keys = sorted([p for p in ksnps])

    with open(prefix, 'w') as o:
        o.write(vcf_header)
        for p in keys:
            o.write(to_vcf_record(ksnps, ids, p, max_alts, contiglens[0][0]) + '\n')



if __name__ == '__main__':
    main()