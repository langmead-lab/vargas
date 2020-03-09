
# Usage: python vargas_preprocess_VCF.py varfile.vcf
# Python 3.7

# prints header + records to stdout in VCF format with comma separated multi-alt-alleles where necessary (remains sorted)

# The input VCF must be:
#  - biallelic
#  - only contain SNPs and INDELs
#  - all INDELs have at least 1 base in REF and ALT fields - from VCF 4.2 specification: "For simple insertions and
#       deletions in which either the REF or one of the ALT alleles would otherwise be null/empty, the REF and ALT
#       Strings must include the base before the event (which must be reflected in the POS field)"
#  - sorted by position (with records from each contig continuous)

import sys
from itertools import combinations

vcf = open(sys.argv[1], 'r')

curr_chrom = None
overlapping_lines = []
overlap_end = 0


def variants_are_compatible(L):
	for (x,y) in combinations(L,2): 
		x_start = int(x[1])
		x_end = x_start + len(x[3]) - 1
		y_start = int(y[1])
		y_end = y_start + len(y[3]) - 1
		
		x_snp = len(x[3]) == 1 and len(x[4]) == 1
		y_snp = len(y[3]) == 1 and len(y[4]) == 1

		x_del = len(x[3]) > 1
		y_del = len(y[3]) > 1

		x_ins = len(x[4]) > 1
		y_ins = len(y[4]) > 1

		# variants of the same type at the same position
		if ((x_snp and y_snp) or (x_del and y_del) or (x_ins and y_ins)) and x_start == y_start: return False
		# deletions overlap		
		# if the last base of the first deletion is the DEL base of the second
		# they are still compatible, hence "y_start < x_end"
		if x_del and y_del and y_start < x_end: return False
		# the SNP base is deleted 
		# note, the first base of a deletion is not actually deleted, hence "y_start > x_start"
		if x_del and y_snp and y_start > x_start and y_start <= x_end: return False
		# INS base is deleted 
		# note, the first base of a deletion region is not deleted, hence "y_start > x_start"
		# if the last base of the deletion is the INS base, they are still compatible, hence "y_start < x_end"
		if x_del and y_ins and y_start > x_start and y_start < x_end: return False
		# for the 2nd and 3rd cases, the start position of y is always >= the start position of x
		# so we don't need to check the opposite direction or the relative start positions.


	return True


def process_overlapping_lines(L):
	if len(L) == 0: return # no variants = nothing to do 
	elif len(L) == 1: # one variant = no overlap = don't change the record
		print("\t".join(L[0]))
	else: # build a multi-alt-allele consensus record for all variants, in all compatible combinations
		overlap_start = int(L[0][1]) # records are in positional order, so the first start position is the start of the interval
		
		# build the ref allele
		refallele = L[0][3] # the prefix
		ref_endpos = overlap_start + len(refallele) - 1
		for rec in L[1:]:
			rec_endpos = int(rec[1]) + len(rec[3]) - 1
			if rec_endpos > ref_endpos: # if the ref allele of this continues beyond the ref allele we have, add some characters
				offset = ref_endpos - int(rec[1]) + 1 # how many characters to skip
				refallele = refallele + rec[3][offset:]
				ref_endpos = rec_endpos

		# assert(this_endpos = overlap_end)
		# assert(len(refallele) == overlap_end - overlap_start + 1)

		# build the alt alleles

		# alt alleles with one variant
		altallele_list = []
		for rec in L:
			rec_startpos = int(rec[1])
			rec_endpos = rec_startpos + len(rec[3]) - 1
			prefix_length = rec_startpos - overlap_start # length of prefix from reference allele
			suffix_length = overlap_end - rec_endpos # length of suffix from reference allele
			if suffix_length == 0:
				a = refallele[:prefix_length] + rec[4]
			else:
				a = refallele[:prefix_length] + rec[4] + refallele[-suffix_length:] 
			altallele_list.append(a)

		# alt alleles with >1 variant
		for r in range(2,len(L)+1):
			for com in combinations(L,r):
				if variants_are_compatible(com):
					# build the alt allele with all these variants combined
					# each character is a list element
					a = list(refallele)
					# add the SNPs before the length changes
					for rec in com:
						if len(rec[3]) == 1 and len(rec[4]) == 1: # is a SNP
							rec_startpos = int(rec[1])
							offset = rec_startpos - overlap_start
							a[offset] = rec[4]

					# add the indels from right-to-left
					last_offset_insertion = -1
					for rec in reversed(com):
						if len(rec[3]) == 1 and len(rec[4]) > 1: # is an insertion
							rec_startpos = int(rec[1])
							offset = rec_startpos - overlap_start
							a.insert(offset+1, rec[4][1:]) #element is inserted before the index supplied to insert()
							last_offset_insertion = offset
						if len(rec[3]) > 1 and len(rec[4]) == 1: # is a deletion
							rec_startpos = int(rec[1])
							offset = rec_startpos - overlap_start
							if last_offset_insertion == offset: offset += 1 # there was an insertion at the same postion and it has already been implemented, so need to adjust which positions will be deleted
							for i in range(len(rec[3])-1):
								a.pop(offset+1) #element at the index supplied to pop() is deleted

					a = ''.join(a)
					altallele_list.append(a)

		# remove possible duplicate sequences
		altallele_list = set(altallele_list)

		# print the consensus record
		print(curr_chrom + "\t" + 
			  str(overlap_start) + "\t.\t" +
			  refallele + "\t" +
			  ",".join(altallele_list) + "\t.\tPASS\t.\t.")


for line in vcf.readlines():
	if line[0] == "#": # header line
		print(line, end='') # don't include an additional newline
		continue
	fields = line.strip().split() # VCF is tab separated
	if fields[0] != curr_chrom:
		if curr_chrom is not None: 
			process_overlapping_lines(overlapping_lines)
			overlapping_lines = []
			overlap_end = 0
		curr_chrom = fields[0]
	if int(fields[1]) <= overlap_end: # start position of this variant is within the current overlap interval -> add to set to be processed
		overlapping_lines.append(fields)
		overlap_end = max(overlap_end, int(fields[1]) + len(fields[3]) - 1) # extend interval bounds
	else:
		process_overlapping_lines(overlapping_lines)
		overlapping_lines = [fields]
		overlap_end = int(fields[1]) + len(fields[3]) - 1 # reset interval bounds

process_overlapping_lines(overlapping_lines)
