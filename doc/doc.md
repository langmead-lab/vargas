[TOC]

_Updated: June 26, 2016_

Vargas aligns short reads to a directed acyclic graph (DAG). Reads are aligned using a SIMD vectorized version of Smith-Waterman, aligning multiple reads at once. The the aligner reports the best positions of reads without traceback to keep memory usage manageable. Additionally, vargas can be used to simulate reads from a DAG while stratifying by various parameters.

# Building {#building}

---

When cloning, use the `--recursive` option to automatically retrieve dependencies.

    git clone --recursive git@github.com:RaviGaddipati/vargas.git

## Dependencies {#dependencies}

Vargas relies on htslib to provide core file processing. Once cloned, the htslib is built with autoconf (version 2.63+).

    cd htslib
    autoconf && ./configure && make

Vargas also uses [libsimdpp](https://github.com/p12tic/libsimdpp) for SIMD support, and [doctest](https://github.com/onqtam/doctest).

## Compiling {#compile}

vargas is built with cmake.

    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release .. && make
    cd ..
    export PATH=${PWD}/bin:$PATH

# Modes of operation {#modes}

---

`vargas -h`

	---------------------- vargas, Jun 24 2016. rgaddip1@jhu.edu ----------------------
	Operating modes 'vargas MODE':
		define      Define a set of graphs for use with sim/align.
		sim         Simulate reads from a set of graphs.
		align       Align reads to a set of graphs.
		test        Run doctests.
		profile     Run profiles.

## define {#define}

`vargas define -h`

	-------------------- vargas define, Jun 24 2016. rgaddip1@jhu.edu --------------------
	-f	--fasta         *<string> Reference filename.
	-v	--var           *<string> VCF/BCF filename.
	-t	--out           *<string> Output filename.
	-g	--region        *<string> Region of graph, format CHR:MIN-MAX.
	-l	--nodelen       <int> Max node length, default 1,000,000
	-i	--ingroup       <int, int...> Percent ingroup subgraphs.
		                Prefix with 'c' to interpret as number of individuals.
	-n	--num           <int> Number of unique graphs to create for all -i. Default 1.

`define` is used to create a graph specification (*.gdf) for `sim` and `align` modes. Various subgraphs deriving off of a base graph are described. With `-i` and `-n` multiple subgraphs based on population filters are made. The original FASTA and variant files are needed in the same relative path when using the GDEF file.

For example:

    vargas define -f hs37d5.fa -v chr22.bcf -g 22:1,000,000-2,000,000 -i 15,30 -n 2 -t test.gdf

Will uniquely identify (including implicit graphs) 10 graphs:

- 2x 15% ingroup graphs
- 2x 30% ingroup graphs
- 2x 15% outgroup graphs
- 2x 30% outgroup graphs
- Reference graph
- Maximum allele frequency graph

The ingroup determines which samples are included in the graph. If a VCF file has 50 samples, an ingroup of up to 100 individuals can be made- each haplotype is treated independently.

## sim {#sim}

`vargas sim -h`

	-------------------- vargas sim, Jun 24 2016. rgaddip1@jhu.edu --------------------
	-g	--gdef          *<string> Graph definition file. Reads are simulated from Ingroups.
	-t	--outfile       *<string> Alignment output file.
	-o	--outgroup      Simulate from outgroup graphs.
	-n	--numreads      <int> Number of reads to simulate from each subgraph, default 1000.
	-m	--muterr        <float, float...> Read mutation error. Default 0.
	-i	--indelerr      <float, float...> Read indel error. Default 0.
	-v	--vnodes        <int, int...> Number of variant nodes, default any (-1).
	-b	--vbases        <int, int...> Number of variant bases, default any (-1).
	-l	--rlen          <int> Read length, default 50.
	-a	--rate          Interpret -m, -i as rates, instead of exact number of errors.
	-j	--threads       <int> Number of threads. 0 for maximum hardware concurrency.

	Comments preceded by '#'.
	-n reads are produced for each -b, -m, -i, -v combination.


`sim` uses a GDEF file and generates `-n` reads of each combination of `-m`, `-i`, `-v`, and `-b`. `-m` Introduces mutation errors, substituting n bases with an alternate base. Likewise, `-i` will delete a base or insert a random base. With `-a`, `-m` and `-i` are interpreted as rates (0 to 1). `-s` controls which subgraphs are used to generate reads. By default, all graphs in the GDEF file are used, and their complements (outgroup).

For example:

	vargas sim -g test.gdf -t reads -o -n 1000 -m 0,1,2 -v 0,1,2 -j 12

will generate 1000 reads for each combination of `-m`, `-v`, for each graph in `test.gdf` and their complements. This results in `3 -m * 3 -v * 1000 reads * 10 subgraphs = 90,000` reads.

## align {#align}

`vargas align -h`

	------------------- vargas align, Jun 24 2016. rgaddip1@jhu.edu -------------------
	-g	--gdef          *<string> Graph definition file.
	-r	--reads         *<string, string...> Read files to align.
	-t	--outfile       *<string> Alignment output file.
	-l	--rlen          <int> Max read length. Default 50.
	-R	                Align to reference graph.
	-X	                Align to maximum allele frequency graph.
	-I	                Align to ingroup graph.
	-O	                Align to outgroup graph.
	-m	--match         <int> Match score, default 2.
	-n	--mismatch      <int> Mismatch penalty, default 2.
	-o	--gap_open      <int> Gap opening penalty, default 3.
	-e	--gap_extend    <int> Gap extend penalty, default 1.
	-j	--threads       <int> Number of threads. 0 for maximum hardware concurrency.

	Lines beginning with '#' are ignored.


Reads are aligned to graphs specified in the GDEF file. `-R` Aligns the reads to the reference sequence, `-X` to the maximum allele frequency graph, `-I` to the ingroup graph, and `-O` to the outgroup graph. If a set of reads is generated from the ingroup and outgroup using `vargas sim -o`, a subsequent alignment using `vargas align -O` will align all reads to the same outgroup graph. A `[rxio]` tag at the beginning of the alignment record indicates which alignment was done.

For example:

	vargas align  -g test.gdf -r reads -t reads.aligns -R -I -j 12

Will align all the reads to both the reference graph and ingroup graphs, resulting in `2 * number_of_reads` alignments.

# File Formats {#formats}

---

## GDEF {#gdef}

GDEF Files are produced with `vargas define`. Below is an example produced by:

`vargas define -f hs37d5_22.fa -v chr22.bcf -r 22:22,000,000-23,000,000 -i c1,c5 -n 2`

	@gdef
	ref=hs37d5_22.fa;var=chr22.bcf;reg=22:22,000,000-23,000,000;nlen=1000000
	i,1,0,0:00010000000000000000000000000000000000
	i,1,1,0:00000000000000000000000000100000000000
	i,5,0,0:00010010000000000010001000100000000000
	i,5,1,0:00010000000110000000000000000001010000

Each gdef file begins with `@gdef`. The following line describes the source FASTA and VCF/BCF files, the region, and the maximum node length. Each following line describes an ingroup graph with the format:

`i,num,id,pct:filter`

where `num` indicates number of individuals if `pct=false` or percent when `pct=true`. The length of filter is equal to `num_samples * 2`

## Reads {#reads}

Reads from `vargas sim` are output in FASTA format.

	>end=22422171 mut=0 indel=0 vnode=-1 vbase=-1 src=i,1,0,0
	AAAAGGAAAAGAGATATCTGGACGCGCAGACAGAAAAAGATCCCAAAGAT
	>end=37431030 mut=0 indel=0 vnode=-1 vbase=-1 src=i,1,0,0
	ATGTTTCGTGTCCTGGTACCAGCACAGCAACTGACATGGATGGGGGGTTA
	>end=44240012 mut=0 indel=0 vnode=-1 vbase=-1 src=i,1,0,0
	GGGGGAGAAATCCCCAGCCTGTCCGCCACAGGAGGCCTCTGTGTTGACTC

Where `end` is the position of the last base in the read, `mut` is the number of mutation errors, `indel` is the number of indel erros, `vnode` is number of variant nodes traversed, and `vbase` is number of variant bases. `src` corresponds to the graph (as defined in the GDEF file, or its complement) the read came from. Parameters correspond to the profiles used to generate them, so `vnode=-1` can mean any number of variant nodes were traversed.

## Alignments {#aligns}

The output alignment format is a CSV line.

	r,o,25,0,1,TTGAGCTACNGAGTGGGGCACAGCNTCCCNCAGCCCAGGCTCCACGTACG,22453471,10,0,2,2,58,22453471,1,38,22453459,32,1
	r,o,25,1,1,NAGGAGCCCTATAAACNGCTTTATGTNGGCCTACTAAGCTGAGGGGCGAG,22761414,10,0,2,2,62,22761414,1,51,22761466,2,1
	i,i,25,0,0,ANAGGNGCAGTGCCTGCCTTNTCCGACTCACGTACCTCNGCCCNATGTGA,22554276,10,0,2,2,70,22554276,1,36,22554404,2,1
	i,i,25,1,0,TTTTGCTATTTTATGGGTCAGNGACCAATCATCATGACNATCAGATTTGA,22590921,10,0,2,2,64,22590921,1,53,22590983,1,1

The fields in the CSV line and there format:

- alignment type [rxio]
- read origin in/outgroup [io]
- read origin graph number [0-9]+
- read origin graph ID [0-9]+
- read sequence [ACGT]+
- read origin position [0-9]+
- read substitution errors [0-9]+
- read indel errors [0-9]+
- read variant nodes [0-9]+
- read variant bases [0-9]+
- best score [0-9]+
- best position [0-9]+
- best score count [0-9]+
- second best score [0-9]+
- second best position [0-9]+
- second best count [0-9]+
- read origin correlation flag [012]

Where the read origin correlation flag is 1 if the read origin position had the best score, 2 if it had the second best score, otherwise 0. 