[TOC]

# Variant Genome Aligner/Simulator {#desc}

Vargas aligns short reads to a directed acyclic graph (DAG). Reads are aligned using a SIMD vectorized version of Smith-Waterman, aligning multiple reads at once. The the aligner reports the best positions of reads without traceback to keep memory usage manageable. Additionally, vargas can be used to simulate reads from a DAG, and reads can be stratified by various parameters.

`-h` can be used in any mode to obtain menus.

# Building vargas {#building}

## Cloning the vargas repo {#cloning}

When cloning, use the `--recursive` option to automatically retrieve dependencies.


    git clone --recursive git@github.com:gaddra/vargas.git

## Dependencies {#dependencies}

Vargas relies on htslib to provide core file processing. Once cloned, the library is built with autoconf, version 2.63 or higher.

    cd htslib
    autoconf && ./configure && make

## Building {#build}

vargas is built with cmake.

    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release .. && make
    cd ..
    export PATH=${PWD}/bin:$PATH

# Modes of operation {#modes}

	---------------------- vargas, Jun 24 2016. rgaddip1@jhu.edu ----------------------
	Operating modes 'vargas MODE':
		define      Define a set of graphs for use with sim/align.
		sim         Simulate reads from a set of graphs.
		align       Align reads to a set of graphs.
		test        Run doctests.
		profile     Run profiles.

## define {#define}

	-------------------- vargas define, Jun 24 2016. rgaddip1@jhu.edu --------------------
	-f	--fasta         *<string> Reference filename.
	-v	--var           *<string> VCF/BCF filename.
	-t	--out           *<string> Output filename.
	-g	--region        *<string> Region of graph, format CHR:MIN-MAX.
	-l	--nodelen       <int> Max node length, default 1,000,000
	-i	--ingroup       <int, int...> Percent ingroup subgraphs.
		                Prefix with 'c' to interpret as number of individuals.
	-n	--num           <int> Number of unique graphs to create for all -i. Default 1.

`define` is used to create a graph specification (*.gdf) for `sim` and `align` modes. Various subgraphs deriving off of a base graph are described. With `-i` and `-n` multiple subgraphs based on population filters are made. Graph definition files are marked with `@gdef`. The original FASTA and variant files are needed in the same relative path when using the GDEF file.

For example:

    vargas define -f hs37d5.fa -v chr22.bcf -g 22:1,000,000-2,000,000 -i 15,30 -n 2 -t test.gdf

Will uniquely identify (including implicit graphs) 10 graphs:

2x 15% ingroup graphs

2x 30% ingroup graphs

2x 15% outgroup graphs

2x 30% outgroup graphs

Reference graph

Maximum allele frequency graph

The ingroup determines which samples are included in the graph. If a VCF file has 50 samples, a ingroup of up to 100 individuals can be included in an ingroup; each haplotype is treated independently.

## sim {#sim}


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
	-n reads are produced for each -s, -m, -i, -v combination.


`sim` uses a GDEF file and generates `-n` reads of each combination of `-m`, `-i`, `-v`, and `-b`. `-m` Introduces mutation errors, substituting n bases with an alternate base. Likewise, `-i` will delete a base or insert a random base. With `-a`, `-m` and `-i` are interpreted as rates (0 to 1). `-s` controls which subgraphs are used to generate reads. By default, all graphs in the GDEF file are used, and their complements (outgroup).

For example:

	vargas sim -g test.gdf -t reads -o -n 1000 -m 0,1,2 -v 0,1,2 -j 12

will generate 1000 reads for each combination of `-m`, `-v`, for each graph in `test.gdf` and their complements. This results in 3 -m * 3 -v * 1000 reads * 10 subgraphs = 90,000 reads.

## align {#align}

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
Output format:
	[RXIO], read origin ingroup, read origin graph number, read sequence,
 read origin pos, read sub errors, read indel errors, read variant nodes, read variant bases,
 best score, best pos, best count, sub score, sub pos, sub count, corflag


Reads are aligned to graphs specified in the GDEF file. `-R` Aligns the reads to the reference sequence, `-X` to the maximum allele frequency graph, `-I` to the graph the read was generated from, and `-O` to the complement of the origin graph. A `[RXIO]` tag at the beginning of the alignment record indicates which alignment was done. The output alignment format is:

`[RXIO], read origin ingroup, read origin graph number, read sequence, read origin pos, read sub errors, read indel errors, read variant nodes, read variant bases, best score, best pos, best count, sub score, sub pos, sub count, corflag`

Where `corflag` is 1 if the read origin position had the best score, 2 if it had the second best score, otherwise 0. 

For example:

	vargas align  -g test.gdf -r reads -t reads.aligns -R -I -j 12

Will align all the reads to both the reference graph and ingroup graphs, resulting in 2 * number_of_reads alignments.


