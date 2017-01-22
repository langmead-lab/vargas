[TOC]

_Updated: June 26, 2016_

[![Build Status](https://travis-ci.org/RaviGaddipati/vargas.svg?branch=master)](https://travis-ci.org/RaviGaddipati/vargas)

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
    cmake -DCMAKE_BUILD_TYPE=Release .. && make -j4
    cd ..
    export PATH=${PWD}/bin:$PATH

To  build for Xeon Phi using an intel compiler:

    cmake -DCMAKE_CXX_COMPILER=icpc -DBUILD_PHI=ON -DCMAKE_BUILD_TYPE=Release ..
    make

Note `icpc` requires up to 35GB memory.

# Modes of operation {#modes}

---

`vargas -h`

```
---------------------- vargas, Nov  5 2016. rgaddip1@jhu.edu ----------------------
Operating modes 'vargas MODE':
	define      Define a set of graphs for use with sim and align.
	sim         Simulate reads from a set of graphs.
	align       Align reads to a set of graphs.
	export      Export graph to DOT format.
	convert     Convert a SAM file to a CSV file.
	query       Pull a region from a GDEF/VCF/FASTA file.
	test        Run unit tests.
	profile     Run profiles (debug).
```

## define {#define}

`vargas define -h`

```
-------------------- vargas define, Nov  5 2016. rgaddip1@jhu.edu --------------------
-f	--fasta         *<string> Reference filename.
-v	--vcf           *<string> VCF or BCF file.
-g	--region        *<string> Region of graph, format CHR:MIN-MAX.
-s	--subgraph      <string> Subgraph definition or filename.
-p	--filter        <string> Filename of sample filter.
-x	--invert        Invert sample filter.
-l	--nodelen       <int> Max node length, default 1,000,000
-t	--out           <string> Output file, default stdout.
-d	--dot           <string> Output DOT graph of subgraph hierarchy to file.
Subgraphs are defined using the format "label=N[%t]".
	 where 'N' is the number of samples / percentage of samples selected.
The samples are selected from the parent graph, scoped with ':'.
The BASE graph is implied as the root for all labels. Example:
	a=50;a:b=10%;~a:c=5
'~' indicates the complement graph. 'BASE' refers the the whole graph.
```

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

```
-------------------- vargas sim, Nov  5 2016. rgaddip1@jhu.edu --------------------
-g	--gdef          *<string> Graph definition file. Default stdin.
-s	--sub           <string(;string)*> list of graphs to simulate from. Default all.
-f	--file          -s specifies a file name.
-t	--out           <string> Output file. Default stdout.
-n	--numreads      <int> Number of reads to simulate from each profile, default 1000.
-m	--muterr        <int/float, int/float...> Read mutation error. Default 0.
-i	--indelerr      <int/float, int/float...> Read indel error. Default 0.
-d	--vnodes        <int, int...> Number of variant nodes, default any (*).
-b	--vbases        <int, int...> Number of variant bases, default any (*).
-l	--rlen          <int> Read length, default 50.
-a	--rate          Interpret -m, -i as rates, instead of exact number of errors.
-j	--threads       <int> Number of threads. 0 for maximum hardware concurrency.

-n reads are produced for each -m, -i, -v, -b combination. If set to '*', any value is accepted.

```


`sim` uses a GDEF file and generates `-n` reads of each combination of `-m`, `-i`, `-v`, and `-b`. `-m` Introduces mutation errors, substituting n bases with an alternate base. Likewise, `-i` will delete a base or insert a random base. With `-a`, `-m` and `-i` are interpreted as rates (0 to 1). `-s` controls which subgraphs are used to generate reads. By default, all graphs in the GDEF file are used, and their complements (outgroup).

For example:

	vargas sim -g test.gdf -t reads -o -n 1000 -m 0,1,2 -v 0,1,2 -j 12

will generate 1000 reads for each combination of `-m`, `-v`, for each graph in `test.gdf` and their complements. This results in `3 -m * 3 -v * 1000 reads * 10 subgraphs = 90,000` reads.

## align {#align}

`vargas align -h`

```
------------------- vargas align, Nov  5 2016. rgaddip1@jhu.edu -------------------
-g	--gdef          *<string> Graph definition file.
-r	--reads         *<string> SAM file to align. Default stdin.
-a	--align         *<string:string> Alignment targets, origin graph : target graph.
-f	--file          -a specifies a file name.
-t	--out           *<string> Alignment output file, default stdout.
-l	--rlen          <int> Max read length. Default 50.
-m	--match         <int> Match score, default 2.
-n	--mismatch      <int> Mismatch penalty, default 2.
-o	--gap_open      <int> Gap opening penalty, default 3.
-e	--gap_extend    <int> Gap extend penalty, default 1.
-c	--tolerance     <int> Count an alignment as correct if within this. Default (read_len/2).
-j	--threads       <int> Number of threads. 0 for maximum hardware concurrency.

```


Reads are aligned to graphs specified in the GDEF file. `-R` Aligns the reads to the reference sequence, `-X` to the maximum allele frequency graph, `-I` to the ingroup graph, and `-O` to the outgroup graph. If a set of reads is generated from the ingroup and outgroup using `vargas sim -o`, a subsequent alignment using `vargas align -O` will align all reads to the same outgroup graph. A `[rxio]` tag at the beginning of the alignment record indicates which alignment was done.

For example:

	vargas align  -g test.gdf -r reads -t reads.aligns -R -I -j 12

Will align all the reads to both the reference graph and ingroup graphs, resulting in `2 * number_of_reads` alignments.

## export {#export}

`vargas export -h`

```
-------------------- vargas export, Nov  5 2016. rgaddip1@jhu.edu --------------------
-t	--out           Output filename.
-g	--graph         Subgraph to export, default BASE

vargas export -g "IN" -t out.dot < graph_definition.gdef
```

Exports the graph in DOT format.

## convert {#convert}

`vargas convert -h`

```
-------------------- vargas convert, Nov  5 2016. rgaddip1@jhu.edu --------------------
-s	--sam          <string> SAM input file. Default stdin.
-f	--format       *<string,string...> Specify tags per column. Case sensitive.

Output printed to stdout.
Requred column names: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
Prefix with "RG:" to obtain a value from the associated read group.
```

Convert a SAM file into a CSV file, outputting the specified fields.

## query {#query}

`vargas query -h`

```
-------------------- vargas query, Nov  5 2016. rgaddip1@jhu.edu --------------------
-g	--region        *<string> Region to export, format CHR:MIN-MAX
-f	--fasta         <string> Get a subsequence.
-v	--vcf           <string> VCF or BCF file.
-d	--gdef          <string> Query a graph, export DOT format.
-t	--out           <string> DOT output file for -g.
-s	--subgraph      <string> Subgraph of GDEF to query, default is the whole graph.
```

Extract a subsequence from a graph, FASTA, or VCF file.

## test {#test}

`vargas test`

Execute unit tests.

## profile {#profile}

`vargas profile -h`

```
---------------------- vargas profile, Nov  5 2016. rgaddip1@jhu.edu ----------------------
-f	--fasta         *<string> Reference filename.
-v	--var           *<string> VCF/BCF filename.
-g	--region        *<string> Region of graph, format CHR:MIN-MAX.
-i	--ingroup       <int> Percent of genotypes to include in alignment.
-s	--string        <string,string..> Include reads in alignment. Rest will be random.
```

Generates a summary of performance. For example,

`vargas profile -f ../../hs37d5_22.fa -v ../../chr22.bcf -g 22:22,000,000-32,000,000`

```
Initial Graph Build:
	93.7215 s, Nodes: 909890
Insertion order traversal:
	0.200349 s
Filtering traversal, 100% ingroup:
	0.653423 s
Filtering traversal, 5% ingroup:
	0.593949 s, Nodes: 689949
REF traversal:
	0.482411 s
MAXAF traversal:
	0.48786 s
Filter constructor:
	0.808029 s
REF constructor:
	0.786603 s
MAXAF constructor:
	0.81553 s
16 read alignment:
	Derived Graph (0.904558 s)
	8.70102 s
```