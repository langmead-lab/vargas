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
    cmake -DCMAKE_BUILD_TYPE=Release .. && make
    cd ..
    export PATH=${PWD}/bin:$PATH

# Modes of operation {#modes}

---

`vargas -h`

```
---------------------- Vargas, Aug 31 2016. rgaddip1@jhu.edu ----------------------
Operating modes 'Vargas MODE':
	define      Define a set of graphs for use with sim/align.
	sim         Simulate reads from a set of graphs.
	align       Align reads to a set of graphs.
	split       Split a SAM file into multiple files.
	merge       Merge SAM files.
	convert     Convert a SAM file to a CSV file.
	test        Run doctests.
	profile     Run profiles.
```

## define {#define}

`vargas define -h`

```
-------------------- Vargas define, Aug 31 2016. rgaddip1@jhu.edu --------------------
-f	--fasta         *<string> Reference filename.
-v	--vcf           <string> VCF or BCF file.
-k	--ksnp          <string> KSNP File. Only one of -v or -k should be defined.
-g	--region        *<string> Region of graph, format CHR:MIN-MAX.
-l	--nodelen       <int> Max node length, default 1,000,000
-s	--subgraph      *<string> Subgraph definition file. Default stdin.
-t	--out           <string> Output file, default stdout.
-d	--dot           <string> Output DOT graph to file.

Subgraphs are defined using the format "label=N[%]".
'N' is the number of samples / percentage of samples selected.
The samples are selected from the parent graph, scoped with ':'.
The base graph is implied as the root for all labels. Example:
a=50;a:b=10;~a:c=5
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
-------------------- Vargas sim, Aug 31 2016. rgaddip1@jhu.edu --------------------
-g	--gdef          *<string> Graph definition file. Default stdin.
-s	--sub           *<string;string; ...> list of graphs to simulate from.
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
------------------- Vargas align, Aug 31 2016. rgaddip1@jhu.edu -------------------
-g	--gdef          *<string> Graph definition file.
-r	--reads         *<string, string...> Read files to align. Default stdin.
-a	--align         *<string> Alignment targets.
-f	--file          -a specifies a file name.
-t	--out           *<string> Alignment output file, default stdout.
-l	--rlen          <int> Max read length. Default 50.
-m	--match         <int> Match score, default 2.
-n	--mismatch      <int> Mismatch penalty, default 2.
-o	--gap_open      <int> Gap opening penalty, default 3.
-e	--gap_extend    <int> Gap extend penalty, default 1.
-j	--threads       <int> Number of threads. 0 for maximum hardware concurrency.
            	            Optimal reads per subgraph: n * j * 16
```


Reads are aligned to graphs specified in the GDEF file. `-R` Aligns the reads to the reference sequence, `-X` to the maximum allele frequency graph, `-I` to the ingroup graph, and `-O` to the outgroup graph. If a set of reads is generated from the ingroup and outgroup using `vargas sim -o`, a subsequent alignment using `vargas align -O` will align all reads to the same outgroup graph. A `[rxio]` tag at the beginning of the alignment record indicates which alignment was done.

For example:

	vargas align  -g test.gdf -r reads -t reads.aligns -R -I -j 12

Will align all the reads to both the reference graph and ingroup graphs, resulting in `2 * number_of_reads` alignments.

## convert {#convert}

`vargas convert -h`

```
-------------------- Vargas convert, Aug 31 2016. rgaddip1@jhu.edu --------------------
-s	--sam          <string> SAM input file. Default stdin.
-f	--format       *<string,string...> Specify tags per column. Case sensitive.

Output printed to stdout.
Requred column names: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
Prefix with "RG:" to obtain a value from the associated read group.
```

Convert a SAM file into a CSV file, outputting the specified fields.
