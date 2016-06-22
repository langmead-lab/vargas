# Variant Genome Aligner/Simulator

Vargas aligns short reads to a directed acyclic graph (DAG). Reads are aligned using a SIMD vectorized version of Smith-Waterman, aligning multiple reads at once. The the aligner reports the best positions of reads without traceback to keep memory usage manageable. Additionally, vargas can be used to simulate reads from a DAG, and reads can be stratified by various parameters.

`-h` can be used in any mode to obtain menus.

## Modes of operation

```
---------------------- vargas, Jun 21 2016. rgaddip1@jhu.edu ----------------------
Operating modes 'vargas MODE':
	define      Define a set of graphs for use with sim/align.
	sim         Simulate reads from a graph.
	align       Align reads to a graph.
	test        Run doctests.
	profile     Run profiles.
```

### define

```
-------------------- vargas define, Jun 21 2016. rgaddip1@jhu.edu --------------------
-f	--fasta         <string> Reference filename.
-v	--var           <string> VCF/BCF filename.
-t	--out           <string> Output filename.
-g	--region        <string> Region of graph, format CHR:MIN-MAX.
-l	--nodelen       <int> Max node length, default 1,000,000
-b	--base          <int> Base graph ingroup. -i graphs derived from this. Default 100%.
-i	--ingroup       <int, int...> Percent ingroup subgraphs.
-n	--num           <int> Number of unique graphs to create for all -i, -x. Default 1.
```

`define` is used to create a graph specification (*.gdf) for `sim` and `align` modes. Various subgraphs deriving off of a base graph are described. With `-i` and `-n` multiple subgraphs based on population filters are made. Graph definition files are marked with `@gdef`. The original FASTA and variant files are needed in the same relative path when using the GDEF file.

For example:

```
vargas define -f hs37d5.fa -v chr22.bcf -g 22:1,000,000-2,000,000 -i 15,30 -n 2
```

Will uniquely identify (including implicit graphs) 10 graphs:

2x 15% ingroup graphs

2x 30% ingroup graphs

2x 15% outgroup graphs

2x 30% outgroup graphs

Reference graph

Maximum allele frequency graph

### sim

```
-------------------- vargas sim, Jun 21 2016. rgaddip1@jhu.edu --------------------
-g	--gdef          <string> Graph definition file.
-t	--outfile       <string> Alignment output file.
-s	--source        <string, string...> Simulate from specified subgraphs, default all.
-n	--numreads      <int> Number of reads to simulate from each subgraph.
-m	--muterr        <float, float...> Read mutation error. Default 0.
-i	--indelerr      <float, float...> Read indel error. Default 0.
-v	--vnodes        <int, int...> Number of variant nodes, default any (-1).
-b	--vbases        <int, int...> Number of variant bases, default any (-1).
-l	--rlen          <int> Read length, default 50.
-a	--rate          Interpret -m, -i as rates, instead of exact number of errors.
-e	--exclude       Exclude comments, include only FASTA lines
-j	--threads       <int> Number of threads. 0 for maximum hardware concurrency.
```

`sim` uses a GDEF file and generates `-n` reads of each combination of `-m`, `-i`, `-v`, and `-b`. `-m` Introduces mutation errors, substituting n bases with an alternate base. Likewise, `-i` will delete a base or insert a random base. With `-a`, `-m` and `-i` are interpreted as rates (0 to 1). `-s` controls which subgraphs are used to generate reads. By default, all graphs in the GDEF file are used, and their complements (outgroup).

### align

```
------------------- vargas align, Jun 21 2016. rgaddip1@jhu.edu -------------------
-g	--gdef          <string> Graph definition file.
-r	--reads         <string, string...> Read files to align
-t	--outfile       <string> Alignment output file.
-l	--rlen          <int> Max read length. Default 50.
-R	                Align to reference graph.
-X	                Align to maximum allele frequency graph.
-I	                Align to ingroup graph.
-O	                Align to outgroup graph.
-m	--match         <int> Match score, default 2
-n	--mismatch      <int> Mismatch penalty, default 2
-o	--gap_open      <int> Gap opening penalty, default 3
-e	--gap_extend    <int> Gap extend penalty, default 1
-j	--threads       <int> Number of threads. 0 for maximum hardware concurrency.
```

Reads are aligned to graphs specified in the GDEF file. `-R` Aligns the reads to the reference sequence, `-X` to the maximum allele frequency graph, `-I` to the graph the read was generated from, and `-O` to the complement of the origin graph. A `[RXIO]` tag at the beginning of the alignment record indicates which alignment was done. The output alignment format is:

`[RXIO], read origin ingroup, read origin graph number, read sequence, read origin pos, read sub errors, read indel errors, read variant nodes, read variant bases, best score, best pos, best count, sub score, sub pos, sub count, corflag`

Where `corflag` is 1 if the read origin position had the best score, 2 if it had the second best score, otherwise 0. 

### Cloning the vargas repo

When cloning, use the `--recursive` option to automatically retrieve dependencies.

```
git clone --recursive git@github.com:gaddra/vargas.git
```

### Building vargas

vargas is built with cmake. `htslib` needs to be configured and built first.

```
cd htslib
autoconf && ./configure && make
cd ..
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
cd ..
export PATH=${PWD}/bin:$PATH
```