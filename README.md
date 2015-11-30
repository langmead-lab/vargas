# Variant Genome Aligner/Simulator
Vargas uses gssw to align reads to a partial order graph without traceback. Alignment information is given as the ending position of the best score, as well as a suboptimal alignment score. There are four main modes of operation.

###Basic Structure
All objects are in the vargas namespace with the exception of `utils` functions.

`vcfstream.h` : vcfstream is a wrapper for a VCF file. It reads, parses, and filters the VCF file. Variants are provided
in a VariantRecord struct. Individual variants not included in the ingroup are removed from the records. Associated allele
frequency values are also provided.

`readsource.h` : Abstract class defining functions for classes that provide reads for alignment.

`readsim.h` : Inherits from `readsource.h`. Reads are simulated from a given `gssw_graph` or `Graph` object. Reads 
can be provided indefinetly or may be written to a set of files as determined by regular expressions. Reads are provided
in a Read struct.

`readfile.h` : Inherits from `readsource.h`. Wrapper for a reads file (such as that produced from `readsim`). Provides reads in a Read struct.

`graph.h` : Contains facilities for building and aligning to a graph. A buildfile is first exported from a `vcfstream` and
a `readsource`. Graphs are built in memory from the buildfile. Graph parameters are defined in a GraphParams struct.
Alignments to the graph are done with `Graph.align(Read)` and returns with an Alignment struct.

`xcoder.h` : Uses masked VByte compression and base64 encoding to reduce the size of the list of individuals for each variant.


###Modes of operation
 ```
---------------------- vargas, Nov 25 2015. rgaddip1@jhu.edu ----------------------
Operating modes 'vargas MODE':
	build     Generate graph build file from reference and VCF files.
	sim       Simulate reads from a graph.
	align     Align reads to a graph.
	export    Export graph in DOT format.
 ```
 
 ```
------------------- vargas build, Nov 25 2015. rgaddip1@jhu.edu -------------------
-v	--vcf           (required) VCF file, uncompressed.
-r	--ref           (required) reference, single record FASTA
-l	--maxlen        Maximum node length
-R	--region        <min:max> Ref region, inclusive. Default is the entire graph.
-m	--maxref        Generate a graph using alleles in the ingroup w/ the highest frequency.
-s	--set           <#,#,..,#> Generate a buildfile for a list of ingroup %'s and their complements.
-c	--complement    <graph.build> Generate a complement of all graphs in -s

Buildfile is output to [s][In/Out].build

 ```
 
 ```
------------------- vargas sim, Nov 25 2015. rgaddip1@jhu.edu -------------------
-b	--buildfile     quick rebuild file, generate with vargas build
-n	--numreads      Number of reads to simulate
-m	--muterr        Simulated read mutation error rate
-i	--indelerr      Simulated read Indel error rate
-l	--readlen       Nominal read length
-e	--regex         <r1 r2 .. r3> Match regex expressions, space delimited. Produces -n of each.
-p	--prefix        Prefix to use for read files generated with -e
-r	--randwalk      Random walk, read may change individuals at branches.

NOTE: End of line anchor may not work in regex depending on C++ version. 
Reads are printed on stdout.
Read Format:
READ#READ_END_POSITION,INDIVIDUAL,NUM_SUB_ERR,NUM_VAR_NODE,NUM_VAR_BASES

 ```
 
 ```
------------------- vargas align, Nov 25 2015. rgaddip1@jhu.edu -------------------
-b	--buildfile     Quick rebuild file.
-m	--match         Match score, default  2
-n	--mismatch      Mismatch score, default 2
-o	--gap_open      Gap opening penalty, default 3
-e	--gap_extend    Gap extend penalty, default 1
-r	--reads         Reads to align.

Alignments output to stdout. Reads read from stdin or -r, 1 per line.
Lines beginning with '#' are ignored.
Output format:
READ,OPTIMAL_SCORE,OPTIMAL_ALIGNMENT_END,NUM_OPTIMAL_ALIGNMENTS,SUBOPTIMAL_SCORE,
SUBOPTIMAL_ALIGNMENT_END,NUM_SUBOPTIMAL_ALIGNMENTS,ALIGNMENT_MATCH

ALIGNMENT_MATCH: 0- optimal match, 1- suboptimal match, 2- no matc
```

```
------------------- vargas export, Nov 25 2015. rgaddip1@jhu.edu -------------------
-b	--buildfile    (required) Graph to export to DOT.

DOT file printed to stdout.

```

### Cloning the vargas repo

When cloning, use the `--recursive` option to automatically retrieve dependencies.  For example:

```
git clone --recursive git@github.com:gaddra/vargas.git
```

### Building vargas

vargas is built with cmake.

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
export PATH=${PWD}/bin:$PATH
```

### Generating a buildfile
 
 A graph is built from a buildfile. A graph buildfile can be made using `vargas build`.
  `-r` and `-v` specify the reference FASTA and VCF respectively. A region can be specified with `-R`.
  The percentage of individuals to include in the graph is specified with `-s` and a comma seperated list of values.
  `-c` will generate complements of each graph.
   For example to generate a buildfile and its complement from `REF` and `VAR` from position `a` to position `b` with 50% of individuals:
 
 `vargas build -r REF.fa -v VAR.vcf -R a:b -s 50 -c > GRAPH.build`
 
### Simulating reads
 
 Reads can be simulated with `vargas sim`. Mutation error and Indel errors can be introduced (default 0%).
 Reads can take random paths through the graph or from a random individual. To generate 1000 reads of length 100 and a 1% error rate from a graph:
 
 `vargas sim -b GRAPH -n 1000 -l 100 -m 0.01 -i 0.01 > SIMREADS`
 
 To target certain kinds of reads, a list of space delimited regular expressions can be provided with `-e`. `-n` reads will be genereated for each regular expression. Any reads that do not match the expression will be discarded.

### Aligning to the graph
 
To process a reads file and write the output to a file:
 
 `vargas align -b GRAPH -d READS > ALIGNMENTS`
 
If the read contains `#` (as in the simulated reads), everything after it is stripped for alignment but preserved in the output.

### File formats

#### Buildfile

The first line of the buildfile lists the individuals included in the graph, where each number represents the column of the individual. Each haplotype is its own column.

Following lines beginning with `##` are comment lines and include information on the configuration used to make the buildfile.

Each line consisting of 3 comma separated values represents a new node, with the format `End of node position, node number, node sequence`.

Each line consisting of 2 comma separated values represents an edge relative to the nodes above it. E.g. `-2,-1` means there is an edge from the second to last node listed thus far in the buildfile to the most recently listed node.

Each line beginning with a `:` is a base64 encoded, masked Vbyte compressed, list of numbers representing the individuals that possess the allele directly before it.

Example:
The reference `A` allele is in individual `9` and the `GA` allele belongs to individual `10`.
```
  G
 / \
A   GA
 \ /
  T
```

```
#9,10,
##<Comments>
1,0,G
2,1,A
:9       <-- This would be encoded
2,2,GA
:10      <-- This would be encoded
-3,-2
-3,-1
3,3,T
-3,-1
-2,-1
```
