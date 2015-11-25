# Variant Genome Aligner/Simulator
Vargas uses gssw to align reads to a partial order graph without traceback. Alignment information is given as the ending position of the best score, as well as a suboptimal alignment score. There are four main modes of operation.

###Basic Structure
All objects are in the vargas namespace with the exception of `utils.cpp`.

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


###Modes of operation
 ```
  ---------------------- vargas, September 2015. rgaddip1@jhu.edu ----------------------
  Operating modes 'vargas MODE':
  build     Generate graph build file from reference and VCF files.
  sim       Simulate reads from a graph.
  align     Align reads to a graph.
  export    Export graph in DOT format.
 ```
 
 ```
 ------------------- vargas build, September 2015. rgaddip1@jhu.edu -------------------
-v     --vcf           (required) VCF file, uncompressed.
-r     --ref           (required) reference single record FASTA
-l     --maxlen        Maximum node length
-R     --region        [min:max] Ref region, inclusive. Default is entire graph.
-g     --ingroup       Percent of individuals to build graph from, default all.
-c     --complement    Generate a complement of a graph (outgroup becomes ingroup).

Buildfile is printed on stdout.
 ```
 
 ```
------------------- vargas sim, September 2015. rgaddip1@jhu.edu -------------------
-b     --buildfile     quick rebuild file, required if -v, -r are not defined.
-n     --numreads      Number of reads to simulate
-m     --muterr        Simulated read mutation error rate
-i     --indelerr      Simulated read Indel error rate
-r     --rand          Use a random mutation error rate, up to the value specified by -m.
-l     --readlen       Nominal read length
-e     --regex         Match regex expressions. Produces -n of each, discard others.
                       List of expressions is space delimited -e "exp1 exp2"
  
Reads are printed on stdout.
Read Format:
READ#READ_END_POSITION,INDIVIDUAL,NUM_SUB_ERR,NUM_VAR_NODE,NUM_VAR_BASES
 ```
 
 ```
 ------------------- vargas align, September 2015. rgaddip1@jhu.edu -------------------
-b     --buildfile     quick rebuild file, required if -v, -r are not defined.
-m     --match         Match score, default 2
-n     --mismatch      Mismatch score, default 2
-o     --gap_open      Gap opening score, default 3
-e     --gap_extend    Gap extend score, default 1
-r     --reads         Reads to align, one per line. Symbols after '#' are ignored.

Alignments output to stdout.
Output format:
READ,OPTIMAL_SCORE,OPTIMAL_ALIGNMENT_END,NUM_OPTIMAL_ALIGNMENTS,SUBOPTIMAL_SCORE,
SUBOPTIMAL_ALIGNMENT_END,NUM_SUBOPTIMAL_ALIGNMENTS,ALIGNMENT_MATCH

ALIGNMENT_MATCH: 0- optimal match, 1- suboptimal match, 2- no match
```

```
------------------- vargas export, September 2015. rgaddip1@jhu.edu -------------------
-b     --buildfile    Graph to export to DOT.

DOT file printed to stdout.
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
 
 A graph is built from a buildfile. A graph buildfile can be made using `vmatch build`.
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
