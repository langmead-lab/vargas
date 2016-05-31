# Variant Genome BatchAligner/Simulator
Vargas uses gssw to align reads to a partial order Graph without traceback. Alignment information is given as the ending position of the best score, as well as a suboptimal alignment score. There are four main modes of operation.

###Basic Structure
All objects are in the vargas namespace with the exception of `utils` functions.

`vcfstream.h` : vcfstream is a wrapper for a VCF file. It reads, parses, and filters the VCF file. Variants are provided
in a VariantRecord struct. Individual variants not included in the ingroup are removed from the records. Associated allele
frequency values are also provided.

`readsim.h` : Inherits from `readsource.h`. Reads are simulated from a given `gssw_graph` or `Graph` object. Reads 
can be provided indefinetly or may be written to a set of files as determined by regular expressions. Reads are provided
in a Read struct.

`_read_file.h` : Inherits from `readsource.h`. Wrapper for a reads file (such as that produced from `readsim`). Provides reads in a Read struct.

`Graph.h` : Contains facilities for building and aligning to a Graph. A buildfile is first exported from a `vcfstream` and
a `readsource`. Graphs are built in memory from the buildfile. Graph parameters are defined in a GraphParams struct.
Alignments to the Graph are done with `Graph.align(Read)` and returns with an Alignment struct.

`xcoder.h` : Uses masked VByte compression and base64 encoding to reduce the size of the list of individuals for each variant in the buildfile. Lists of individuals are stored in nodes using VByte compression.


###Modes of operation
 ```
---------------------- vargas, Apr 23 2016. rgaddip1@jhu.edu ----------------------
Operating modes 'vargas MODE':
	build     Generate Graph build file from reference FASTA and VCF files.
	sim       Simulate reads from a Graph.
	align     Align reads to a Graph.
	stat      Count nodes and edges of a given Graph.
	export    Export Graph in DOT format.

 ```
 
 ```
------------------- vargas build, Apr 23 2016. rgaddip1@jhu.edu -------------------
-v	--vcf           <string> VCF file, uncompressed.
-r	--ref           <string> reference, single record FASTA
-l	--maxlen        <int> Maximum node length, default 50000
-R	--region        <<int>:<int>> Ref region, inclusive. Default is the entire Graph.
-m	--maxref        Generate a linear Graph using maximum allele frequency nodes.
-e	--exref         Exclude the list of individuals from the reference alleles.
-s	--set           <<int>,<int>,..,<int>> Generate a buildfile for a list of ingroup percents.
-c	--complement    <string> Generate a complement of all graphs in -s, or of provided Graph.

--maxref is applied after ingroup filter.
Buildfile is output to [s][In Out].build

 ```
 
 ```
-------------------- vargas sim, Apr 23 2016. rgaddip1@jhu.edu --------------------
-b	--buildfile     <string> Graph build file, generate with 'vargas build'
-n	--numreads      <int> Number of reads to simulate, default 10000
-m	--muterr        <float> Read mutation error rate, default 0
-i	--indelerr      <float> Read Indel error rate, default 0
-l	--readlen       <int> Read length, default 100
-e	--profile       <p1 p2 .. p3> Space delimited read profiles. Produces -n of each
-p	--prefix        Prefix to use for read files, default 'sim'
-r	--randwalk      Random walk, read may change individuals at branches
-a	--ambiguity     Max number of ambiguous bases to allow in reads, default 0

Outputs to '[prefix][n].reads' where [n] is the profile number.
Read Profile format (use '*' for any): 
	sub_err,indel_err,var_nodes,var_bases
	Example: Any read with 1 substitution error and 1 variant node.
		vargas sim -b BUILD -e "1,*,1,*"
Read Format:
	READ#READ_END_POSITION,INDIVIDUAL,NUM_SUB_ERR,NUM_INDEL_ERR,NUM_VAR_NODE,NUM_VAR_BASES

 ```
 
 ```
------------------- vargas align, Apr 23 2016. rgaddip1@jhu.edu -------------------
-b	--buildfile     <string> Graph build file
-m	--match         <int> Match score, default 2
-n	--mismatch      <int> Mismatch score, default 2
-o	--gap_open      <int> Gap opening penalty, default 3
-e	--gap_extend    <int> Gap extend penalty, default 1
-r	--reads         <string> Reads to align. Use stdin if not defined.
-f	--outfile       <string> Alignment output file. If not defined, use stdout.
Lines beginning with '#' are ignored.
Output format:
	READ,OPTIMAL_SCORE,OPTIMAL_ALIGNMENT_END,NUM_OPTIMAL_ALIGNMENTS,SUBOPTIMAL_SCORE,SUBOPTIMAL_ALIGNMENT_END,NUM_SUBOPTIMAL_ALIGNMENTS,ALIGNMENT_MATCH

ALIGNMENT_MATCH:
	0- optimal match, 1- suboptimal match, 2- no match

```

------------------ vargas export, Apr 23 2016. rgaddip1@jhu.edu -------------------
-b	--buildfile    <string> Graph to export to DOT.
-c	--context      <string> Export the local context Graph of these alignments.

DOT file printed to stdout.

```

### Cloning the vargas repo

When cloning, use the `--recursive` option to automatically retrieve dependencies.

```
git clone --recursive git@github.com:gaddra/vargas.git
```

### Building vargas

vargas is built with cmake.

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
cd ..
export PATH=${PWD}/bin:$PATH
```

### Generating a buildfile
 
 A Graph is built from a buildfile. A Graph buildfile can be made using `vargas build`.
  `-r` and `-v` specify the reference FASTA and VCF respectively. A region can be specified with `-R`.
  The percentage of individuals to include in the Graph is specified with `-s` and a comma separated list of values.
  `-c` will generate complements of each Graph.
   For example to generate a buildfile and its complement from `REF` and `VAR` from position `a` to position `b` with 50% of individuals:
 
 `vargas build -r REF.fa -v VAR.vcf -R a:b -s 50 -c > GRAPH.build`
 
### Simulating reads
 
 Reads can be simulated with `vargas sim`. Mutation error and Indel errors can be introduced.
 Reads can take random paths through the Graph or from a random individual. To generate 1000 reads of length 100 and a 1% error rate from a Graph:
 
 `vargas sim -b GRAPH -n 1000 -l 100 -m 0.01 -i 0.01`
 
 To target certain kinds of reads, a list of space delimited profiles can be provided with `-e`. `-n` reads will be generated for each profile. Any reads that do not match the profile will be discarded.
 
 Note: Simulating reads may consume a large amount of memory as it loads the (compressed) list of individuals into each variant node. `hs37d5 chr22 -R 22000000:52000000` was at ~10GB usage.

### Aligning to the Graph
 
To process a reads file and write the output to a file:
 
 `vargas align -b GRAPH -d READS > ALIGNMENTS`
 
If the read contains `#` (as in the simulated reads), everything after it is stripped for alignment but preserved in the output.

If a tie for a best score is found, the alignment closer to the read origin is preserved.

### File formats

#### Aligns

The format of an alignment result is

```
READ,OPTIMAL_SCORE,OPTIMAL_ALIGNMENT_END,NUM_OPTIMAL_ALIGNMENTS,SUBOPTIMAL_SCORE,SUBOPTIMAL_ALIGNMENT_END,NUM_SUBOPTIMAL_ALIGNMENTS,ALIGNMENT_MATCH
```

where `READ` is in the format described below. `AlIGNMENT_MATCH` is a flag that determines the best alignment relative to the simulated position. `0` indicates a match with the simulated origin, `1` is a second-best alignment match with the simulated origin, `2` is other. 

#### Reads

The beginning of a simulated reads file will contain a comment indicating the profile used to generate reads, where `-1` means anything. The second line includes the random seed used to generate the reads, as well as parameters passed to the simulator.

```
#indiv=-1,sub_err=0,var_nodes=2,var_bases=-1
#Seed: 1451149022, Mut err: 0.02, Indel Err: 0, read Len: 64, max reads: 10000, random walk? 0
```

The format of a read is below, where `-1` in the `INDIVIDUAL` column indicates the entire read was from the reference.

```
READ_SEQUENCE#READ_END_POSITION,INDIVIDUAL,NUM_SUB_ERR,NUM_INDEL_ERR,NUM_VAR_NODE,NUM_VAR_BASES
```

#### Buildfile

The first line of the buildfile lists the individuals included in the Graph, where each number represents the column of the individual. Each haplotype is its own column.

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

Currently only `CN` tags are supported. If a VCF record specifies a reference of `[seq]` and a variant of `<CNn>` where `n` is an integer, the variant node will copy `[seq]` `n` times. VCF records containing other tags (such as mobile elements) are ignored.
