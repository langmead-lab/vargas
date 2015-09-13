# VMatch
 VMatch uses gssw to align reads to a partial order graph without traceback. Alignment information is given as the ending position of the best score, as well as a suboptimal alignment score. There are four main modes of operation.
 ```
  ---------------------- VMatch, September 2015. rgaddip1@jhu.edu ----------------------
  Operating modes 'vmatch MODE':
  build     Generate graph build file from reference and VCF files.
  sim       Simulate reads from a graph.
  align     Align reads to a graph.
  export    Export graph in DOT format.
 ```
 
 ```
 ------------------- VMatch build, September 2015. rgaddip1@jhu.edu -------------------
-v     --vcf           (required) VCF file, uncompressed.
-r     --ref           (required) reference single record FASTA
-l     --maxlen        Maximum node length
-R     --region        [min:max] Ref region, inclusive. Default is entire graph.
-g     --ingroup       Percent of individuals to build graph from, default all.
-c     --complement    Generate a complement of a graph (outgroup becomes ingroup).

Buildfile is printed on stdout.
 ```
 
 ```
------------------- VMatch sim, September 2015. rgaddip1@jhu.edu -------------------
-b     --buildfile     quick rebuild file, required if -v, -r are not defined.
-n     --numreads      Number of reads to simulate
-m     --muterr        Simulated read mutation error rate
-i     --indelerr      Simulated read Indel error rate
-l     --readlen       Nominal read length
  
Reads are printed on stdout.
Read Format:
READ#READ_END_POSITION,INDIVIDUAL,NUM_SUB_ERR,NUM_VAR_NODE,NUM_VAR_BASES
 ```
 
 ```
 ------------------- VMatch align, September 2015. rgaddip1@jhu.edu -------------------
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
------------------- VMatch export, September 2015. rgaddip1@jhu.edu -------------------
-b     --buildfile    Graph to export to DOT.

DOT file printed to stdout.
```
### Building VMatch

VMatch is built with cmake.

```
mkdir build && cd build
cmake ..
make
export PATH=path_to_VMatch/build:$PATH
```

### Generating a graph
 
 A graph is built on every launch of the program. A graph can be built using the `-v` and `-r` arguments, specifying the VCF and reference FASTA resepectivly and is used to generate a build file. A region can be specified with `-R`. The percentage of individuals ot include in the graph is specified with `-g` For example to generate a graph from `REF` and `VAR` from position `a` to position `b` with 50% of individuals:
 
 `vmatch build -r REF.fa -v VAR.vcf -R a:b -g 50 > GRAPH.build`
 
### Simulating reads
 
 Reads can be simulated with `vmatch sim`. Mutation error and Indel errors can be introduced (default 1%). Reads take random paths through the graph. To generate 1000 reads of length 100 and a 1% error rate from a graph:
 
 `vmatch sim -b GRAPH -n 1000 -l 100 -m 0.01 -i 0.01 > SIMREADS`

### Aligning to the graph
 
To process a reads file and write the output to a file:
 
 `vmatch align -b GRAPH -d READS > ALIGNMENTS`
 
If the read contains `#` (as in the simulated reads), everything after it is stripped for alignment but preserved in the output.
