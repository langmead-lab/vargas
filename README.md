# VMatch
 VMatch uses gssw to align reads to a partial order graph without traceback. Alignment information is given as the position of the best score, as well as a suboptimal alignment score.
 ```
-v     --vcf           VCF file, uncompressed
-r     --ref           reference FASTA
-b     --buildfile     quick rebuild file, required if -v, -r are not defined. Takes priority.
-B     --NVbuildfile   quick rebuild file for no variant graph, use with -D.
-m     --match         Match score, default 2
-n     --mismatch      Mismatch score, default 2
-o     --gap_open      Gap opening score, default 3
-e     --gap_extend    Gap extend score, default 1
-t     --outfile       Output file for quick rebuild of graph
-d     --reads         Reads to align, one per line. set equal to -T to align after sim
-s     --string        Align a single string to stdout, overrides reads file
-a     --aligns        Outputfile for alignments
-R     --region        Ref region, inclusive. min:max
-p     --noprint       Disable stdout printing
-x     --novar         Generate and align to a no-variant graph
-D     --dual          Align to both variant and non-variant graphs.
-i     --simreads      Simulate reads and write to the file specified by -T
-T     --readout       File to output reads to
-N     --numreads      Number of reads to simulate, default
-M     --muterr        read mutation error rate, default
-I     --indelerr      read Indel error rate, default
-L     --readlen       nominal read length, default
 ```

### Building VMatch

VMatch is built with cmake.

```
mkdir build && cd build
cmake ..
make
export PATH=path_to_VMatch/build:$PATH
```

For an optimizied build a higher version of GCC is needed, such as gcc-4.9.3.

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
export PATH=path_to_VMatch/build:$PATH
```

Before usage of vmatch,

`export LD_LIBRARY_PATH=PATH_TO_GCC-4.9.3/lib64:$LD_LIBRARY_PATH`

Using chr22 in the region of 21,000,000:51,000,000 with an optimized build, it takes ~10s to align to a variant and a no variant graph using the `-D` option in a 3GB memory footprint.

### Generating a graph
 
 A graph is built on every launch of the program. A graph can be built using the `-v` and `-r` arguments, specifying the VCF and reference FASTA resepectivly. For large graphs the the graph can be output with the `-t` option. This allows for a quick-rebuild of the graph in future alignments using the `-b` argument. A region can be specified with `-R`. A graph with no variants is built using the `-x` option. For example to generate a graph from `REF` and `VAR` from position `a` to position `b` without variants:
 
 `vmatch -r REF.fa -v VAR.vcf -t GRAPH -R a:b -x`
 
### Simulating reads
 
 Reads can be simulated with `-i`. The reads are written to the file specified by `-T`. Mutation error and Indel errors can be introduced (default 1%). Reads take random paths through the graph. To generate 1000 reads of length 100 and a 1% error rate from a previously built graph:
 
 `vmatch -b GRAPH -i -N 1000 -L 100 -M 0.01 -I 0.01 -T SIMREADS`
 
 Each read is appended with `#Node ID, node max position, end of read position in node`. 
 
### Aligning to the graph
 
 Alignment can be done from a reads file or specified with `-s`. Aligning a single read will print to stdout. To process a reads file and write the output to a file:
 
 `vmatch -b GRAPH -d READS -a ALIGNMENTS`
 
 If the read contains `#`, everything after it is stripped for alignment but preserved in the output. If the `-d` argument matches the `-T` argument, it is equivilant to generating random reads and aligning them. To build a graph, generate 1000 reads and align them:
 
 `vmatch -r REF.fa -v VAR.vcf -t GRAPH -R a:b -i -N 1000 -T SIMREADS -d SIMREADS -a ALIGNMENTS`
 
 Reads can also be aligned to two graphs (intended for a variant and a no varant graph) using the `-D` option and specifying two build files. To generate 10000 reads and align it to both graphs:
 
 `vmatch -D -b GRAPH -B GRAPH_NV -d READS -a ALIGNS -i -N 10000 -T READS`
 
 An example output with the corresponding format is below.
 
 Sim read format: `READ`#`NODE_ID`, `NODE_LEN`, `NODE_MAX_POSITION`, `READ_END_POSITION`
 
 Alignment format: #`READ`; `1_NODE_ID`, `NODE_LEN`, `1_NODE_MAX`, `1_SCORE`, `1_END_POS`; `2_NODE_ID`, `2_NODE_LEN`, `2_NODE_MAX`, `2_SCORE`, `2_END_POS`
 
```
CAGACAAATCTGGGTTCAAGTCCTCACTTT#54,13,217,10;54,13,217,60,9;0,8,8,10,4
GCTGCTCTCTTCTTGTCAGATCGTATTCTC#197,14,915,5;197,14,915,56,4;0,8,8,6,4
GAATCTTTCCAGAACCTGCTCTTTCCTCA#178,30,857,15;178,30,857,55,14;0,8,8,6,4
GTGAAGCTGAGGGAATAGTGCCTGGCATAG#93,23,396,1;93,23,396,60,0;0,8,8,8,5
GTTACTGTTATTTACTATGAATCCTCACCT#104,57,465,19;104,57,465,60,18;0,8,8,6,4
```
