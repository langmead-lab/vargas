# VMatch
 VMatch uses gssw to align reads to a partial order graph without traceback. Alignment information is given as the position of the best score, as well as a suboptimal alignment score.
 ```
-v     --vcf           VCF file, uncompressed
-r     --ref           reference FASTA
-b     --buildfile     quick rebuild file, required if -v, -r are not defined. Takes priority.
-m     --match         Match score, default  " << matc
-n     --mismatch      Mismatch score, default " << mismatc
-o     --gap_open      Gap opening score, default " << int32_t(gap_open
-e     --gap_extend    Gap extend score, default " << int32_t(gap_extension
-t     --outfile       Output file for quick rebuild of graph
-d     --reads         Reads to align, one per line. set equal to -T to align after sim
-s     --string        Align a single string to stdout, overrides reads file
-a     --aligns        Outputfile for alignments
-R     --region        Ref region, inclusive. min:max
-p     --noprint       Disable stdout printing
-x     --novar         Generate and align to a no-variant graph
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
 
 An example output with the corresponding format is below.
 
 `Simulated read` # `Node ID`,`node max position`,`end of read position in node`;`node ID`,`max position in node`,`best score`,`best score position in node`,`suboptimal score`,`suboptimal end in node`
 
 ```
>GGAGG#110,488,2;12,33,8,12,0,-1
>AGATC#45,186,0;44,178,10,0,0,-1
>ATCGC#122,552,4;93,396,8,0,0,-1
>ATCTC#66,286,4;60,271,10,35,0,-1
```
