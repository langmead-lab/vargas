# Vargas
_Updated: May 11, 2017_

[![Build Status](https://travis-ci.org/RaviGaddipati/vargas.svg?branch=master)](https://travis-ci.org/RaviGaddipati/vargas)

Vargas aligns short reads to a directed acyclic graph (DAG). Reads are aligned using a SIMD vectorized version of Smith-Waterman, aligning multiple reads at once. The the aligner reports the maximal scoring positions. Vargas can also be used to simulate reads from a DAG.


# Building

When cloning, use the `--recursive` option to automatically retrieve dependencies.

    git clone --recursive git@github.com:RaviGaddipati/vargas.git


Vargas relies on htslib to provide core file processing. Once cloned, the htslib is built with autoconf (version 2.63+).

    cd htslib
    autoconf && ./configure && make -j4

Vargas is built with cmake, and targets SSE4.1 for SIMD support. AVX512 is targeted on Phi platforms.

    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release .. && make -j4
    cd ..
    export PATH=${PWD}/bin:$PATH

To  build for Xeon Phi using an Intel compiler:

    cmake -DCMAKE_CXX_COMPILER=icpc -DBUILD_PHI=ON -DCMAKE_BUILD_TYPE=Release ..
    make -j4

On systems like Stampede, required modules will need to be loaded:
 
    module load cmake intel cxx11

# Modes of Operation

`vargas -h`

```
vargas, May 11 2017. rgaddip1@jhu.edu
        define          Define a set of graphs for use with sim and align.
        sim             Simulate reads from a set of graphs.
        align           Align reads to a set of graphs.
        convert         Convert a SAM file to a CSV file.
        query           Convert a graph to DOT format.
        test            Run unit tests.
```

## define

`vargas define -h`

```
Define subgraphs deriving from a reference and VCF file.
Usage:
  vargas define [OPTION...]

  -h, --help  Display this message.

 Input options:
  -f, --fasta arg  <str> *Reference FASTA filename.

 Optional options:
  -v, --vcf arg       <str> Variant file (vcf, vcf.gz, or bcf).
  -t, --out arg       <str> Output filename. (default: stdout)
  -g, --region arg    <CHR[:MIN-MAX];...> list of regions. (default: all)
  -s, --subgraph arg  <str> Subgraph definitions, see below.
  -p, --filter arg    <str> Filter by sample names in file.
  -n, --limvar arg    <N> Limit to the first N variant records
  -c, --notcontig     VCF records for a given contig are not contiguous.


Subgraphs are defined using the format "label=N[%]",
where 'N' is the number of samples or percentage of samples to select.
The samples are selected from the parent graph, scoped with ':'.
The BASE graph is implied as the root for all labels. Example:
        a=50;a:b=10%;a:c=5
```

See [Define documentation](doc/define.md).

## align

`vargas align -h`

```
Align reads to a graph.
Usage:
  vargas align [OPTION...]

  -h, --help  Display this message.

 Input options:
  -g, --gdef arg   <str> *Graph definition file.
  -U, --reads arg  <str> *Unpaired reads in SAM, FASTQ, or FASTA format.

 Optional options:
  -S, --sam arg            <str> Output file.
      --msonly             Only report max score.
      --phred64            Qualities are Phred+64, not Phred+33.
  -p, --subsample arg      <N> Sample N random reads, 0 for all. (default: 0)
  -a, --alignto arg        <str> Target graph, or SAM Read Group -> graph
                           mapping."(RG:ID:<group>,<target_graph>;)+|<graph>"
  -s, --assess [=arg(=.)]  [ID] Use score profile from a previous alignment.
  -c, --tolerance arg      <N> Correct if within readlen/N. (default: 4)
  -f, --forward            Only align to forward strand.

 Scoring options:
      --ete      End to end alignment.
      --ma arg   <N> Match bonus. (default: 2)
      --mp arg   <MX,MN> Mismatch penalty. Lower qual=lower penalty.
                 (default: 6,2)
      --np arg   <N> Penalty for non-A/C/G/T. (default: 1)
      --rdg arg  <N1,N2> Read gap open/extension penalty. (default: 1,3)
      --rfg arg  <N1,N2> Ref gap open/extension penalty. (default: 1,3)

 Threading options:
  -j, --threads arg  <N> Number of threads. (default: 1)
  -u, --chunk arg    <N> Partition into tasks of max size N. (default: 64)
```

Reads are aligned to graphs specified in the GDEF file. `-x` will preform end to end alignment and is generally faster than full local alignment. The memory usage increase in marginal for high numbers of threads. As a result, as many threads as available should be used (271 on Xeon Phi KNL).

For example:

    vargas align  -g test.gdef -r reads.fa -t reads.sam -x

See the [Alignment documentation](doc/align.md) for more information.

## convert

`vargas convert -h`

```
Export a SAM file as a CSV file.
Usage:
  vargas convert [OPTION...]

  -f, --format arg  <str> *Output format.
  -s, --sam arg     <str,...> SAM files. Default stdin.
  -h, --help        Display this message.


Required column names:
        QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
Prefix with "RG:" to obtain a value from the associated read group.
```

Convert a SAM file into a CSV file, outputting the specified fields. Any of the SAM required fields or any aux tags can be output. For example,

```
vargas convert -s <samfile> -f "RG:ID,mp,ms"
```
will report the corresponding read group ID, max score position, and max score for each alignment. If multiple SAM files are provided, field 1 will be the file name. See [vargas align](doc/align.md) for tag information.

## sim

`vargas sim -h`

```
Simulate reads from genome graphs.
Usage:
  vargas sim [OPTION...]

 Required options:
  -g, --graph arg  <str> *Graph definition file.

 Optional options:
  -t, --out arg       <str> Output file. (default: stdout)
  -s, --sub arg       <S1,S2..> Subgraphs to simulate from. (default: base)
  -f, --file          -s specifies a filename.
  -l, --rlen arg      <N> Read length. (default: 50)
  -n, --numreads arg  <N> Number of reads to generate. (default: 1000)
  -j, --threads arg   <N> Number of threads. (default: 1)

 Stratum options:
  -d, --vnodes arg  <N1,N2...> Number of variant nodes. '*' for any. (default: *)
  -b, --vbases arg  <N1,N2...> Number of variant bases. '*' for any. (default: *)
  -m, --mut arg     <N1,N2...> Number of mutations. '*' for any. (default: 0)
  -i, --indel arg   <N1,N2...> Number of insertions/deletions. '*' for any. (default: 0)
  -a, --rate        Interpret -m, -i as error rates.

  -h, --help  Display this message.
```


`sim` generates `-n` reads of each combination of `-m`, `-i`, `-v`, and `-b`. `-m` Introduces mutation errors, substituting _N_ bases with an alternate base. Likewise, `-i` will delete a base or insert a random base. With `-a`, `-m` and `-i` are interpreted as rates (0.0 to 1.0). `-s` controls which subgraphs are used to generate reads.

For example:

    vargas sim -g test.gdef -t reads -n 1000 -m 0,1,2 -v 0,1,2 -j 12

will generate 1000 reads for each combination of `-m`, `-v`, for each graph in `test.gdef`.

Provided SAM tags:

- `ro` Unmutated read
- `nd` VCF sample simulated from
- `se` Number of substitution errors
- `vd` Number of variant nodes traversed
- `vb` Number of variant bases traversed
- `ni` Number of indel errors
- `gd` Read Group tag. Graph simulated from.
- `rt` Read Group tag. Rates were used.
- `ph` Read Group tag. GDEF file.

Reads are randomly sampled by weighting graph nodes by their length. As a result this isn't feasible for large graphs or specific conditions that are rare.

## query

`vargas query -h`

```
Query a graph and export a DOT graph.
Usage:
  vargas query [OPTION...]

  -g, --graph arg           *<str> Graph file to query.
  -d, --dot arg             <str> Subgraph to export as a DOT graph.
  -t, --out arg             <str> DOT output file. (default: stdout)
  -a, --stat [=arg(=base)]  <str> Print statistics about a subgraph.
  -h, --help                Display this message.
```

Export a subgraph to a DOT graph, or get graph statistics.

## Other

`vargas test` executes unit tests.

`vargas profile` Generates a summary of performance.

```
Run profiles. 
Usage:
  vargas profile [OPTION...]

  -f, --fasta arg    <str> *Reference FASTA.
  -v, --vcf arg      <str> *Variant file (vcf, vcf.gz, or bcf)
  -g, --region arg   <str> *Region of format "CHR:MIN-MAX". "CHR:0-0" for all.
  -i, --ingroup arg  <N> Ingroup percentage. (default: 100)
  -n, --nreads arg   <N> Number of reads. (default: 32)
  -l, --len arg      <N> Number of reads. (default: 50)
  -h, --help         Display this message.
```