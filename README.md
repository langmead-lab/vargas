# Vargas
_Updated: May 11, 2017_

[![Build Status](https://travis-ci.org/RaviGaddipati/vargas.svg?branch=master)](https://travis-ci.org/RaviGaddipati/vargas)

Vargas aligns short reads to a directed acyclic graph (DAG). Reads are aligned using a SIMD vectorized version of Smith-Waterman, aligning multiple reads at once. The the aligner reports the best positions of reads without traceback to keep memory usage manageable. Additionally, vargas can be used to simulate reads from a DAG, and reads can be stratified by various parameters.


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

# Modes of operation

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

 Required options:
  -f, --fasta arg  <str> *Reference FASTA filename.

 Optional options:
  -v, --vcf arg       <str> Variant file (vcf, vcf.gz, or bcf).
  -t, --out arg       <str> Output filename. (default: stdout)
  -g, --region arg    <CHR[:MIN-MAX];...> CSV list of regions. (default: all)
  -s, --subgraph arg  <str> Subgraph definitions, see below.
  -p, --filter arg    <str> Filter by sample names in file.

  -h, --help  Display this message.
```

See [Define documentation](doc/define.md).

## align

`vargas align -h`

```
Align reads to a graph.
Usage:
  vargas align [OPTION...]

 Required options:
  -g, --gdef arg   <str> *Graph definition file.
  -r, --reads arg  <str> *Unpaired reads in SAM, FASTQ, or FASTA format.

 Optional options:
  -t, --out arg            <str> Output file.
  -p, --subsample arg      <N> Sample N random reads, 0 for all. (default: 0)
  -a, --alignto arg        <str> Target graph, or SAM Read Group -> graph mapping.
                           "(RG:ID:<group>,<target_graph>;)+|<graph>"
  -s, --assess [=arg(=-)]  [ID] Use score profile from a previous alignment,
                           and target nearby alignments.
  -c, --tolerance arg      <N> Correct if within readlen/N. (default: 4)

 Scoring options:
  -m, --match arg       <N> Match score. (default: 2)
  -n, --mismatch arg    <N> Mismatch penalty. (default: 2)
  -o, --gap_open arg    <N> Gap opening penalty. (default: 3)
  -e, --gap_extend arg  <N> Gap extension penalty. (default: 1)
  -x, --endtoend        Perform end to end alignment

 Threading options:
  -j, --threads arg  <N> Number of threads. (default: 1)
  -u, --chunk arg    <N> Partition tasks into chunks with max size N. (default: 64)

  -h, --help  Display this message.
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