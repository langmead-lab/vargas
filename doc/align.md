# Aligning
_Updated: May 11, 2017_

Before any alignment, a graph needs to be generated with [vargas define](define.md). To use default parameters and produce a SAM file,

```
vargas align -g <graph_def> -U <reads> -S <aligns_out.sam>
```

where `<reads>` can be a FASTA/Q or SAM file.

## Scoring options

- `--ete` Use end to end alignment. This is faster than full local.
- `--ma <INT>` Match bonus.
- `--mp <INT>` Mismatch penalty (provide a positive number).
- `--rdg <INT,INT>` Read gap penalty <Open,Extend>.
- `--rfg <INT,INT>` Reference gap penalty <Open,Extend>.
- `--np <INT>` Penalty for non-A/C/G/T bases.

## Alignment targets

By default, reads will align to the base graph. If the input is in SAM format, specific read groups can be aligned to specific graphs with the `-a` option. The argument can either (a) be a graph to align to, or (b) A list of target mappings with the format:

```
vargas align -g <graph_def> -r <reads> -t <out.sam> -a "RG:<tag>:<val>,<subgraph>; ..."
```

where `<tag>` can be `ID` or some auxiliary tag in the SAM Read Group line. RG's with `<val>` in the `<tag>` field will be aligned to `<subgraph>`. Reads that are not associated with a read group are assigned to `VAUGRP`.

Using a SAM input where an alignment is already defined will enable the reporting of the `cf` and `ts` flags.

## Assess

If a SAM read file is provided, `-s` can attempt to match a previous scoring function. Currently Bowtie2, HISAT2, and BWA MEM are supported.


## Output

Alignments are written as SAM files. Note CIGAR alignments are not provided. Fields relevant to Vargas are listed below.

- `POS` Maximum scoring alignment, only known for end to end alignments.
- `RNAME` Maximum scoring sequence.
- `FLAG` Reverse complement flag bit if aligned to opposite strand.
- `AS` Maximum score.
- `mp` Position of the maximum scoring cell.
- `ss` Second best score.
- `st` Strand of second best score.
- `sp` Position of the second highest scoring cell.
- `su` Sequence name of the second best score.
- `mc` Number of max-score occurrences.
- `sc` Number of second-best score occurrences.
- `gd` Read group tag. Subgraph aligned to.

`vargas convert` can be used to extract these fields into a CSV file.

## Coordinates

When a graph diverges, parallel nodes can have different lengths. To project the position onto the reference, nodes are anchored to the end of the node. For example:

```
    012 1234 5678
        ||||
        ACAC
       /    \
    AAC--GT--AAAA
         ||
    012  34  5678
```
Note how `C` and `T`, the ends of the nodes, have the same position.
