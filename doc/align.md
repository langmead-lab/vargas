# Aligning
_Updated: May 11, 2017_

Before any alignment, a graph needs to be generated with [vargas define](define.md). To use default parameters and produce a SAM file,

```
vargas align -g <graph_def> -r <reads> -t <aligns_out.sam>
```

where `<reads>` can be a FASTA/Q or SAM file.

## Scoring options

- `-x` Use end to end alignment. This is faster than full local.
- `-m` Match bonus.
- `-n` Mismatch penalty (provide a positive number)
- `-o` Gap opening penalty (provide a positive number)
- `-e` Gap extension penalty (provide a positive number)

## Alignment targets

By default, reads will align to the base graph. If the input is in SAM format, specific read groups can be aligned to specific graphs with the `-a` option. The argument can either (a) be a graph to align to, or (b) A list of target mappings with the format

```
vargas align -g <graph_def> -r <reads> -t <out.sam> -a "RG:<tag>:<val>,<subgraph>; ..."
```

where `<tag>` can be `ID` or some auxiliary tag in the SAM Read Group line. RG's with `<val>` in the `<tag>` field will be aligned to `<subgraph>`. Reads that are not associated with a read group are assigned to `VAUGRP`.

## Assess

If a SAM read file is provided, `-s` can attempt to match a previous scoring function. Currently Bowtie2 and BWA MEM are supported. This also enables a few ancillary behaviors:

- If multiple positions are tied, the one closest to the previous alignment will be reported.
- `cf` and `ts` will report, as described below.
- `-c` controls the tolerance for deciding what is "correct"

## Output

Alignments are written as SAM files. Alignments are not written as primary alignments, but as aux tags. Note CIGAR alignments are not provided, just positions.

- `mp` Maximum scoring position.
- `sp` Second best scoring position.
- `ms` Maximum score.
- `ss` Second best score.
- `mc` Number of max score occurrences
- `sc` Number of second-best score occurrences
- `cf` Only in assess mode: 1 if alignment matches the best score, 2 for the second best, else 0.
- `ts` Only in assess mode: Score at the given alignment position.
- `pr` Scoring profile used.
- `mu` Sequence name of the best score.
- `su` Sequence name of the second best score.
- `gd` Subgraph aligned to.

`vargas convert` can be used to extract these fields into a CSV file.

## Coordinates

When a graph diverges, parallel nodes can have different lengths. To project the position onto the reference, nodes are anchored at the ends of nodes. For example:

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
