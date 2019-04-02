# markdup.py

## Purpose

This is yet another markdup tool (YAMT) that serves much the same purpose as 
`Picard MarkDuplicates` or `samtools markdup` -- flag reads that are likely 
the result of PCR amplification.

For each alignment in an aligned, sorted BAM, I store reference id, alignment
start, alignment end, read quality (priority: `rq` tag, mean of quality string, 
or dropped if missing), number of passes (`np` tag, dropped if missing), query
length, and the md5 hash of the query name.

Reads are classified as duplicate if:

1. reference id matches
2. alignment starts and ends within aln_wiggle bp
3. read length is within len_wiggle %

I choose a single alignment to represent each set of duplicates by the following
rules in descending priority:

1. highest read quality
2. highest number of passes
3. lowest query name md5 hash

The duplicate read frequency is output to `stdout`.  If `--output <target>` is 
specified, a BAM is written to `<target>` with duplicate reads flagged with the 
`0x400` bit.

## Requirements

```bash
pip install pysam pandas numpy
```

## Usage

```bash
usage: markdup.py [-h] [--outBAM OUTBAM]
                  [--aln_wiggle ALN_WIGGLE]
                  [--len_wiggle LEN_WIGGLE]
                  inBAM

Reads with identical alignments within aln_wiggle on both ends
and with length within len_wiggle % are marked as duplicates.
Alignment with highest read quality, number of passes, and 
query name md5 hash (in that order) is chosen as primary.

positional arguments:
  inBAM                 aligned bam

optional arguments:
  -h, --help                show this help message and exit
  --outBAM OUTBAM           output bam
  --aln_wiggle ALN_WIGGLE   bp alignment wiggle at ends, default 2
  --len_wiggle LEN_WIGGLE   percent read len wiggle, default 10
```
