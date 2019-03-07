# markdup.py

## Purpose:

This script serves much the same purpose as `Picard MarkDuplicates` or `samtools
markdup` -- flag reads that are likely the result of PCR amplification.  I wrote
yet another markdup tool (YAMT) to handle duplicates in PacBio BAMs. I allow
for a few bases of "wiggle" on the ends of the alignment, with a default value
of 2 bp. If an `rq` tag is available, as in aligned CCS BAMs, I choose the read
with the highest `rq`.  Otherwise, as in aligned subreads BAMs, I choose the
first read that the algorithm encounters. The duplicate read frequency is output
to `stdout`.  If `--output <target>` is specified, a BAM is written to
`<target>` with duplicate reads flagged with the `0x400` bit.

## Requirements:

```bash
pip install pysam pandas
```

## Usage:

```bash
usage: markdup.py [-h] [--outBAM OUTBAM] [--wiggle WIGGLE] inBAM

Identify and flag PCR duplicate reads in an aligned PacBio BAM. For an aligned
CCS BAM, the alignment with the highest rq is chosen. For an aligned subreads
BAM, the first alignment encountered is chosen.

positional arguments:
  inBAM            aligned bam

 optional arguments:
   -h, --help       show this help message and exit
   --outBAM OUTBAM  output bam
   --wiggle WIGGLE  bp alignment wiggle at ends (2)
```
