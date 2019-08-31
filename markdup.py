#!/usr/bin/env python
"""Identify and flag PCR duplicate reads in an aligned PacBio BAM.

Reads with identical alignments within aln_wiggle on both ends
and with length within len_wiggle % are marked as duplicates.
Alignment with highest read quality, number of passes, and
query name md5 hash (in that order) is chosen as primary.
"""
import argparse
import hashlib
import pysam
import pandas as pd
import numpy as np


def mean_read_qual(query_qualities):
    "convert phred scaled qual array to mean read quality"
    return 1 - 10 ** (np.mean(query_qualities)/-10)


def match_ref(df, i):
    "i is aligned to same reference as i-1"
    return df.loc[i, 'ref_id'] == df.loc[i-1, 'ref_id']


def start_match(df, i, aln_wiggle):
    "i starts within aln_wiggle bp of i-1"
    return abs(df.loc[i-1, 'ref_start'] - df.loc[i, 'ref_start']) <= aln_wiggle


def end_match(df, i, aln_wiggle):
    "i ends within aln_wiggle bp of i-1"
    return abs(df.loc[i, 'ref_end'] - df.loc[i-1, 'ref_end']) <= aln_wiggle


def length_match(df, i, len_wiggle):
    "i query_length within len_wiggle bp of i-1"
    return abs(df.loc[i-1, 'query_length'] - df.loc[i, 'query_length']) <= \
        len_wiggle/100.0 * df.loc[i, 'query_length']


def set_primary(df, dup_index):
    ""
    # if rq and np, use top rq, top np, top md5(name)
    # if rq and not np, use top rq, top md5(name)
    # else use top md5(name)
    if 'read_qual' in fields and 'num_passes' in fields:
        by = ['read_qual', 'num_passes', 'query_md5']
        ascending = [False, False, True]
    elif 'read_qual' in fields:
        by = ['read_qual', 'query_md5']
        ascending = [False, True]
    else:
        by = ['query_md5']
        ascending = [True]
    # set the top sorted alignment as the primary
    df.loc[df[df['dup_index'] == dup_index]
           .sort_values(by=by, ascending=ascending).index[0], 'dup'] = False


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('inBAM', help='aligned bam', type=str)
parser.add_argument('--outBAM', help='output bam', type=str)
parser.add_argument('--aln_wiggle', help='bp alignment wiggle at ends, default 2',
                    type=int, default=2)
parser.add_argument('--len_wiggle', help='percent read len wiggle, default 10',
                    type=int, default=10)
args = parser.parse_args()

# fields to store for each alignment
fields = ['query_md5',
          'ref_id',
          'ref_start',
          'ref_end',
          'query_length',
          'read_qual',
          'num_passes',
          'dup',
          'dup_index']
data = {f: list() for f in fields}

# iterate through reads, saving coordinates and read quality
with pysam.AlignmentFile(args.inBAM, 'rb') as infile:
    for read in infile.fetch():
        data['query_md5'].append(hashlib.md5(read.query_name).hexdigest())
        data['ref_id'].append(int(read.reference_id))
        data['ref_start'].append(int(read.reference_start))
        data['ref_end'].append(int(read.reference_end))
        data['query_length'].append(int(read.infer_query_length()))
        if read.has_tag('rq'):
            data['read_qual'].append(float(read.get_tag('rq')))
        elif read.query_qualities:
            data['read_qual']\
                .append(mean_read_qual(read.query_qualities))
        if read.has_tag('np'):
            data['num_passes'].append(int(read.get_tag('np')))
        data['dup'].append(False)
        data['dup_index'].append(None)
if not data['read_qual']:
    data.pop('read_qual', None)
    fields.remove('read_qual')
if not data['num_passes']:
    data.pop('num_passes', None)
    fields.remove('num_passes')
df = pd.DataFrame(data)

# iterate through the dataframe, marking duplicate alignments
dup_index = 0
dup_state = False
for i in range(1, len(df)):
    # two reads aligned to the same ref chrom are duplicate if:
    # - alignments start within aln_wiggle bp of each other
    # - alignments end within aln_wiggle bp of each other
    # - read lengths with len_wiggle bp of each other
    if match_ref(df, i) and \
       start_match(df, i, args.aln_wiggle) and \
       end_match(df, i, args.aln_wiggle) and \
       length_match(df, i, args.len_wiggle):
        df.loc[[i, i-1], 'dup'] = True  # mark this read and previous read
        df.loc[[i, i-1], 'dup_index'] = dup_index  # a block of duplicates
        dup_state = True  # in the middle of a block of duplicates
    elif dup_state:
        dup_state = False
        dup_index += 1

# iterate through duplicate sets, and unmark the primary alignment
for i in range(0, dup_index):
    set_primary(df, i)

duplicates = set(df[df['dup']].index.values)
# iterate through bamfile again, marking all duplicates with dup flag
if args.outBAM:
    with pysam.AlignmentFile(args.inBAM, 'rb') as infile, \
         pysam.AlignmentFile(args.outBAM, 'wb', template=infile) as outfile:
        for i, read in enumerate(infile.fetch()):
            if i in duplicates:
                read.is_duplicate = True
            outfile.write(read)

print('{:.08f}'.format((len(duplicates)/float(len(df)))))
