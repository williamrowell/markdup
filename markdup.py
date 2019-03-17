#!/usr/bin/env python
"""Identify and flag PCR duplicate reads in an aligned PacBio BAM.

For an aligned CCS BAM, the alignment with the highest rq is chosen.
For an aligned subreads BAM, the first alignment encountered is chosen.
"""
import argparse
import pysam
import pandas as pd


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('inBAM', help='aligned bam', type=str)
parser.add_argument('--outBAM', help='output bam', type=str)
parser.add_argument('--wiggle', help='bp alignment wiggle at ends (2)',
                    type=int, default=2)
args = parser.parse_args()

# iterate through reads, saving coordinates and read quality
data = {'ref_id': list(),
        'ref_start': list(),
        'ref_end': list(),
        'read_qual': list(),
        'dup': list(),
        'dup_index': list()}

with pysam.AlignmentFile(args.inBAM, 'rb') as infile:
    for read in infile.fetch():
        data['ref_id'].append(int(read.reference_id))
        data['ref_start'].append(int(read.reference_start))
        data['ref_end'].append(int(read.reference_end))
        if read.has_tag('rq'):  data['read_qual'].append(float(read.get_tag('rq')))
        data['dup'].append(False)
        data['dup_index'].append(None)
if not data['read_qual']: data.pop('read_qual', None)
df = pd.DataFrame(data)

# iterate through this dataframe, marking duplicate reads
dup_index = 0
dup_state = False
for i in range(1, len(df)):
    if df.loc[i, 'ref_id'] == df.loc[i-1, 'ref_id'] and \
       abs(df.loc[i-1, 'ref_start'] - df.loc[i, 'ref_start']) <= args.wiggle and \
       abs(df.loc[i, 'ref_end'] - df.loc[i-1, 'ref_end']) <= args.wiggle:
        df.loc[[i, i-1], 'dup'] = True  # mark this read and previous read
        df.loc[[i, i-1], 'dup_index'] = dup_index  # a block of duplicate alignments
        dup_state = True  # we're in the middle of a block of duplicates
    elif dup_state:
        dup_state = False
        dup_index += 1

# iterate through duplicate sets, and unmark the read with the highest quality
# if quality is not available (subreads), choose one at random
if 'read_qual' in data:
    for i in range(0, dup_index):
        df.loc[df[df['dup_index'] == i]['read_qual'].idxmax(), 'dup'] = False
else:
    for i in range(0, dup_index):
        df[df['dup_index'] == i].sample(axis=0)['dup'] = False

duplicates = set(df[df['dup']].index.values)
# iterate through bamfile again, marking all duplicates with dup flag
if args.outBAM:
    with pysam.AlignmentFile(args.inBAM, 'rb') as infile, \
         pysam.AlignmentFile(args.outBAM, 'wb', template=infile) as outfile:
        for i, read in enumerate(infile.fetch()):
            if i in duplicates:
                read.is_duplicate = True
            outfile.write(read)

print '{:.08f}'.format((len(duplicates)/float(len(df))))
