#!/usr/bin/env python

import os
import pandas as pd
from pybedtools import BedTool
from itertools import combinations

# user-defined variables
window = 5000  #(pads either side)

# file paths
workdir = "/home/dad5925/workdir/"
scriptdir = "/home/dad5925/scripts/"
exon_bed = os.path.join(workdir, "hexevent-all-chr-cassette-0-1-100-100-strict-20210622.bed.sorted")
alus_bed = os.path.join(workdir, "hg38_fasta/family/SINE_Alu")

# start with single window size, then expand to other sizes
processed_dir = "/home/dad5925/processing/20210709-bedops/"
pairs_bed = os.path.join(processed_dir, "closest-features--dist-exon-alu.tsv")
pairs_window_bed = os.path.join(processed_dir, "closest-features--dist--no-overlaps-exon-alu-5000.tsv")
inverted_pairs_window_bed = os.path.join(processed_dir, "closest-features--dist--no-overlaps-exon-alu-5000-opposite-strands.tsv")
#closest_exon_window_bed = os.path.join(processed_dir, "closest-features--dist--closest-alu-exon-5000.tsv")
closest_exon_bed = os.path.join(processed_dir, "closest-features--dist--closest-alu-exon.tsv")

# read in with multiple field delimiters
bed6_cols = ['chr', 'start', 'end', 'gene', 'score', 'strand']
alu_bed_cols = ['chr', 'start', 'end', 'subfamily', 'sw_score', 'strand', 'percent_substitution', 'percent_deleted', 'percent_inserted', 'num_bases_past_end', 'family', 'p1', 'p2', 'p3', 'id', 'dist']
pairs_df_header = ["exon_{}".format(i) for i in bed6_cols] + ["upstream_alu_{}".format(i) for i in alu_bed_cols] + ["downstream_alu_{}".format(i) for i in alu_bed_cols]
closest_exon_header = ["alu_{}".format(i) for i in alu_bed_cols[:-1]] + ["exon_{}".format(i) for i in bed6_cols] + ["dist"] 
dtypes_cols = ["exon_{}".format(i) for i in ['start', 'end']] + \
              ["upstream_alu_{}".format(i) for i in ['start', 'end', 'dist', 'id']] + \
              ["downstream_alu_{}".format(i) for i in ['start', 'end', 'dist', 'id']]
dtypes_dct = { i:"int" for i in dtypes_cols}  #XXX not working... probably bc lots of NaN for chrUn etc

pairs_df = pd.read_csv(pairs_bed, sep='[\t|]', engine='python', header=None, names=pairs_df_header)  #dtype=dtypes_dct
pairs_window_df = pd.read_csv(pairs_window_bed, sep='[\t|]', engine='python', header=None, names=pairs_df_header)
inverted_pairs_window_df = pd.read_csv(inverted_pairs_window_bed, sep='[\t|]', engine='python', header=None, names=pairs_df_header) 
#closest_exon_window_df = pd.read_csv(closest_exon_window_bed, sep='[\t|]', engine='python', header=None, names=closest_exon_header)
closest_exon_df = pd.read_csv(closest_exon_bed, sep='[\t|]', engine='python', header=None, names=closest_exon_header)

# find *all* Alus that are within window_size of exon, and try all combinations of pairs
#closest_exon_df[closest_exon_df['dist'].between(-window,window)].groupby(['exon_chr','exon_start','exon_end']).filter(lambda x: x['alu_strand'] != x['exon_strand']).size()

window_df = closest_exon_df[closest_exon_df['dist'].between(-window,window)].dropna(how='any')  # dropna doesn't change anything once window is applied
window_df = window_df.rename_axis('idx').reset_index().astype({'dist': 'int'})

keep_exon_cols = ['exon_chr','exon_start','exon_end', 'exon_gene', 'exon_score', 'exon_strand']
# also below: enumerate(x.values)
grouped_df = window_df.groupby(keep_exon_cols, sort=False)[['alu_strand', 'idx']].apply(lambda x: list(combinations(x.values,2))).apply(pd.Series).stack().reset_index(name='strand_pairs')
grouped_df[['upstream', 'downstream']] = pd.DataFrame(grouped_df['strand_pairs'].tolist(), columns=['upstream', 'downstream'])  # XXX is it for sure upstream/downstream order
# check if bedops takes into account exon strand when computing signed distance
grouped_df[['upstream_strand', 'upstream_id']] = pd.DataFrame(grouped_df['upstream'].tolist(), columns=['upstream_strand', 'upstream_idx'])
grouped_df[['downstream_strand', 'downstream_id']] = pd.DataFrame(grouped_df['downstream'].tolist(), columns=['downstream_strand', 'downstream_idx'])
grouped_df = grouped_df.astype({'exon_start': 'int', 'exon_end': 'int'})

# get dist and other cols based off of idx
ir_df = grouped_df[(grouped_df.upstream_strand != grouped_df.downstream_strand)]

cols = ['alu_subfamily','alu_id', 'dist']
rename_upstream_cols = dict(zip(cols,["upstream_{}".format(i) for i in cols]))
rename_downstream_cols = dict(zip(cols,["downstream_{}".format(i) for i in cols]))

ir_df = ir_df.merge(window_df[['idx'] + cols]
                    .rename(columns={'idx': 'upstream_id'} | rename_upstream_cols),
                    on='upstream_id', how='left')
ir_df = ir_df.merge(window_df[['idx'] + cols]
                    .rename(columns={'idx': 'downstream_id'} | rename_downstream_cols),
                    on='downstream_id', how='left')

# filter on pairs with flank the exon
flank_ir_df = ir_df[(ir_df['upstream_dist']>0) & (ir_df['downstream_dist']<0)]

# Notes/Todo:
# to check: window_df[window_df['exon_start'] == 92090.0]
# e.g. 4 entries in df (same exon start stop), resulting in 6 combinations (4 choose 2 pairs)

# filter on pairs that are within window_size and inverted
### are inverted enriched near skipped exons compared to non-inverted pairs? across window sizes
### are inverted Alus enriched compared to other types of inversions?
# make plots (plots should go in different script and be called in here)


