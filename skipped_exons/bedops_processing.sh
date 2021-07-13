#!/usr/bin/env bash

# Script for BEDOPS processing of skipped exons and Alus

# Note: Load the "leelab" Conda environment beforehand
#conda deactivate && conda activate leelab

# TODO: add nicer argparse
WINDOW=50  #added upstream AND downstream to exon for search

WORKDIR="/home/dad5925/workdir/"
EXON_BED="${WORKDIR}/hexevent-all-chr-cassette-0-1-100-100-strict-20210622.bed.sorted"  # reference file
ALU_BED="${WORKDIR}/hg38_fasta/family/SINE_Alu"  # query file
BEDOPS_OUTDIR="/home/dad5925/processing/20210709-bedops/"

# find nearest Alu pair that lies both upstream and downstream
sort-bed $EXON_BED | closest-features --dist - $ALU_BED > ${BEDOPS_OUTDIR}/closest-features--dist-exon-alu.tsv  # total: 7439, same with --no-overlaps

# find nearest Alu pair that lies both upstream and downstream and that is within 10kb window
# $1: ref, $2: upstream Alu, $3: signed distance (upstream), $4: downstream Alu, $5: signed distance (downstream)
sort-bed $EXON_BED | closest-features --no-overlaps --dist - $ALU_BED | awk -v var="$WINDOW" -F'|' '$3 >= -var && $3 <= var' | awk -v var="$WINDOW" -F'|' '$5 >= -var && $5 <= var' > ${BEDOPS_OUTDIR}/closest-features--dist--no-overlaps-exon-alu-${WINDOW}.tsv  # for window of 5000, total: 5287

# strand (Alu orientation)
# if first entry + and second - or vice versa (since always printed in upstream/downstream fashion)
sort-bed $EXON_BED | closest-features --no-overlaps --dist - $ALU_BED | awk -v var="$WINDOW" -F'|' '$3 >= -var && $3 <= var' | awk -v var="$WINDOW" -F'|' '$5 >= -var && $5 <= var' | awk -F'[|\t]' '($12 == "+" && $28 == "-") || ($12 == "-" && $28 == "+")' > ${BEDOPS_OUTDIR}/closest-features--dist--no-overlaps-exon-alu-${WINDOW}-opposite-strands.tsv  # total: 2521
# XXX check orientation relative to exon (gene)

# IR Alus quick total counts:
### 5kb window: 2521
### 4kb window: 2271
### 3kb window: 1945
### 2kb window: 1457
### 1kb window: 741
### 500b window: 318
### 250b window: 101
### 100b window: 12
### 50b window: 3

# find *all* Alus within window upstream
# find *all* Alus within window downstream

# closest-features --dist --nearest alu exon | cutoff window # each alu will have one associated exon # check that each alu to exon signed distance is less than window
sort-bed $ALU_BED | closest-features --dist --closest - $EXON_BED | awk -F'|' -v var="$WINDOW" '$3 >= -var && $3 <= var' > ${BEDOPS_OUTDIR}/closest-features--dist--closest-alu-exon-${WINDOW}.tsv  # total: 38247; --no-overlaps total: 37431
wc -l $ALU_BED  # total: 1262524
sort-bed $ALU_BED | closest-features --dist --closest - $EXON_BED > ${BEDOPS_OUTDIR}/closest-features--dist--closest-alu-exon.tsv

# above, opposite strands

