#!/usr/bin/env bash

# BEAST analysis needs a time signal. Test if we find a significant one.

tree=../results/IIa/iqtree/model_comparison.rooted.treefile
outdir=../results/IIa/time_signal

mkdir -p "$outdir"

python ../utils/test_time_signal.py \
    --tree "$tree" \
    --meta ../data/IIa/meta_all.tsv \
    --plot_out "${outdir}/time_signal_test.png" \
    | tee "${outdir}/time_signal_significance.txt"
