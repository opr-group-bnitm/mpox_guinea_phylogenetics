#!/usr/bin/env bash


fname_seqs=../data/IIb/sequences.fasta
outdir=../results/IIb/squirrel

mkdir -p "$outdir"

squirrel "$fname_seqs" \
    --clade cladeiib \
    -o $outdir \
    --threads 15

echo "DONE"
