#!/usr/bin/env bash

# compute a multiple sequence alignment with squirrel
set -euo pipefail

outgroup="NC_063383.1"
outdir=../results/IIa/squirrel
input_dir=../data/IIa

mkdir -p "${outdir}"
prefix=sequences_all

squirrel "${input_dir}/${prefix}.fasta" \
    --clade cladeiia \
    -o "$outdir" \
    --outgroup "$outgroup" \
    --threads 15
