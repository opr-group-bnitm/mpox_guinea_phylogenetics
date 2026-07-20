#!/usr/bin/env bash

set -euo pipefail
set -x


fname_msa=../results/IIb/squirrel/sequences.aln.fasta
outgroup="KJ642615|human|Nigeria||1978"

outdir=../results/IIb/iqtree
prefix=$outdir/model_comparison

mkdir -p $outdir

iqtree3 \
    -pre $prefix -s $fname_msa \
    -m MFP+MERGE -mset JC,HKY,TN93,GTR -mrate G,E -mfreq FU,F \
    -bb 1000 -bnni \
    --polytomy -czb -blmin 1e-10 \
    --ancestral \
    -T AUTO

# reroote without removing the outgroup for the supplement plot
gotree reroot outgroup \
    -i $prefix.treefile $outgroup > $prefix.rooted.with_outgroup.treefile

sed -E 's/Node([0-9]+)\/([0-9]+)/Node\1_\2/g' \
    $prefix.rooted.with_outgroup.treefile \
    > $prefix.4baltic.with_outgroup.proto.treefile

python ../utils/fix_for_baltic.py \
    --in $prefix.4baltic.with_outgroup.proto.treefile \
    --out $prefix.4baltic.with_outgroup.treefile

# reroot with outgroup + afterwards remove it from the tree
gotree reroot outgroup \
    -i $prefix.treefile "${outgroup}" \
    --remove-outgroup > $prefix.rooted.treefile

sed -E 's/Node([0-9]+)\/([0-9]+)/Node\1_\2/g' \
    $prefix.rooted.treefile \
    > $prefix.4baltic.proto.treefile

python ../utils/fix_for_baltic.py \
    --in $prefix.4baltic.proto.treefile \
    --out $prefix.4baltic.treefile

python ../utils/assign_snps.py \
    --ancestral "$prefix.state" \
    --tree "$prefix.4baltic.treefile" \
    --msa "$fname_msa" \
    --output "$prefix.apobec_counts.csv" \
    --table "$prefix.snps.csv"
