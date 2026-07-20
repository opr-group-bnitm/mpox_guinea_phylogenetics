#!/usr/bin/env bash

set -euo pipefail
set -x

cat <<'EOF'
IMPORTANT:
This script runs Beast with the alignment that we previously prepared.
You can of course use your own alignment. Run the previous steps first.
The alignment can then e found in ../results/IIa/beast_in/cladeIIa.annotated.aln.fasta.
Load it into beauti and prepare a config file. Then adapt the path of the beast configs in this
script (fname_beast_config) and run this script.
EOF

config_path=beast_configs

for casename in strict_clock relaxed_clock
do
    fname_beast_config="${config_path}/${casename}.xml"

    for seed in 3 5 7 11 13 17
    do
        echo "Running BEAST for ${casename} with seed ${seed}"

        outdir="../results/IIa/beast/${casename}/SEED${seed}"
        prefix="$outdir/SEED${seed}"

        mkdir -p "$outdir"
        mkdir -p "$prefix"

        beast \
            -seed $seed \
            -overwrite \
            -threads 8 \
            -prefix "$prefix" \
            "$fname_beast_config"
    done
done
