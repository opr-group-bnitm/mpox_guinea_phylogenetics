#!/usr/bin/env bash

set -euo pipefail
set -x


# Checked convergence before with tracer and determined burn-in rate to
# get large enough ESS for all 6 runs.

for case in strict_clock relaxed_clock
do
    beast_results_dir="../results/IIa/beast/${case}"

    burnin=10

    bash ../utils/combine_logs.sh \
        -i "${beast_results_dir}" \
        -s cladeIIa \
        -b $burnin \
        -d 3,5,7,11,13,17

    # build maximum clade credibility tree with common ancestor heights

    heights_measure=median

    export _JAVA_OPTIONS='-Xms2g -Xmx40g'

    treeannotator \
        -burnin 0 \
        -heights $heights_measure \
        "${beast_results_dir}/combined.trees" \
        "${beast_results_dir}/mcc_${heights_measure}.tree"
done
