#!/usr/bin/env bash

usage() {
    echo "Usage: $0 -i INPUT_DIR -s SUFFIX -b BURNIN -d SEEDS"
    echo
    echo "Options:"
    echo "  -i    Input/output base directory (e.g. /path/to/data)"
    echo "  -s    Suffix used in filenames (e.g. mpox_cladeIIa_2025_03_05)"
    echo "  -b    Burn-in value for logcombiner (e.g. 10)"
    echo "  -d    Comma-separated list of seed numbers (e.g. 3,5,7)"
    echo
    exit 1
}

# Parse options
while getopts ":i:s:b:d:" opt; do
    case $opt in
        i) datadir="$OPTARG" ;;
        s) suffix="$OPTARG" ;;
        b) burnin="$OPTARG" ;;
        d) IFS=',' read -ra seeds <<< "$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if all required arguments are provided
if [[ -z "$datadir" || -z "$suffix" || -z "$burnin" || ${#seeds[@]} -eq 0 ]]; then
    usage
fi

combine_logs() {
    local datadir="$1"
    local suffix="$2"
    local burnin="$3"
    shift 3
    local seeds=("$@")

    local log_cmd="logcombiner -b $burnin"
    local trees_cmd="logcombiner -b $burnin"

    for seed in "${seeds[@]}"; do
        log_cmd+=" -log $datadir/SEED${seed}/SEED${seed}${suffix}.log"
        trees_cmd+=" -log $datadir/SEED${seed}/SEED${seed}${suffix}.trees"
    done

    log_cmd+=" -o $datadir/combined.log"
    trees_cmd+=" -o $datadir/combined.trees"

    eval "$log_cmd"
    eval "$trees_cmd"
}

combine_logs "$datadir" "$suffix" "$burnin" "${seeds[@]}"
