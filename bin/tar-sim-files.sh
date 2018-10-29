#!/bin/bash

bin_dir="$(pwd)"

for batch_dir in ../simulations/validation/03pops-dpp-root-0[05][05][05]*/batch*
do
    cd "$batch_dir"
    echo "$(pwd)"
    rm simcoevolity-sim-* var-only-simcoevolity-sim-*
    tar czf sim-files-run-1.tar.gz run-1-* && rm run-1-*
    tar czf sim-files-run-2.tar.gz run-2-* && rm run-2-*
    tar czf sim-files-run-3.tar.gz run-3-* && rm run-3-*
    tar czf sim-files-run-4.tar.gz run-4-* && rm run-4-*

    cd "$bin_dir"
done
