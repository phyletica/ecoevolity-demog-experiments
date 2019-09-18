#!/bin/bash

set -e

bin_dir="$(pwd)"

for sim_dir in ../simulations/validation/0[36]p*/batch*
do
    cd "$sim_dir"
    if [ -e "simcoevolity-model-used-for-sims.yml" ]
    then
        pwd
        rm *.nex
        rm *.sh
        rm *.sh.o*
        rm *simcoevolity-sim-*yml
        tar czf sim-files-run-1.tar.gz run-1-* && rm run-1-*
        tar czf sim-files-run-2.tar.gz run-2-* && rm run-2-*
        tar czf sim-files-run-3.tar.gz run-3-* && rm run-3-*
        tar czf sim-files-run-4.tar.gz run-4-* && rm run-4-*
        tar czf sim-files-true-values.tar.gz simcoevolity-sim-???-true-values.txt simcoevolity-model-used-for-sims.yml && rm simcoevolity-sim-???-true-values.txt simcoevolity-model-used-for-sims.yml
    else
        echo "Skipping $(pwd)"
    fi
    cd "$bin_dir"
done

cd "$bin_dir"
