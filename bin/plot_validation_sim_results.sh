#! /bin/bash

set -e -x

current_dir="$(pwd)"

./plot_validation_sim_results.py

cd ../results

for p in grid-*.tex;
do
    latexmk -C "$p"
    latexmk -pdf "$p"
    pdfcrop "$p" "cropped-${p}"
    latexmk -C "$p"
done

cd "$current_dir"
