#! /bin/bash

set -e -x

current_dir="$(pwd)"

./plot_validation_sim_results.py

./rasterize-model-plots.sh

cd ../results

for p in grid-*.tex;
do
    path_prefix="${p%.*}"
    pdf_path="${path_prefix}.pdf"
    cropped_pdf_path="cropped-${pdf_path}"
    latexmk -C "$p"
    latexmk -pdf "$p"
    pdfcrop "$pdf_path" "$cropped_pdf_path"
    latexmk -C "$p"
done

cd "$current_dir"
