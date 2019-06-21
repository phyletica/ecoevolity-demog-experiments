#! /bin/bash

set -e

paths="$(ls ../results/plots/*-nevents.pdf)"

for f in $paths
do
    filesuffix="${f##*-}"
    if [ "$filesuffix" = "rasterized.pdf" ]
    then
        continue
    fi
    
    echo "Rasterizing and compressing $f"
    n=${f/\.pdf/-rasterized\.pdf}
    convert -density 600 -compress jpeg -quality 80 $f $n
done
