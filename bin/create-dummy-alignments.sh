#!/bin/bash

NSPECIES=1
NGENOMES=20
NCHARS=500000

alignment_dir="../alignments"
mkdir -p "$alignment_dir"

# Generate 500k dummy alignments
i=1
while [ "$i" -lt 4 ]
do
    prefix="c${i}sp"
    outfile="${alignment_dir}/comp0${i}-${NSPECIES}species-${NGENOMES}genomes-${NCHARS}chars.nex"
    ./generate-dummy-biallelic-alignment.py \
        --nspecies "$NSPECIES" \
        --ngenomes "$NGENOMES" \
        --ncharacters "$NCHARS" \
        --prefix "$prefix" \
        > "$outfile"
    i=`expr $i + 1`
done

# Generate 100k dummy alignments
NCHARS=100000

i=1
while [ "$i" -lt 4 ]
do
    prefix="c${i}sp"
    outfile="${alignment_dir}/comp0${i}-${NSPECIES}species-${NGENOMES}genomes-${NCHARS}chars.nex"
    ./generate-dummy-biallelic-alignment.py \
        --nspecies "$NSPECIES" \
        --ngenomes "$NGENOMES" \
        --ncharacters "$NCHARS" \
        --prefix "$prefix" \
        > "$outfile"
    i=`expr $i + 1`
done

# Generate 500k dummy alignments with missing data
NCHARS=500000

for m in 0.1 0.25 0.5
do
    mlabel=${m/\./}
    i=1
    while [ "$i" -lt 4 ]
    do
        prefix="c${i}sp"
        outfile="${alignment_dir}/comp0${i}-${NSPECIES}species-${NGENOMES}genomes-${NCHARS}chars-${mlabel}missing.nex"
        ./generate-dummy-biallelic-alignment.py \
            --nspecies "$NSPECIES" \
            --ngenomes "$NGENOMES" \
            --ncharacters "$NCHARS" \
            --prefix "$prefix" \
            --missing-probability "$m" \
            > "$outfile"
        i=`expr $i + 1`
    done
done
