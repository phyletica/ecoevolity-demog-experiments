#! /bin/sh

username="$USER"
if [ "$username" == "aubjro" ]
then
    module load gcc/6.1.0
fi

if [ -n "$PBS_JOBNAME" ]
then
    source ${PBS_O_HOME}/.bash_profile
    cd $PBS_O_WORKDIR

    module load gcc/5.3.0
fi

ssprob="0.40"
simname="03pops-dpp-root-0005-004-t0002-500k"
cfgpath="../../configs/config-${simname}.yml"
outputdir="../../simulations/validation/${simname}-040s/batch001"
rngseed=747229760
nreps=500

mkdir -p "$outputdir"

simcoevolity --seed="$rngseed" --singleton-sample-probability "$ssprob" -n "$nreps" -o "$outputdir" "$cfgpath"
