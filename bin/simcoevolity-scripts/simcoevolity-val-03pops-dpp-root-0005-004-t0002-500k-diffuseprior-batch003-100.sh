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

simname="03pops-dpp-root-0005-004-t0002-500k"
cfgpath="../../configs/config-${simname}.yml"
priorcfgpath="../../configs/config-03pops-diffuse-prior-500k.yml"
outputdir="../../simulations/validation/${simname}-diffuseprior/batch003"
rngseed=528988595
nreps=100

mkdir -p "$outputdir"

simcoevolity --seed="$rngseed" -p "$priorcfgpath" -n "$nreps" -o "$outputdir" "$cfgpath"


simname="03pops-dpp-root-0005-004-3_8-t0002-500k"
cfgpath="../../configs/config-${simname}.yml"
outputdir="../../simulations/validation/${simname}-diffuseprior/batch003"
rngseed=634077095

mkdir -p "$outputdir"

simcoevolity --seed="$rngseed" -p "$priorcfgpath" -n "$nreps" -o "$outputdir" "$cfgpath"
