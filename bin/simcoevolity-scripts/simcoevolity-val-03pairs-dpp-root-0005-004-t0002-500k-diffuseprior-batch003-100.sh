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

simname="03pairs-dpp-root-0005-004-t0002-500k"
cfgpath="../../configs/config-${simname}.yml"
priorcfgpath="../../configs/config-03pairs-diffuse-prior-500k.yml"
outputdir="../../simulations/validation/${simname}-diffuseprior/batch003"
rngseed=528988595
nreps=100

mkdir -p "$outputdir"

simcoevolity --seed="$rngseed" -p "$priorcfgpath" -n "$nreps" -o "$outputdir" "$cfgpath"
