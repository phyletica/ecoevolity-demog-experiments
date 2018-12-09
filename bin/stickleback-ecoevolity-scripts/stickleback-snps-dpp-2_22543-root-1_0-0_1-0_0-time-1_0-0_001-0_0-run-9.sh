#! /bin/sh

if [ -n "$PBS_JOBNAME" ]
then
    source "${PBS_O_HOME}/.bash_profile"
    cd "$PBS_O_WORKDIR"
    module load gcc/5.3.0
fi

run_number=9
seed=73416868
config_prefix="stickleback-snps-dpp-2_22543-root-1_0-0_1-0_0-time-1_0-0_001-0_0"

output_dir="../../stickleback-ecoevolity-output"

if [ ! -d "$output_dir" ]
then
    mkdir "$output_dir"
fi

config_path="../../configs/${config_prefix}.yml"
prefix="${output_dir}/run-${run_number}-"

ecoevolity --seed "$seed"  --prefix "$prefix" --relax-missing-sites --relax-constant-sites --relax-triallelic-sites "$config_path" 1>"${prefix}${config_prefix}.out" 2>&1
