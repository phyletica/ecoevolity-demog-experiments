#!/bin/bash

set -e -x

if [ -n "$PBS_JOBNAME" ]
then
    source ${PBS_O_HOME}/.bash_profile
    cd $PBS_O_WORKDIR
fi

current_dir="$(pwd)"

project_dir="$(python ../project_util.py)"
ecoevolity_output_dir="${project_dir}/stickleback-ecoevolity-output"
cd "$ecoevolity_output_dir"

for datatype in "allsites" "snps"
do
    for dpp in "0_3725" "2_22543" "13_0"
    do
        for t in "1_0-0_001-0_0" "1_0-0_0005-0_0" "1_0-0_005-0_0"
        do
            for root in "1_0-1_0-0_0" "1_0-0_5-0_0" "1_0-0_1-0_0"
            do
                gzip -d -k run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log.gz
                echo "pyco-sumchains ${datatype}run-?-africa-codiv${suffix}dpp-${dpp}-time-${t}-state-run-1.log"
                pyco-sumchains run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log 1>pyco-sumchains-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-table.txt 2>pyco-sumchains-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-stderr.txt
                rm run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
            done
        done
    done
done

cd "$current_dir"
