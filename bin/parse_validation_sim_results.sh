#! /bin/sh

if [ -n "$PBS_JOBNAME" ]
then
    source ${PBS_O_HOME}/.bash_profile
    cd $PBS_O_WORKDIR
fi

./parse_validation_sim_results.py -r 3 -s 1501 --burnin 401 1> parse_validation_sim_results.out 2>&1
