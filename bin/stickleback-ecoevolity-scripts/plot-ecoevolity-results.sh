#!/bin/bash

set -e -x

if [ -n "$PBS_JOBNAME" ]
then
    source ${PBS_O_HOME}/.bash_profile
    cd $PBS_O_WORKDIR
fi

label_array=()
convert_labels_to_array() {
    local concat=""
    local t=""
    label_array=()

    for word in $@
    do
        local len=`expr "$word" : '.*"'`

        [ "$len" -eq 1 ] && concat="true"

        if [ "$concat" ]
        then
            t+=" $word"
        else
            word=${word#\"}
            word=${word%\"}
            label_array+=("$word")
        fi

        if [ "$concat" -a "$len" -gt 1 ]
        then
            t=${t# }
            t=${t#\"}
            t=${t%\"}
            label_array+=("$t")
            t=""
            concat=""
        fi
    done
}

burnin=101

current_dir="$(pwd)"

project_dir="$(python ../project_util.py)"
ecoevolity_output_dir="${project_dir}/ecoevolity-output"
cd "$ecoevolity_output_dir"

plot_dir="../../results/stickleback-plots"
mkdir -p "$plot_dir"

labels='-l "bear" "Bear Paw Lake"
-l "root-bear" "Bear Paw Lake Ancestor"
-l "boot" "Boot Lake"
-l "root-boot" "Boot Lake Ancestor"
-l "mud" "Mud Lake"
-l "root-mud" "Mud Lake Ancestor"
-l "rabbit" "Rabbit Slough"
-l "root-rabbit" "Rabbit Slough Ancestor"
-l "resurrection" "Resurrection Bay"
-l "root-resurrection" "Resurrection Bay Ancestor"'

convert_labels_to_array $labels

time_ylabel="Population"
size_ylabel="Population"

upper_time=0.0225
upper_size=0.0075

bf_font_size=2.8

size_base_font=8.0

for datatype in "allsites" "snps"
do
    for dpp in "0_3725" "2_22543" "13_0"
    do
        for t in "1_0-0_001-0_0" "1_0-0_0005-0_0" "1_0-0_005-0_0"
        do
            for root in "1_0-1_0-0_0" "1_0-0_5-0_0" "1_0-0_1-0_0"
            do
                gzip -d -k run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log.gz
                pyco-sumtimes -f -y "$time_ylabel" -b $burnin "${label_array[@]}" -p "${plot_dir}/pyco-sumtimes-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
                pyco-sumtimes -f -x "" -y "" --x-limits 0.0 $upper_time -b $burnin "${label_array[@]}" -p "${plot_dir}/pyco-sumtimes-bare-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
                pyco-sumsizes -f --base-font-size $size_base_font -y "$size_ylabel" -b $burnin "${label_array[@]}" -p "${plot_dir}/pyco-sumsizes-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
                pyco-sumsizes -f --base-font-size $size_base_font -x "" -y "" --x-limits 0.0 $upper_size -b $burnin "${label_array[@]}" -p "${plot_dir}/pyco-sumsizes-bare-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
                sumcoevolity -b $burnin -n 1000000 -p "${plot_dir}/sumcoevolity-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" -c "../ecoevolity-configs/stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}.yml" run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
                pyco-sumevents -f --bf-font-size $bf_font_size -p "${plot_dir}/pyco-sumevents-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" --no-legend "${plot_dir}/sumcoevolity-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-sumcoevolity-results-nevents.txt"
                pyco-sumevents -f -x "" -y "" --bf-font-size $bf_font_size -p "${plot_dir}/pyco-sumevents-bare-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" --no-legend "${plot_dir}/sumcoevolity-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-sumcoevolity-results-nevents.txt"
                rm run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
            done
        done
    done
done

cd "$plot_dir"

for p in pyco-*.pdf
do
    pdfcrop "$p" "$p"
done

for p in grid-*.tex
do
    latexmk -pdf "$p"
done

for p in grid-*.pdf
do
    pdfcrop "$p" "$p"
done

cd "$current_dir"
