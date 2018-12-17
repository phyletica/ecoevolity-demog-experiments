#!/bin/bash

function run_summary_tools () {
    upper_time=0.006
    upper_size=0.07
    # if [ "$datatype" = "snps" ]
    # then
    #     upper_time=0.000225
    #     upper_size=0.012
    # fi
    lower_time="$(echo "scale=8;-${upper_time}*0.02" | bc)"
    lower_size="$(echo "scale=8;-${upper_size}*0.02" | bc)"
    echo "time: ${lower_time}-${upper_time}"
    echo "size: ${lower_size}-${upper_size}"
    gzip -d -k run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log.gz
    pyco-sumtimes -f -w "$plot_width" --violin -y "$time_ylabel" -b $burnin "${label_array[@]}" -p "${plot_dir}/pyco-sumtimes-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
    pyco-sumtimes -f -w "$plot_width" --violin -x "" -y "" --x-limits $lower_time $upper_time -b $burnin "${label_array[@]}" -p "${plot_dir}/pyco-sumtimes-bare-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
    pyco-sumsizes -f -w "$plot_width" --violin --base-font-size $size_base_font -y "$size_ylabel" -b $burnin "${label_array[@]}" -p "${plot_dir}/pyco-sumsizes-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
    pyco-sumsizes -f -w "$plot_width" --violin --base-font-size $size_base_font -x "" -y "" --x-limits $lower_size $upper_size -b $burnin "${label_array[@]}" -p "${plot_dir}/pyco-sumsizes-bare-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
    if [ -e "${plot_dir}/sumcoevolity-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-sumcoevolity-results-nevents.txt" ]
    then
        rm "${plot_dir}/sumcoevolity-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-sumcoevolity-results-nevents.txt"
    fi
    if [ -e "${plot_dir}/sumcoevolity-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-sumcoevolity-results-model.txt" ]
    then
        rm "${plot_dir}/sumcoevolity-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-sumcoevolity-results-model.txt"
    fi
    sumcoevolity -b $burnin -n 1000000 -p "${plot_dir}/sumcoevolity-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" -c "${config_dir}/stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}.yml" run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
    pyco-sumevents -f -w "$plot_width" --bf-font-size $bf_font_size -p "${plot_dir}/pyco-sumevents-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" --no-legend "${plot_dir}/sumcoevolity-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-sumcoevolity-results-nevents.txt"
    pyco-sumevents -f -w "$plot_width" -x "" -y "" --bf-font-size $bf_font_size -p "${plot_dir}/pyco-sumevents-bare-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-" --no-legend "${plot_dir}/sumcoevolity-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-sumcoevolity-results-nevents.txt"
    rm run-?-stickleback-${datatype}-dpp-${dpp}-root-${root}-time-${t}-state-run-1.log
}

set -e

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
ecoevolity_output_dir="${project_dir}/stickleback-ecoevolity-output"
config_dir="${project_dir}/configs"

plot_dir="${project_dir}/results/stickleback-plots"
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

bf_font_size=4.0

size_base_font=9.0

plot_width=6.0

cd "$ecoevolity_output_dir"

for datatype in "allsites" "snps"
do
    for dpp in "2_22543"
    do
        for t in "1_0-0_001-0_0"
        do
            for root in "1_0-1_0-0_0" "1_0-0_5-0_0" "1_0-0_1-0_0"
            do
                run_summary_tools
            done
        done
    done
done

for datatype in "allsites" "snps"
do
    for dpp in "2_22543"
    do
        for t in "1_0-0_0005-0_0" "1_0-0_005-0_0"
        do
            for root in "1_0-1_0-0_0"
            do
                run_summary_tools
            done
        done
    done
done

for datatype in "allsites" "snps"
do
    for dpp in "0_3725" "13_0"
    do
        for t in "1_0-0_001-0_0"
        do
            for root in "1_0-1_0-0_0"
            do
                run_summary_tools
            done
        done
    done
done


cd "$plot_dir"

for p in pyco-*.pdf
do
    pdfcrop "$p" "$p"
done

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
