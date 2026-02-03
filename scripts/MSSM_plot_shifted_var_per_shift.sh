#!/bin/bash
source ./common_run3_MSSM.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        # --configs $config
        --config $config
        --processes $processes
        --datasets $datasets
        --version $version
        --categories $categories
        --variables $variables
        --hist-hooks qcd
        --shift-sources pu_weight,top_pt_weight,Trigger_SF_weight,electron_weight,muon_weight,zpt_weight
        --file-types png
        --general-settings "cms-label=pw,yscale=log" #yscale=log,
        "${@:2}"
    )
echo run cf.PlotShiftedVariablesPerShift1D "${args[@]}"
law run cf.PlotShiftedVariablesPerShift1D "${args[@]}"