#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --version $version
        --datasets $datasets
        --workflow local
        #--selector-steps "trigger,dilepton_veto,has_proper_tau_decay_products,b_veto,has_lep_pair_b4_trigmatch,has_lep_pair_after_trigmatch" 
        "${@:2}"
    )
echo run cf.PlotCutflow "${args[@]}"
law run cf.PlotCutflow "${args[@]}"