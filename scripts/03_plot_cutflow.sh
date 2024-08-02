#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --version $version
        --datasets $datasets
        --workflow htcondor
        --selector-steps "trigger,met_filter,b_veto,trigobj_prematch,trigobj_postmatch,dilepton_veto,extra_lepton_veto,single_hcand,has_proper_tau_decay_products"
        "${@:2}"
    )
echo run cf.PlotCutflow "${args[@]}"
law run cf.PlotCutflow "${args[@]}"