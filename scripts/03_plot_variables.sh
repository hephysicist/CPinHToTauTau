#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --version $version
        --categories mutau,etau,tautau
        --variables muon_1_pt,muon_1_eta,muon_1_phi,TauSpinner_weight_cp_0,TauSpinner_weight_cp_0p5
        --general-settings "cms-label=pw"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"
