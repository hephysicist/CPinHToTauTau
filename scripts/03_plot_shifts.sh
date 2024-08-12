#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --datasets $datasets
        --version signal_test
        --categories mutau
        --cf.CalibrateEvents-workflow htcondor
        --cf.SelectEvents-workflow htcondor
        --cf.ReduceEvents-workflow htcondor
        --variables phi_cp_mu_pi,phi_cp_mu_rho
        --shift-sources tauspinner
        --process-settings "tauspinner_up,label=cp_odd:tauspinner_down,label=cp_even"
        --general-settings "cms-label=pw"
        "${@:2}"
    )
echo law run cf.PlotShiftedVariables1D "${args[@]}"
law run cf.PlotShiftedVariables1D "${args[@]}"
