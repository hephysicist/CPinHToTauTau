#!/bin/bash
source ./commcommon_run3_MSSM.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --datasets $datasets
        --version $version
        --categories mutau #'mutau_signal_reg,tautau_signal_reg'
        --cf.CalibrateEvents-workflow $workflow
        --cf.SelectEvents-workflow $workflow
        --cf.ReduceEvents-workflow $workflow
        --cf.MergeReducedEvents-workflow local
        --variables  puppi_met_pt,jets_pt,jet_1_pt #hcand_leg1_pt,hcand_leg2_pt #phi_cp_mu_pi,phi_cp_mu_rho,phi_cp_mu_a1_1pr,phi_cp_rho_rho
        --shift-sources jer #tauspinner
        --general-settings "cms-label=simpw"
        --cf.PlotShiftedVariables1D-hide-errors True
        --hide-errors True
        "${@:2}"
    )
echo law run cf.PlotShiftedVariables1D "${args[@]}"
law run cf.PlotShiftedVariables1D "${args[@]}"