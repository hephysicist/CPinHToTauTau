#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --version $version
        --categories mutau,etau,tautau
        --variables tau_1_pt,hcand_mass #,mutau_mass_no_tes,muon_pt,muon_eta,muon_phi,tau_pt,tau_eta,tau_phi,muon_mT,mutau_mass,met_pt,met_phi
        #--variables muon_iso #,muon_eta_1bin,mutau_mass_1bin
        --general-settings "cms-label=pw"
        "${@:2}"
    )
#echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"
