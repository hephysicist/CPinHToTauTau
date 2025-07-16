#!/bin/bash
source ./common_run3.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes data
        --datasets  'data_mu_E,data_mu_F,data_mu_G'
        --categories  'cat_mutau_jet_veto_check_sr'

        --cf.CalibrateEvents-workflow $workflow
        --cf.CalibrateEvents-version $version
        
        --cf.SelectEvents-workflow $workflow
        --cf.SelectEvents-version $version
        
	--cf.ReduceEvents-workflow $workflow
	--cf.ReduceEvents-version $version

        --cf.MergeReducedEvents-workflow $workflow
        --cf.MergeReducedEvents-version $version
        
        --cf.MergeSelectionStats-version $version
        --cf.ProvideReducedEvents-version $version
        --version debug_include_tes 
        --variables jets_eta-jets_phi
        --file-types pdf,png
        --general-settings "cms-label=sim"
        "${@:2}"
    )
echo law run cf.PlotVariables2D "${args[@]}"
law run cf.PlotVariables2D "${args[@]}"
