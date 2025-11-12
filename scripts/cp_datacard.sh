#!/bin/bash
source ./common_run3_MSSM.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
prod_version=bdt_check
args=(
        --config $config
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
        
        
        
        --cf.ProduceColumns-version $prod_version
        --cf.CreateHistograms-version $prod_version
        --cf.MergeHistograms-version $prod_version
        
        --version $prod_version
       
        --inference-model hcp_model
        --hist-hooks good_old_abcd
        "${@:2}"
    )
echo law run cf.CreateDatacards "${args[@]}"
law run cf.CreateDatacards "${args[@]}"
