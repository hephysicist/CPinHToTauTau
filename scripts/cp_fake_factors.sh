#!/bin/bash
source ./common_run3.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --configs $configs
        --processes $processes
        --datasets $datasets
        
       

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
        
        --cf.MergeSelectionStats-version $version
        --cf.ProvideReducedEvents-version $version

        --cf.ProduceColumns-workflow $workflow
        --cf.ProduceColumns-version 'cf0p3_ff_test'
        --cf.CreateHistograms-workflow $workflow
       
        --cf.ComputeFakeFactors-version 'cf0p3_4bins'

        --categories  'cat_mutau_sr'
        --variables 'tau_pt-tau_dm_pnet-n_jets'
        --version cf0p3_4bins
        --hist-producer ff_hist_producer
        
        "${@:2}"
    )
echo law run cf.ComputeFakeFactors "${args[@]}"
law run cf.ComputeFakeFactors "${args[@]}" 
