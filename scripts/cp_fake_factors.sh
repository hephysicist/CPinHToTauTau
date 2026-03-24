#!/bin/bash
source ./cp_samples_2025_v1.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
prod_version=prod_march_19
hist_version=mar_23
args=(
        --configs $configs
        --processes $processes
        --datasets  $datasets
        --categories  cat_mutau_sr__fake_incl__bdt_incl
        
        --cf.CalibrateEvents-workflow $workflow
        --cf.CalibrateEvents-version $version

        --cf.SelectEvents-workflow $workflow
        --cf.SelectEvents-version $version
        
        --cf.ReduceEvents-workflow $workflow
        --cf.ReduceEvents-version $version
        --cf.ReduceEvents-tasks-per-job 3 

        --cf.MergeReducedEvents-workflow $workflow
        --cf.MergeReducedEvents-version $version

        
        --cf.MergeSelectionStats-version $version
        --cf.ProvideReducedEvents-version $version

        --cf.ProduceColumns-workflow $workflow
        --cf.ProduceColumns-version  $prod_version
        --cf.CreateHistograms-workflow $workflow
        --cf.CreateHistograms-version $hist_version
        --cf.MergeHistograms-version $hist_version
        --cf.MergeHistograms-workflow local

        --variables 'tau_pt-tau_dm_pnet-n_j'
        --version $hist_version
        --hist-producer ff_hist_producer
        
        "${@:2}"
    )
echo law run cf.ComputeFakeFactors "${args[@]}"
law run cf.ComputeFakeFactors "${args[@]}" 
