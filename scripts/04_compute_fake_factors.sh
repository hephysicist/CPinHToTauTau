#!/bin/bash
source ./common_run3.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
ff_version='debug_include_tes'
main_category='cat_mutau_sr' #Specify the category for which the fake factors should be calculated
args=(
        --config $config 
        --datasets $datasets
        --processes $processes
        --cf.CalibrateEvents-workflow $workflow
        --cf.CalibrateEvents-version $version
        
        --cf.SelectEvents-workflow $workflow
        --cf.SelectEvents-version $version

        --cf.ReduceEvents-workflow $workflow
        --cf.ReduceEvents-version $version
        
        --cf.MergeReducedEvents-workflow $workflow
        --cf.MergeReducedEvents-version $version

        --cf.ProvideReducedEvents-version $version

        --cf.ProduceColumns-producer 'ff_method'
        --cf.ProduceColumns-version  $ff_version

        --cf.PrepareFakeFactorHistograms-version $ff_version
        --cf.PrepareFakeFactorHistograms-categories $main_category

        --cf.MergeFakeFactorHistograms-version $ff_version

        --cf.ComputeFakeFactors-version  "met_recoil" 
        --cf.ComputeFakeFactors-categories $main_category
        
        "${@:2}"
    )
echo law run cf.ComputeFakeFactors "${args[@]}"
law run cf.ComputeFakeFactors "${args[@]}" 