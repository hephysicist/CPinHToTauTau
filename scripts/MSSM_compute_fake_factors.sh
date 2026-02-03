#!/bin/bash
source ./common_run3_MSSM.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
ff_version='desy_dev'
main_category='cat_emu_sr' #Specify the category for which the fake factors should be calculated
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

        --cf.ProduceColumns-producer 'main'
        --cf.ProduceColumns-workflow $workflow
        --cf.ProduceColumns-version  $ff_version

        --cf.PrepareFakeFactorHistograms-version $ff_version
        --cf.PrepareFakeFactorHistograms-categories $main_category

        --cf.MergeFakeFactorHistograms-version $ff_version

        --cf.ComputeFakeFactors-version  "desy_dev" 
        --cf.ComputeFakeFactors-categories $main_category
        
        "${@:2}"
    )
echo law run cf.ComputeFakeFactors "${args[@]}"
law run cf.ComputeFakeFactors "${args[@]}" 