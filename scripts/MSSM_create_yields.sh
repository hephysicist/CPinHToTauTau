#!/bin/bash
source ./common_run3_MSSM.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --datasets $datasets
        --version $version
        --categories 'cat_emu_sr'
        --cf.CalibrateEvents-workflow $workflow
        --cf.SelectEvents-workflow $workflow
        --cf.ReduceEvents-workflow $workflow
        --cf.MergeReducedEvents-workflow $workflow
        "${@:2}"
    )
echo law run cf.CreateYieldTable "${args[@]}"
law run cf.CreateYieldTable "${args[@]}" 
