#!/bin/bash
source ./common_run3_MSSM.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
wrapper_args=(
    --configs $config
    --datasets $datasets
    --version $version
    
   
    --cf.CalibrateEvents-workflow $workflow
    --cf.CalibrateEvents-htcondor-memory 4096 
    --cf.CalibrateEvents-max-runtime 2
    --cf.SelectEvents-workflow $workflow
    --cf.SelectEvents-htcondor-memory 4096
    --cf.SelectEvents-max-runtime 2
    --cf.ReduceEvents-workflow $workflow
    --cf.ReduceEvents-htcondor-memory 8192
    --cf.ReduceEvents-max-runtime 2
    --cf.MergeReducedEvents-workflow $workflow
    "${@:2}"
    )
echo law run cf.MergeReducedEventsWrapper "${wrapper_args[@]}"
law run cf.MergeReducedEventsWrapper "${wrapper_args[@]}"
