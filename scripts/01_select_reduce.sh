#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
wrapper_args=(
    --configs $config
    --datasets $datasets
    --version $version
    --cf.CalibrateEvents-workflow $workflow
    --cf.SelectEvents-workflow $workflow
    --cf.ReduceEvents-workflow $workflow
    --cf.MergeReducedEvents-workflow local
    "${@:2}"
    )
echo law run cf.MergeReducedEventsWrapper "${wrapper_args[@]}"
law run cf.MergeReducedEventsWrapper "${wrapper_args[@]}"
