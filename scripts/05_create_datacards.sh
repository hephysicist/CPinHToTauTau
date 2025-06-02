#!/bin/bash
source ./common_run3_MSSM.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config

        --cf.CalibrateEvents-workflow $workflow
        --cf.CalibrateEvents-version $version
        --cf.SelectEvents-workflow $workflow
        --cf.SelectEvents-version $version
        --cf.MergeReducedEvents-workflow $workflow
        --cf.MergeReducedEvents-version $version
        --version $version
        --cf.MergeSelectionStats-version $version
        --cf.ProvideReducedEvents-version $version
        --inference-model example
        --hist-hooks good_old_abcd
        "${@:2}"
    )
echo law run cf.CreateDatacards "${args[@]}"
law run cf.CreateDatacards "${args[@]}"