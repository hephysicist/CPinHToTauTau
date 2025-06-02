#!/bin/bash
source ./common_run3_MSSM.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --datasets $datasets
        --version $version
        --categories $categories
        --cf.CalibrateEvents-workflow $workflow
        --cf.SelectEvents-workflow $workflow
        --cf.ReduceEvents-workflow $workflow
        --cf.MergeReducedEvents-workflow $workflow
        --variables $variables
        --file-types pdf,png
	    --hist-hooks good_old_abcd
        --general-settings "cms-label=pw"
        --process-settings "h_ggf_htt,unstack,scale=stack"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"
 