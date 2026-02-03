#!/bin/bash
source ./common_run3_MSSM.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        # --configs $config
        --config $config
        --processes $processes
        --datasets $datasets
        --version $version
        --categories $categories
        --cf.CalibrateEvents-workflow $workflow
        --cf.SelectEvents-workflow $workflow
        --cf.ReduceEvents-workflow $workflow
        --cf.MergeReducedEvents-workflow $workflow
        --cf.ProduceColumns-workflow $workflow
        --cf.CreateHistograms-workflow $workflow
        --cf.MergeHistograms-workflow local
        --variables $variables
        --pilot True
        --file-types pdf,png
	--hist-hooks qcd
        --general-settings "cms-label=pw" #yscale=log,
        --process-settings "dy_lep,color=#FFFF00:h_ggf_htt_100,unstack,scale=100,color=#FF0000:bbh_htt_100,unstack,scale=10000,color=#0000FF"
        #"h_ggf_htt_80,unstack,scale=stack,color=#FF0000:h_ggf_htt_100,unstack,scale=stack,color=#0000FF:h_ggf_htt_120,unstack,scale=stack"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"