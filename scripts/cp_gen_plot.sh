#!/bin/bash
source ./cp_samples_2025_v1.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --configs $configs
        --processes $processes
        --datasets  $datasets
        --categories cat_mutau_sr__fake_incl
        
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
        --cf.ProduceColumns-version add_gen_info
        --cf.CreateHistograms-workflow $workflow
        --cf.MergeHistograms-version add_gen_info
        --cf.MergeHistograms-workflow local
        --cf.PlotVariables1D-version add_gen_info
        
        --variables 'gen_lep0_pt,gen_lep0_pt,gen_lep0_phi,gen_lep0_IPx,gen_lep0_IPy,gen_lep0_IPz,'`
        `'gen_lep1_pt,gen_lep1_pt,gen_lep1_phi,gen_lep1_IPx,gen_lep1_IPy,gen_lep1_IPz,'

        --file-types pdf #,png
        --general-settings "cms-label=pw" #,yscale=log" # 
        --process-settings "h_ggf_htt_sm_prod_sm,unstack,scale=stack,color=#0000FF:h_ggf_htt_cpo_prod_sm,unstack,scale=stack,color=#FF0000:h_ggf_htt_mm_prod_sm,unstack,scale=stack,color=#00FF00"
        #"h_ggf_htt_cpo,unstack,scale=1,color=#FF0000:h_ggf_htt_sm,unstack,scale=30,color=#0000FF"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"