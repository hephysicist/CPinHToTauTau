#!/bin/bash
source ./cp_samples_2025_v1.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --configs $configs
        --processes $processes
        --datasets  $datasets
        --categories 'cat_mutau_sr__fake_incl__hig__cat_incl__tau2pi,cat_mutau_sr__fake_incl__hig__cat_incl__tau2rho,cat_mutau_sr__fake_incl__hig__cat_incl__tau2a1,cat_mutau_sr__fake_incl__hig__cat_incl__tau2a1_3pr'
        
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
        --cf.ProduceColumns-version gen_check8
        --cf.CreateHistograms-workflow $workflow
        --cf.CreateHistograms-version gen_check8
        --cf.MergeHistograms-version gen_check8
        --cf.MergeHistograms-workflow local
        --cf.PlotVariables1D-version gen_check8
        
        --variables 'phi_cp_mu_pi,phi_cp_mu_pi_gen,'`
        `'phi_cp_mu_rho,phi_cp_mu_rho_gen,'`
        `'phi_cp_mu_a1_1pr,phi_cp_mu_a1_1pr_gen,'`
        `'phi_cp_mu_a1_3pr_dp,phi_cp_mu_a1_3pr_dp_gen,'`
        `'phi_cp_mu_a1_3pr_pv_gef,phi_cp_mu_a1_3pr_pv_gen,phi_cp_mu_a1_3pr_pv_mtt,phi_cp_mu_a1_3pr_pv' 
        
        --file-types pdf #,png
        --general-settings "cms-label=pw" #,yscale=log" # 
        #--hist-hooks incl
        --skip-ratio
        --process-settings 'h_ggf_htt_sm_prod_sm,unstack,scale=2,color=#0000FF:h_ggf_htt_cpo_prod_sm,unstack,scale=2,color=#FF0000:h_ggf_htt_mm_prod_sm,unstack,scale=2,color=#00FF00:'`
        `'h_vbf_htt_sm,unstack,scale=10,color=#023047:h_vbf_htt_cpo,unstack,scale=10,color=#FF6F61:h_vbf_htt_mm,unstack,scale=10,color=#2EC4B6'
        
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"