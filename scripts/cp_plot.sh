#!/bin/bash
source ./cp_samples_2025_v1.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --configs $configs
        --processes $processes
        --datasets  $datasets
        --categories  cat_mutau_sr__fake_incl__tau2a1_3pr
        #'cat_mutau_sr__fake_incl__tau2pi,cat_mutau_sr__fake_incl__tau2rho,'`
       #`'cat_mutau_sr__fake_incl__tau2a1,cat_mutau_sr__fake_incl__tau2a1_3pr'

        #'cat_mutau_sr__fake_incl,cat_mutau_sr__fake_incl__tau2pi,cat_mutau_sr__fake_incl__tau2rho,'`
        #`'cat_mutau_sr__fake_incl__tau2a1,cat_mutau_sr__fake_incl__tau2a1_3pr' 
        
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
        --cf.ProduceColumns-version merge_check
        --cf.CreateHistograms-workflow $workflow
        --cf.MergeHistograms-version merge_check
        --cf.MergeHistograms-workflow local
        --cf.PlotVariables1D-version merge_check
        
        --variables phi_cp_mu_a1_3pr_dp_reco,phi_cp_mu_a1_3pr_pv_gef,theta_gj_mu_a1_3pr_pv_gef,theta_max_mu_a1_3pr_pv_gef,theta_rot_mu_a1_3pr_pv_gef
        #'mutau_lep1_decayModePNet,mutau_mvis,mutau_lep0_pt,mutau_lep1_pt,mutau_pt_vis,puppi_met_pt,leading_jet_pt,subleading_jet_pt,n_j,'`
        #`'mutau_pt_vis_coarse,mutau_lep0_pt_coarse,mutau_lep1_pt_coarse,puppi_met_pt_coarse,mutau_mvis_coarse,mutau_mt0'
        #'mutau_lep1_pt,n_j,mutau_lep1_decayModePNet'
        --file-types pdf #,png
       # 1. ff_method: general fake factor method that requires ff weights to be present at the events tree
       # 2. closure_test: Calclulate fake contribution and apply it to the dr_num regions for the closure tests
       # 3. good_old_abcd: estimates QCD contribution by taking events from same sign region and transfer factors from inv. lep iso
         
       --hist-hooks qcd,order #ff_method_dr_closure_test,qcd,incl,ff_method
        #--skip-ratio
        --general-settings "cms-label=pw" #,yscale=log" # 
        --process-settings "h_ggf_htt_sm_prod_sm,unstack,scale=stack,color=#0000FF:h_ggf_htt_cpo_prod_sm,unstack,scale=stack,color=#FF0000:h_ggf_htt_mm_prod_sm,unstack,scale=stack,color=#00FF00"
        #"h_ggf_htt_cpo,unstack,scale=1,color=#FF0000:h_ggf_htt_sm,unstack,scale=30,color=#0000FF"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"