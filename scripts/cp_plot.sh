#!/bin/bash
source ./cp_samples_2025_v1.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --configs $configs
        --processes $processes
        --datasets  $datasets
        --categories  cat_mutau_sr__fake_incl__bdt_incl__tau2a1_3pr
        # 'cat_mutau_sr__fake_incl__hig__cat2__tau2a1_3pr'
        #`'cat_mutau_sr__fake_incl__tau2pi,cat_mutau_sr__fake_incl__tau2rho,'`
        #`'cat_mutau_sr__fake_incl__tau2a1,cat_mutau_sr__fake_incl__tau2a1_3pr'`
        # 'cat_mutau_sr__fake_incl,cat_mutau_sr__fake_incl__tau2pi,cat_mutau_sr__fake_incl__tau2rho,'`
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
        ## all variables
        # 'npvs,
        #  mutau_lep0_pt,mutau_lep0_eta,mutau_lep0_phi,mutau_lep0_mass,mutau_lep0_iso,
        #  mutau_lep1_pt,mutau_lep1_eta,mutau_lep1_phi,mutau_lep1_mass,mutau_lep1_iso,
        #  mutau_lep1_decayMode,mutau_lep1_decayModePNet,hcand_mutau_delta_r,
        #  mutau_pt_vis,leading_jet_pt,subleading_jet_pt,n_j,mutau_pt_vis_coarse,mutau_lep0_pt_coarse,mutau_lep1_pt_coarse,puppi_met_pt_coarse,mutau_mvis_coarse,mutau_mt0,
        #  mutau_mvis,mutau_mt,mutau_delta_r,puppi_met_pt,puppi_met_phi,
        #  mutau_lep0_IPx,mutau_lep0_IPy,mutau_lep0_IPz,mutau_lep0_ip_sig,
        #  mutau_lep1_IPx,mutau_lep1_IPy,mutau_lep1_IPz,mutau_lep1_ip_sig,
        #  bdt_raw_score_gtau,bdt_raw_score_higgs,bdt_raw_score_fake,bdt_cat,
        #  phi_cp_mu_a1_3pr_dp_reco,phi_cp_mu_a1_3pr_pv_reco,phi_cp_mu_a1_3pr_pv_mtt,phi_cp_mu_a1_3pr_pv_gef,
        #  hcand_mutau_fastMTT_mass,hcand_mutau_fastMTT_lep0_pt,hcand_mutau_fastMTT_lep0_eta,hcand_mutau_fastMTT_lep0_phi,hcand_mutau_fastMTT_lep0_mass,hcand_mutau_fastMTT_lep1_pt,hcand_mutau_fastMTT_lep1_eta,hcand_mutau_fastMTT_lep1_phi,hcand_mutau_fastMTT_lep1_mass,
        #  hcand_mutau_fastMTT_BW_mass,hcand_mutau_fastMTT_BW_lep0_pt,hcand_mutau_fastMTT_BW_lep0_eta,hcand_mutau_fastMTT_BW_lep0_phi,hcand_mutau_fastMTT_BW_lep0_mass,hcand_mutau_fastMTT_BW_lep1_pt,hcand_mutau_fastMTT_BW_lep1_eta,hcand_mutau_fastMTT_BW_lep1_phi,hcand_mutau_fastMTT_BW_lep1_mass,
        #  hcand_mutau_fastMTT_cons_mass,hcand_mutau_fastMTT_cons_lep0_pt,hcand_mutau_fastMTT_cons_lep0_eta,hcand_mutau_fastMTT_cons_lep0_phi,hcand_mutau_fastMTT_cons_lep0_mass,hcand_mutau_fastMTT_cons_lep1_pt,hcand_mutau_fastMTT_cons_lep1_eta,hcand_mutau_fastMTT_cons_lep1_phi,hcand_mutau_fastMTT_cons_lep1_mass,
        #  theta_gj_mu_a1_3pr_pv_gef,theta_max_mu_a1_3pr_pv_gef,theta_rot_mu_a1_3pr_pv_gef,
        #  theta_gj_scaled_mu_a1_3pr_pv_gef,theta_max_scaled_mu_a1_3pr_pv_gef,theta_rot_scaled_mu_a1_3pr_pv_gef'
        
        --file-types pdf,png

        --hist-hooks qcd,order #ff_method_dr_closure_test,qcd,incl,ff_method
        # 1. ff_method: general fake factor method that requires ff weights to be present at the events tree
        # 2. closure_test: Calclulate fake contribution and apply it to the dr_num regions for the closure tests
        # 3. good_old_abcd: estimates QCD contribution by taking events from same sign region and transfer factors from inv. lep iso

        #--skip-ratio
        --general-settings "cms-label=pw" #,yscale=log" # 
        --process-settings "h_ggf_htt_cpo_prod_sm,unstack,scale=stack,color=#28348e:h_ggf_htt_mm_prod_sm,unstack,scale=stack,color=#2b663c:h_ggf_htt_sm_prod_sm,unstack,scale=stack,color=#d62839:h_vbf_htt_cpo,unstack,scale=stack,color=#1f77b4:h_vbf_htt_sm,unstack,scale=stack,color=#ff7f0e:h_vbf_htt_mm,unstack,scale=stack,color=#2ca02c"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"