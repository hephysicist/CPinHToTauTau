#!/bin/bash
source ./cp_samples_2025_v1.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --configs $config
        --processes $processes
        --datasets $datasets
        --categories   'cat_mutau_sr,cat_mutau_sr__hig__tau2a1_3pr' #'cat_mutau_ar_yields,cat_mutau_ar_yields_fakes,cat_mutau_dr_num_qcd,cat_mutau_dr_den_qcd,cat_mutau_dr_den_wj,cat_mutau_dr_num_wj,cat_mutau_ar_qcd,cat_mutau_ar_yields' #,cat_mutau_sr__tau2pi,cat_mutau_sr__tau2rho,cat_mutau_sr__tau2a1,cat_mutau_sr__tau2pi_3pr' #__tau2rho,cat_mutau_sr__tau2rho__hig_cat_0,cat_mutau_sr__tau2rho__hig_cat_1,cat_mutau_sr__tau2rho__hig_cat_2'
        #'cat_mutau_sr__tau2rho,cat_mutau_sr__tau2rho__hig_cat_0,cat_mutau_sr__tau2rho__hig_cat_1,cat_mutau_sr__tau2rho__hig_cat_2'
        
        --cf.CalibrateEvents-workflow $workflow
        --cf.CalibrateEvents-version $version

        --cf.SelectEvents-workflow $workflow
        --cf.SelectEvents-version $version
        
        --cf.ReduceEvents-workflow $workflow
        --cf.ReduceEvents-version $version

        --cf.MergeReducedEvents-workflow $workflow
        --cf.MergeReducedEvents-version $version

        
        --cf.MergeSelectionStats-version $version
        --cf.ProvideReducedEvents-version $version

        --cf.ProduceColumns-workflow $workflow
        --cf.CreateHistograms-workflow $workflow
        --cf.MergeHistograms-workflow local
        --version bugfree_samples_full_kine_MTT_GJ_3_nov
        --variables $variables
        
        #mutau_mt0, mutau_lep0_pt,mutau_lep1_pt,mutau_lep0_eta,mutau_lep1_eta,mutau_lep0_phi,mutau_lep1_phi,mutau_mt0,mutau_lep0_iso #bdt_raw_score_fake,bdt_raw_score_gtau,bdt_raw_score_higgs
        
        #'mutau_lep0_ipx,mutau_lep0_ipx_qm,mutau_lep1_ipx,mutau_lep1_ipx_qm,'`
        #`'mutau_lep0_ipy,mutau_lep0_ipy_qm,mutau_lep1_ipy,mutau_lep1_ipy_qm,'`
        #`'mutau_lep0_ipz,mutau_lep0_ipz_qm,mutau_lep1_ipz,mutau_lep1_ipz_qm'
        #'mutau_mvis,mutau_lep0_pt,mutau_lep1_pt,n_j,mutau_lep1_decayModePNet,mutau_lep0_eta,mutau_lep1_eta,puppi_met_pt' #,mutau_pt_vis,mutau_pt_vis,mutau_mt_ll,'`
        
        #`'mutau_mt0,mutau_mt1,mutau_mt_tot,mutau_delta_eta,mutau_delta_r,'`
        #`'leading_jet_pt,subleading_jet_pt,leading_jet_eta,subleading_jet_eta,leading_jet_phi,subleading_jet_phi,'`
        #`'dijet_delta_eta,mjj,N_jets_pT_20_eta_4_7_Tight,N_b_jets,'`
        #`'mutau_delta_phi,mutau_delta_phi_0_met,mutau_delta_phi_1_met,'`
        #'bdt_raw_score_gtau,bdt_raw_score_higgs,bdt_raw_score_fake,bdt_cat'
        #'phi_cp_mu_rho'
        #'mutau_lep1_pt,mutau_mvis,mutau_mvis_fine,mutau_lep0_ipx,mutau_lep0_ipx_qm,mutau_lep1_ipx,mutau_lep1_ipx_qm,'`
        #`'mutau_lep0_ipy,mutau_lep0_ipy_qm,mutau_lep1_ipy,mutau_lep1_ipy_qm,'`
        #`'mutau_lep0_ipz,mutau_lep0_ipz_qm,mutau_lep1_ipz,mutau_lep1_ipz_qm'
        --file-types pdf,png
       # Currently there are three methods that are implemented as hist hooks:
       # 1. ff_method: general fake factor method that requires ff weights to be present at the events tree
       # 2. ff_method_dr_closure_test: Calclulate fake contribution and apply it to the dr_num regions for the closure tests
       # 3. good_old_abcd: estimates QCD contribution by taking events from same sign region and transfer factors from inv. lep iso
         
       --hist-hooks qcd,incl,order #ff_method_dr_closure_test,ff_method,incl,
        #--skip-ratio
        --general-settings "cms-label=pw" #,yscale=log" # 
        --process-settings "h_ggf_htt_cpo_prod_sm,unstack,scale=stack,color=#28348e:h_ggf_htt_mm_prod_sm,unstack,scale=stack,color=#2b663c:h_ggf_htt_sm_prod_sm,unstack,scale=stack,color=#d62839:h_vbf_htt_cpo,unstack,scale=stack,color=#1f77b4:h_vbf_htt_mm,unstack,scale=stack,color=#ff7f0e:h_vbf_htt_sm,unstack,scale=stack,color=#2ca02c"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"
