#!/bin/bash
source ./common_run3.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --datasets $datasets
        
        --categories  'cat_mutau_sr__tau2rho,cat_mutau_sr__tau2rho__hig_cat_0,cat_mutau_sr__tau2rho__hig_cat_1,cat_mutau_sr__tau2rho__hig_cat_2'
        
        #'cat_mutau_sr__tau2rho,cat_mutau_sr__tau2rho__hig_cat_0,cat_mutau_sr__tau2rho__hig_cat_1,cat_mutau_sr__tau2rho__hig_cat_2'
        
        #################################
        ### SIGNAL REGION  CATEGORIES ###
        #################################

        #'cat_mutau_sr' 
        #`',cat_mutau_sr__nj0_dm0,cat_mutau_sr__nj1_dm0,cat_mutau_sr__nj2_dm0,'`
        # #DM1
        # `'cat_mutau_sr__nj0_dm1,cat_mutau_sr__nj1_dm1,cat_mutau_sr__nj2_dm1,'`
        # #DM2
        # `'cat_mutau_sr__nj0_dm2,cat_mutau_sr__nj1_dm2,cat_mutau_sr__nj2_dm2,'`
        # #DM10
        # `'cat_mutau_sr__nj0_dm10,cat_mutau_sr__nj1_dm10,cat_mutau_sr__nj2_dm10,'`
        # #DM11
        # `'cat_mutau_sr__nj0_dm11,cat_mutau_sr__nj1_dm11,cat_mutau_sr__nj2_dm11,'


        #########################
        ### DR DEN CATEGORIES ###
        #########################

        #'cat_mutau_dr_den_qcd,cat_mutau_dr_den_wj,'
        # #dm1
        # `'cat_mutau_dr_den_qcd__nj0_dm0,cat_mutau_dr_den_wj__nj0_dm0,'`
        # `'cat_mutau_dr_den_qcd__nj1_dm0,cat_mutau_dr_den_wj__nj1_dm0,'`
        # `'cat_mutau_dr_den_qcd__nj2_dm0,cat_mutau_dr_den_wj__nj2_dm0,'`
        # #dm1
        # `'cat_mutau_dr_den_qcd__nj0_dm1,cat_mutau_dr_den_wj__nj0_dm1,'`
        # `'cat_mutau_dr_den_qcd__nj1_dm1,cat_mutau_dr_den_wj__nj1_dm1,'`
        # `'cat_mutau_dr_den_qcd__nj2_dm1,cat_mutau_dr_den_wj__nj2_dm1,'`
        #  #dm2
        # `'cat_mutau_dr_den_qcd__nj0_dm2,cat_mutau_dr_den_wj__nj0_dm2,'`
        # `'cat_mutau_dr_den_qcd__nj1_dm2,cat_mutau_dr_den_wj__nj1_dm2,'`
        # `'cat_mutau_dr_den_qcd__nj2_dm2,cat_mutau_dr_den_wj__nj2_dm2,'`
        #  #dm10
        # `'cat_mutau_dr_den_qcd__nj0_dm10,cat_mutau_dr_den_wj__nj0_dm10,'`
        # `'cat_mutau_dr_den_qcd__nj1_dm10,cat_mutau_dr_den_wj__nj1_dm10,'`
        # `'cat_mutau_dr_den_qcd__nj2_dm10,cat_mutau_dr_den_wj__nj2_dm10,'`
        #  #dm11
        # `'cat_mutau_dr_den_qcd__nj0_dm11,cat_mutau_dr_den_wj__nj0_dm11,'`
        # `'cat_mutau_dr_den_qcd__nj1_dm11,cat_mutau_dr_den_wj__nj1_dm11,'`
        # `'cat_mutau_dr_den_qcd__nj2_dm11,cat_mutau_dr_den_wj__nj2_dm11,'`

        # nj_inclusive cats
        # `'cat_mutau_dr_den_qcd__dm0,cat_mutau_dr_den_wj__dm0,'`
        # `'cat_mutau_dr_den_qcd__dm1,cat_mutau_dr_den_wj__dm1,'`
        # `'cat_mutau_dr_den_qcd__dm2,cat_mutau_dr_den_wj__dm2,'`
        # `'cat_mutau_dr_den_qcd__dm10,cat_mutau_dr_den_wj__dm10,'`
        # `'cat_mutau_dr_den_qcd__dm11,cat_mutau_dr_den_wj__dm11,'


        #########################
        ### DR NUM CATEGORIES ###
        #########################
        
        #'cat_mutau_dr_num_qcd,cat_mutau_dr_num_wj,'
        # `'cat_mutau_dr_num_qcd__nj0_dm0,cat_mutau_dr_num_wj__nj0_dm0,'`
        # `'cat_mutau_dr_num_qcd__nj1_dm0,cat_mutau_dr_num_wj__nj1_dm0,'`
        # `'cat_mutau_dr_num_qcd__nj2_dm0,cat_mutau_dr_num_wj__nj2_dm0,'`
        # #dm1
        # `'cat_mutau_dr_num_qcd__nj0_dm1,cat_mutau_dr_num_wj__nj0_dm1,'`
        # `'cat_mutau_dr_num_qcd__nj1_dm1,cat_mutau_dr_num_wj__nj1_dm1,'`
        # `'cat_mutau_dr_num_qcd__nj2_dm1,cat_mutau_dr_num_wj__nj2_dm1,'`
        # #  #dm2
        # `'cat_mutau_dr_num_qcd__nj0_dm2,cat_mutau_dr_num_wj__nj0_dm2,'`
        # `'cat_mutau_dr_num_qcd__nj1_dm2,cat_mutau_dr_num_wj__nj1_dm2,'`
        # `'cat_mutau_dr_num_qcd__nj2_dm2,cat_mutau_dr_num_wj__nj2_dm2,'`
        #  #dm10
        # `'cat_mutau_dr_num_qcd__nj0_dm10,cat_mutau_dr_num_wj__nj0_dm10,'`
        # `'cat_mutau_dr_num_qcd__nj1_dm10,cat_mutau_dr_num_wj__nj1_dm10,'`
        # `'cat_mutau_dr_num_qcd__nj2_dm10,cat_mutau_dr_num_wj__nj2_dm10,'`
        #  #dm11
        # `'cat_mutau_dr_num_qcd__nj0_dm11,cat_mutau_dr_num_wj__nj0_dm11,'`
        # `'cat_mutau_dr_num_qcd__nj1_dm11,cat_mutau_dr_num_wj__nj1_dm11,'`
        # `'cat_mutau_dr_num_qcd__nj2_dm11,cat_mutau_dr_num_wj__nj2_dm11,'`
        # nj_inclusive cats
        # 'cat_mutau_dr_num_qcd__dm0,cat_mutau_dr_num_wj__dm0,'`
        # `'cat_mutau_dr_num_qcd__dm1,cat_mutau_dr_num_wj__dm1,'`
        # `'cat_mutau_dr_num_qcd__dm2,cat_mutau_dr_num_wj__dm2,'`
        # `'cat_mutau_dr_num_qcd__dm10,cat_mutau_dr_num_wj__dm10,'`
        # `'cat_mutau_dr_num_qcd__dm11,cat_mutau_dr_num_wj__dm11,'

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
        --version $version
        
        --cf.MergeSelectionStats-version $version
        --cf.ProvideReducedEvents-version $version
        --version bdt_check
        --variables 
        #'mutau_mvis,mutau_pt_vis,mutau_pt_vis,mutau_mt_ll,'`
        #`'mutau_mt0,mutau_mt1,mutau_mt_tot,mutau_delta_eta,mutau_delta_r,'`
        #`'leading_jet_pt,subleading_jet_pt,leading_jet_eta,subleading_jet_eta,leading_jet_phi,subleading_jet_phi,'`
        #`'dijet_delta_eta,mjj,N_jets_pT_20_eta_4_7_Tight,N_b_jets,'`
        #`'mutau_delta_phi,mutau_delta_phi_0_met,mutau_delta_phi_1_met,'`
        #'bdt_raw_score_gtau,bdt_raw_score_higgs,bdt_raw_score_fake,bdt_cat'
        'phi_cp_mu_rho'
        #'mutau_lep1_pt,mutau_mvis,mutau_mvis_fine,mutau_lep0_ipx,mutau_lep0_ipx_qm,mutau_lep1_ipx,mutau_lep1_ipx_qm,'`
        #`'mutau_lep0_ipy,mutau_lep0_ipy_qm,mutau_lep1_ipy,mutau_lep1_ipy_qm,'`
        #`'mutau_lep0_ipz,mutau_lep0_ipz_qm,mutau_lep1_ipz,mutau_lep1_ipz_qm'
        --file-types pdf #,png
       # Currently there are three methods that are implemented as hist hooks:
       # 1. ff_method: general fake factor method that requires ff weights to be present at the events tree
       # 2. ff_method_dr_closure_test: Calclulate fake contribution and apply it to the dr_num regions for the closure tests
       # 3. good_old_abcd: estimates QCD contribution by taking events from same sign region and transfer factors from inv. lep iso
         
        #--hist-hooks good_old_abcd #ff_method_dr_closure_test
        --general-settings "yscale=log,cms-label=pw"
        --process-settings "h_ggf_htt_cpo,unstack,scale=1,color=#FF0000:h_ggf_htt_sm,unstack,scale=1,color=#0000FF"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"
