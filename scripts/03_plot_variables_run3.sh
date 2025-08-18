#!/bin/bash
source ./common_run3.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
prod_version=ic_sync_a1
args=(
        --config $config
        --processes 'dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data' #no signals 
        #--processes 'dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,h_ggf_htt_sm_prod_sm,h_ggf_htt_mm_prod_sm,h_ggf_htt_cpo_prod_sm,data' #sm, cpo and mm
        #--processes 'dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,h_ggf_htt_sm_prod_sm,h_vbf_htt_sm,zh_htt_sm,data' #with only sm signals 
        --datasets $datasets
        
        --categories 
        cat_mutau_sr
        #'cat_mutau_sr__hig__cat0__tau2rho,cat_mutau_sr__hig__cat1__tau2rho,cat_mutau_sr__hig__cat2__tau2rho,'`
        #'cat_mutau_sr__hig__cat0__tau2a1_3pr,cat_mutau_sr__hig__cat1__tau2a1_3pr,cat_mutau_sr__hig__cat2__tau2a1_3pr'
        #'cat_mutau_sr__hig__cat0__tau2a1,cat_mutau_sr__hig__cat1__tau2a1,cat_mutau_sr__hig__cat2__tau2a1'
        #'cat_mutau_sr__hig__cat0__tau2pi,cat_mutau_sr__hig__cat1__tau2pi,cat_mutau_sr__hig__cat2__tau2pi'
        
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
        
        --cf.ProduceColumns-version $prod_version
        --cf.CreateHistograms-version $prod_version
        --cf.MergeHistograms-version $prod_version
        
        --version $prod_version
        --variables mutau_fastMTT_mass #'mutau_mvis,mutau_fastMTT_mass,mutau_lep1_pt,mutau_lep0_pt,'`
        #'mutau_lep0_ipx_qm,mutau_lep0_ipx,mutau_lep1_ipx,mutau_lep1_ipx_qm,'`
        #`'mutau_lep0_ipy_qm,mutau_lep0_ipy,mutau_lep1_ipy,mutau_lep1_ipy_qm,'`
        #`'mutau_lep0_ipz_qm,mutau_lep0_ipz,mutau_lep1_ipz,mutau_lep1_ipz_qm,'

        #`'mutau_mt0,mutau_mt1,mutau_mt_tot,mutau_delta_eta,mutau_delta_r,'`
        #`'leading_jet_pt,subleading_jet_pt,leading_jet_eta,subleading_jet_eta,leading_jet_phi,subleading_jet_phi,'`
        #`'dijet_delta_eta,mjj,N_jets_pT_20_eta_4_7_Tight,N_b_jets,'`
        #`'mutau_delta_phi,mutau_delta_phi_0_met,mutau_delta_phi_1_met,'`
        #'bdt_raw_score_gtau,bdt_raw_score_higgs,bdt_raw_score_fake'
        #'phi_cp_mu_a1_3pr'
        #'phi_cp_mu_a1_1pr'
        #phi_cp_mu_pi
        #'mutau_lep1_pt,mutau_mvis,mutau_mvis_fine,mutau_lep0_ipx,mutau_lep0_ipx_qm,mutau_lep1_ipx,mutau_lep1_ipx_qm,'`
        #`'mutau_lep0_ipy,mutau_lep0_ipy_qm,mutau_lep1_ipy,mutau_lep1_ipy_qm,'`
        #`'mutau_lep0_ipz,mutau_lep0_ipz_qm,mutau_lep1_ipz,mutau_lep1_ipz_qm'
        #--file-types pdf #,png
    
        --hist-hooks good_old_abcd
        #--hist-hooks symmetrize_signal,flatten_dy,good_old_abcd,blind_sr #ff_method_dr_closure_test
        --general-settings "cms-label=pw" #,yscale=log"
        #--process-settings 'h_ggf_htt_sm_prod_sm,unstack,scale=stack,color=#0000FF:h_vbf_htt_sm,unstack,scale=stack,color=#00FFFF:zh_htt_sm,unstack,scale=stack,color=#FF00FF'
        #--process-settings 'h_ggf_htt_cpo_prod_sm,unstack,scale=20,color=#FF0000:h_ggf_htt_sm_prod_sm,unstack,scale=20,color=#0000FF:h_ggf_htt_mm_prod_sm,unstack,scale=20,color=#00FF00'
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"
