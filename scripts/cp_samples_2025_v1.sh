#!/bin/bash

set_common_vars() {

version="v1_1510"
    
categories='cat_mutau_sr,cat_mutau_sr__tau2a1_3pr,cat_mutau_sr__prompt__hig__cat2__tau2a1_3pr'
variables='hcand_mutau_fastMTT_BW_mass,hcand_mutau_fastMTT_BW_lep0_pt,hcand_mutau_fastMTT_BW_lep0_eta,hcand_mutau_fastMTT_BW_lep0_phi,hcand_mutau_fastMTT_BW_lep0_mass,hcand_mutau_fastMTT_BW_lep1_pt,hcand_mutau_fastMTT_BW_lep1_eta,hcand_mutau_fastMTT_BW_lep1_phi,hcand_mutau_fastMTT_BW_lep1_mass,theta_gj_mu_a1_3pr_pv_gef,theta_max_mu_a1_3pr_pv_gef,theta_rot_mu_a1_3pr_pv_gef,phi_cp_mu_a1_3pr_dp_reco,phi_cp_mu_a1_3pr_pv_reco,phi_cp_mu_a1_3pr_pv_mtt,phi_cp_mu_a1_3pr_pv_gef,bdt_raw_score_higgs,bdt_cat'

## all variables
# variables='npvs,
#            mutau_lep0_pt,mutau_lep0_eta,mutau_lep0_phi,mutau_lep0_mass,mutau_lep0_iso,
#            mutau_lep1_pt,mutau_lep1_eta,mutau_lep1_phi,mutau_lep1_mass,mutau_lep1_iso,
#            mutau_lep1_decayMode,mutau_lep1_decayModePNet,hcand_mutau_delta_r
#            mutau_mvis,mutau_mt,mutau_delta_r,puppi_met_pt,puppi_met_phi,
#            mutau_lep0_IPx,mutau_lep0_IPy,mutau_lep0_IPz,mutau_lep0_ip_sig,
#            mutau_lep1_IPx,mutau_lep1_IPy,mutau_lep1_IPz,mutau_lep1_ip_sig,
#            bdt_raw_score_gtau,bdt_raw_score_higgs,bdt_raw_score_fake,bdt_cat,
#            phi_cp_mu_a1_3pr_dp_reco,phi_cp_mu_a1_3pr_pv_reco,phi_cp_mu_a1_3pr_pv_mtt,phi_cp_mu_a1_3pr_pv_gef,
#            hcand_mutau_fastMTT_mass,hcand_mutau_fastMTT_lep0_pt,hcand_mutau_fastMTT_lep0_eta,hcand_mutau_fastMTT_lep0_phi,hcand_mutau_fastMTT_lep0_mass,hcand_mutau_fastMTT_lep1_pt,hcand_mutau_fastMTT_lep1_eta,hcand_mutau_fastMTT_lep1_phi,hcand_mutau_fastMTT_lep1_mass,
#            hcand_mutau_fastMTT_BW_mass,hcand_mutau_fastMTT_BW_lep0_pt,hcand_mutau_fastMTT_BW_lep0_eta,hcand_mutau_fastMTT_BW_lep0_phi,hcand_mutau_fastMTT_BW_lep0_mass,hcand_mutau_fastMTT_BW_lep1_pt,hcand_mutau_fastMTT_BW_lep1_eta,hcand_mutau_fastMTT_BW_lep1_phi,hcand_mutau_fastMTT_BW_lep1_mass,
#            hcand_mutau_fastMTT_cons_mass,hcand_mutau_fastMTT_cons_lep0_pt,hcand_mutau_fastMTT_cons_lep0_eta,hcand_mutau_fastMTT_cons_lep0_phi,hcand_mutau_fastMTT_cons_lep0_mass,hcand_mutau_fastMTT_cons_lep1_pt,hcand_mutau_fastMTT_cons_lep1_eta,hcand_mutau_fastMTT_cons_lep1_phi,hcand_mutau_fastMTT_cons_lep1_mass,
#            theta_gj_mu_a1_3pr_pv_gef,theta_max_mu_a1_3pr_pv_gef,theta_rot_mu_a1_3pr_pv_gef,
#            theta_gj_scaled_mu_a1_3pr_pv_gef,theta_max_scaled_mu_a1_3pr_pv_gef,theta_rot_scaled_mu_a1_3pr_pv_gef'
# --shape-norm for the last line

data_muoneg_2022preEE='data_muoneg_C,data_muoneg_D,'
data_mu_2022preEE='data_mu_C,data_mu_D,data_singlemu_C,'

data_egamma_2022postEE='data_egamma_E,data_egamma_F,data_egamma_G,'
data_muoneg_2022postEE='data_muoneg_E,data_muoneg_F,data_muoneg_G,'
data_mu_2022postEE='data_mu_E,data_mu_F,data_mu_G,'
bkg_dy='DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,'
bkg_wj='WtoLNu_madgraphMLM,WtoLNu_1J_madgraphMLM,WtoLNu_2J_madgraphMLM,WtoLNu_3J_madgraphMLM,WtoLNu_4J_madgraphMLM,'
bkg_vv='WW,WZ,ZZ,'
bkg_top='TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,'
bkg_ttbar='TTto2L2Nu,TTto4Q,TTtoLNu2Q,'
signal='h_ggf_htt_sm_prod_sm_filtered,h_ggf_htt_mm_prod_sm_filtered,h_ggf_htt_cpo_prod_sm_filtered,h_vbf_htt_sm_filtered,h_vbf_htt_mm_filtered,h_vbf_htt_cpo_filtered'

data_egamma_2023preBPix='data_egamma_Cv123,data_egamma_Cv4,'
data_egamma_2023postBPix='data_egamma_D,'
data_muoneg_2023preBPix='data_muoneg_Cv123,data_muoneg_Cv4,'
data_muoneg_2023postBPix='data_muoneg_D,'
data_mu_2023preBPix='data_mu_Cv123,data_mu_Cv4,'
data_mu_2023postBPix='data_mu_D,'


signal_procs='h_ggf_htt_sm_prod_sm,h_ggf_htt_mm_prod_sm,h_ggf_htt_cpo_prod_sm,h_vbf_htt_sm,h_vbf_htt_mm,h_vbf_htt_cpo'

case $1 in

##############################
####### Combined datasets ############
###############################
"22and23_mt")
        configs="run3_2022_preEE_mutau,run3_2022_postEE_mutau,run3_2023_preBPix_mutau,run3_2023_postBPix_mutau"
        bkgs="$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$signal"
        datasets="$data_mu_2022preEE${bkgs}:"
        datasets="${datasets}$data_mu_2022postEE${bkgs}:"
        datasets="${datasets}$data_mu_2023preBPix${bkgs}:"
        datasets="${datasets}$data_mu_2023postBPix${bkgs}"
        processes="$signal_procs",'dy_tt_m50,dy_ll_m50,vv,tt,st,wj,data'
        workflow='htcondor'
      ;;
"22and23_mt_sig")
        configs="run3_2022_preEE_mutau,run3_2022_postEE_mutau,run3_2023_preBPix_mutau,run3_2023_postBPix_mutau"
        datasets="$signal:$signal:$signal:$signal"
        processes="$signal_procs"
        workflow='htcondor'
      ;;
##############################
####### 2022preEE ############
##############################
"22pre_mt")
        config="run3_2022_preEE_mutau"
        bkgs="$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$signal"
        datasets="$data_mu_2022preEE${bkgs}"
        processes="$signal_procs",'dy_tt_m50,dy_ll_m50,vv,tt,st,wj,data'
        workflow='htcondor'
      ;; 
"22pre_mt_lim")
        config="run3_2022_preEE_mutau_limited"	
        datasets='WtoLNu_madgraphMLM'
        processes='wj'
        workflow='local'
      ;;
"22pre_mt_sig")
        configs="run3_2022_preEE_mutau"
        datasets="$signal"
        processes="$signal_procs"
        workflow='htcondor'
      ;;
###############################
######## 2022postEE ###########
###############################
"22post_mt")
        config="run3_2022_postEE_mutau"
        bkgs="$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$signal"
        datasets="$data_mu_2022postEE${bkgs}"
        processes="$signal_procs",'dy_tt_m50,dy_ll_m50,vv,tt,st,wj,data'
        workflow='htcondor'
      ;; 
"22post_mt_lim")
        config="run3_2022_postEE_mutau_limited"	
        datasets='WtoLNu_madgraphMLM' 
        processes='wj' 
        workflow='local'
      ;;
##############################
####### 2023preBPix ##########
##############################
"23pre_mt")
        config="run3_2023_preBPix_mutau"
        bkgs="$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$signal"
        datasets="$data_mu_2023preBPix${bkgs}"
        processes="$signal_procs",'dy_tt_m50,dy_ll_m50,vv,tt,st,wj,data'
        workflow='htcondor'
    ;;
"23pre_mt_lim")
        config="run3_2023_preBPix_mutau_limited"	
        datasets='WtoLNu_madgraphMLM' 
        processes='wj' 
        workflow='local'
    ;;
# ##############################
# ####### 2023postBPix #########
# ##############################
"23post_mt")
        config="run3_2023_postBPix_mutau"
        bkgs="$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$signal"
        datasets="$data_mu_2023postBPix${bkgs}"
        processes="$signal_procs",'dy_tt_m50,dy_ll_m50,vv,tt,st,wj,data'
        workflow='htcondor'
      ;;
"23post_mt_lim")
        config="run3_2023_postBPix_mutau_limited"
        datasets='WtoLNu_madgraphMLM'
        processes='wj'
        workflow='local'
    ;;
 
    *)
    echo "Unknown run argument!"
    exit

esac
}
