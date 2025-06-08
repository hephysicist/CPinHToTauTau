#!/bin/bash

set_common_vars() {


version="jet_veto_fix" #"njets_included" 
    
categories_mutau='cat_mutau_sr'
variables_mutau='mutau_lep1_pt,N_jets_pT_20_eta_4_7_Tight,mutau_mvis,mutau_delta_r,mutau_lep0_iso,mutau_lep0_eta,mutau_lep0_phi,mutau_lep0_ip_sig,mutau_lep0_pt,mutau_lep1_eta,mutau_lep1_phi,mutau_lep1_mass,mutau_lep1_decayModePNet,mutau_lep1_decayMode,mutau_mt,mutau_mvis,mutau_delta_r,mutau_pt,puppi_met_pt,puppi_met_phi,phi_cp_incl,phi_cp_mu_pi,phi_cp_mu_rho'

categories_etau="etau_signal_reg,etau_signal_reg_no_mt,etau_signal_reg_endcap_tau,etau_signal_reg_barrel_tau,etau_signal_reg_endcap_tau_no_mt,etau_signal_reg_barrel_tau_no_mt" 
#"etau_signal_reg_0_bjets,etau_signal_reg_0_bjets_endcap_tau,etau_signal_reg_0_bjets_barrel_tau,etau_signal_reg_0_bjets_no_mt,etau_signal_reg_0_bjets_endcap_tau_no_mt,etau_signal_reg_0_bjets_barrel_tau_no_mt,etau_signal_reg_1_bjets,etau_signal_reg_1_bjets_endcap_tau,etau_signal_reg_1_bjets_barrel_tau,etau_signal_reg_1_bjets_no_mt,etau_signal_reg_1_bjets_endcap_tau_no_mt,etau_signal_reg_1_bjets_barrel_tau_no_mt,etau_signal_reg_2_bjets,etau_signal_reg_2_bjets_endcap_tau,etau_signal_reg_2_bjets_barrel_tau,etau_signal_reg_2_bjets_no_mt,etau_signal_reg_2_bjets_endcap_tau_no_mt,etau_signal_reg_2_bjets_barrel_tau_no_mt"
variables_etau='mjj,N_jets_pT_30_eta_4_7_Tight,N_jets_pT_20_eta_2_5_Tight,Leading_jet_pt,Subleading_jet_pt,delta_eta_jj,jet_1_pt,etau_lep0_pt,etau_lep0_eta,etau_lep0_phi,etau_lep0_ip_sig,etau_lep1_pt,etau_lep1_eta,etau_lep1_phi,etau_lep1_mass,etau_lep1_decayModePNet,etau_lep1_decayMode,etau_mt,etau_mvis,etau_delta_r,etau_pt,puppi_met_pt,puppi_met_phi'

data_egamma_2022preEE='data_egamma_C,data_egamma_D,'
data_muoneg_2022preEE='data_muoneg_C,data_muoneg_D,'
data_mu_2022preEE='data_mu_C,data_mu_D,data_singlemu_C,'

data_egamma_2022postEE='data_egamma_E,data_egamma_F,data_egamma_G,'
data_muoneg_2022postEE='data_muoneg_E,data_muoneg_F,data_muoneg_G,'
data_mu_2022postEE='data_mu_E,data_mu_F,data_mu_G,'

bkg_ewk='wj_incl_madgraph,ww,wz,zz,dy_lep_madgraph,'
bkg_top='st_twchannel_t_sl,st_twchannel_t_dl,st_twchannel_tbar_sl,st_twchannel_tbar_dl,st_tchannel_tbar,st_tchannel_t,st_schannel_t_lep,st_schannel_tbar_lep,'
bkg_ttbar='tt_sl,tt_dl,tt_fh,'
signal='h_ggf_htt_cpo_filtered,h_ggf_htt_mm_filtered,h_ggf_htt_sm_filtered,'

data_egamma_2023preBPix='data_egamma_Cv123,data_egamma_Cv4,'
data_egamma_2023postBPix='data_egamma_D,'
data_muoneg_2023preBPix='data_muoneg_Cv123,data_muoneg_Cv4,'
data_muoneg_2023postBPix='data_muoneg_D,'
data_mu_2023preBPix='data_mu_Cv4,data_mu_Cv123,'
data_mu_2023postBPix='data_mu_D,'


case $1 in
##############################
####### 2022preEE ############
##############################
    "run3_2022preEE_etau_lim")
        config="run3_2022_preEE_etau_limited"	
        datasets='dy_lep_madgraph'
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee'
	    categories=$categories_etau
	    variables=$variables_etau
        workflow='local'
    ;;
    "run3_2022preEE_mutau_lim")
        config="run3_2022_preEE_mutau_limited"	
        datasets='wj_incl_madgraph'
        processes='wj'
	    categories=$categories_mutau
	    variables=$variables_mutau
        workflow='local'
    ;;
    "run3_2022preEE_etau")
        config="run3_2022_preEE_etau"
        data=$data_e_2022preEE
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar$signal"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,h_ggf_htt_cpo,h_ggf_htt_mm,h_ggf_htt_sm,data,'
	    categories=$categories_etau
	    variables=$variables_etau
        workflow='htcondor'
    ;;
    "run3_2022preEE_mutau")
        config="run3_2022_preEE_mutau"
        data=$data_mu_2022preEE
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar$signal"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,h_ggf_htt_cpo,h_ggf_htt_mm,h_ggf_htt_sm,data'
	    categories=$categories_mutau
	    variables=$variables_mutau
        workflow='htcondor'
     ;;
###############################
######## 2022postEE ###########
###############################    
    "run3_2022postEE_etau_lim")
        config="run3_2022_postEE_etau_limited"	
        datasets='data_e_E,dy_lep_madgraph' 
        processes='data,dy_z2tautau,dy_z2mumu,dy_z2ee' 
	    categories=$categories_etau
	    variables=$variables_etau
        workflow='local'
    ;;
    "run3_2022postEE_mutau_lim")
        config="run3_2022_postEE_mutau_limited"	
        datasets='wj_incl_madgraph' 
        processes='wj' 
        workflow='local'
    ;;
    "run3_2022postEE_etau")
        config="run3_2022_postEE_etau"
        data=$data_e_2022postEE
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets=$data$bkg_ewk$bkg_top$bkg_ttbar
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
	    categories=$categories_etau
	    variables=$variables_etau
	    workflow='htcondor'
    ;;
    "run3_2022postEE_mutau")
        config="run3_2022_postEE_mutau"
        data=$data_mu_2022postEE
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
      	datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
	    categories=$categories_mutau
	    variables=$variables_mutau
	    workflow='htcondor'
    ;;
##############################
####### 2023preBPix ############
##############################

    "run3_2023preBPix_etau")
        config="run3_2023_preBPix_etau"
        data=$data_e_2023preBPix
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
    categories=$categories_etau
    variables=$variables_etau
        workflow='htcondor'
    ;;
    "run3_2023preBPix_mutau")
        config="run3_2023_preBPix_mutau"
        data=$data_mu_2023preBPix
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
    categories=$categories_mutau
    variables=$variables_mutau
        workflow='htcondor'
     ;;

    "run3_2023preBPix_mutau_lim")
        config="run3_2023_preBPix_mutau_limited"
        datasets='dy_lep_madgraph'
        processes='dy_z2tautau'
    categories='mutau_signal_reg'
    variables=$variables_mutau
        workflow='local'
    ;;
##############################
####### 2023postBPix ############
##############################

    "run3_2023postBPix_etau")
        config="run3_2023_postBPix_etau"
        data=$data_e_2023postBPix
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
    categories=$categories_etau
    variables=$variables_etau
        workflow='htcondor'
    ;;
    "run3_2023postBPix_etau_lim")
        config="run3_2023_postBPix_etau_limited"
        data=$data_e_2023postBPix
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
    categories=$categories_etau
    variables=$variables_etau
        workflow='local'
    ;;

    "run3_2023postBPix_mutau")
        config="run3_2023_postBPix_mutau"
        data=$data_mu_2023postBPix
        bkg_ewk=$bkg_ewk
        bkg_top=$bkg_top
        bkg_ttbar=$bkg_ttbar
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_z2tautau,dy_z2mumu,dy_z2ee,vv,tt,st,wj,data'
    categories=$categories_mutau
    variables=$variables_mutau
        workflow='htcondor'
     ;;

    "run3_2023postBPix_mutau_lim")
        config="run3_2023_postBPix_mutau_limited"
        datasets='data_mu_D'
        processes='data'
    categories='mutau_signal_reg'
    variables=$variables_mutau
        workflow='local'
    ;;
    *)
    echo "Unknown run argument!"
    exit

esac
}
