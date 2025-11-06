#!/bin/bash

set_common_vars() {


version="v1_25"
    
categories_mutau='cat_mutau_sr'

data_muoneg_2022preEE='data_muoneg_C,data_muoneg_D,'
data_mu_2022preEE='data_mu_C,data_mu_D,data_singlemu_C,'

data_egamma_2022postEE='data_egamma_E,data_egamma_F,data_egamma_G,'
data_muoneg_2022postEE='data_muoneg_E,data_muoneg_F,data_muoneg_G,'
data_mu_2022postEE='data_mu_E,data_mu_F,data_mu_G,'
bkg_dy='DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,'
bkg_wj='WtoLNu_1J_madgraphMLM,WtoLNu_2J_madgraphMLM,WtoLNu_3J_madgraphMLM,WtoLNu_4J_madgraphMLM,WtoLNu_madgraphMLM,'
bkg_vv='WW,WZ,ZZ,'
#DYto2Tau_MLL_50_0J_Filtered_amcatnloFXFX,DYto2Tau_MLL_50_1J_Filtered_amcatnloFXFX,DYto2Tau_MLL_50_2J_Filtered_amcatnloFXFX
bkg_top='TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,'
bkg_ttbar='TTto2L2Nu,TTto4Q,TTtoLNu2Q,'
signal='h_ggf_htt_sm_prod_sm_filtered' #,h_ggf_htt_mm_prod_sm_filtered,h_ggf_htt_cpo_prod_sm_filtered,h_vbf_htt_sm_filtered,h_vbf_htt_mm_filtered,h_vbf_htt_cpo_filtered'

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
##############################
"22and23_mt_")
        configs="run3_2022_preEE_mutau,run3_2022_postEE_mutau,run3_2023_preBPix_mutau,run3_2023_postBPix_mutau"
        bkgs="$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar"
        datasets="$data_mu_2022preEE${bkgs}:"
        datasets="${datasets}$data_mu_2022postEE${bkgs}:"
        datasets="${datasets}$data_mu_2022preBPix${bkgs}:"
        datasets="${datasets}$data_mu_2022postBPix${bkgs}"
        processes='dy_tt_m50,dy_ll_m50,vv,tt,st,wj,data'
        workflow='htcondor'
     ;;
"22and23_mt")
        configs="run3_2022_preEE_mutau,run3_2022_postEE_mutau,run3_2023_preBPix_mutau,run3_2023_postBPix_mutau"
        bkgs="$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$signal"
        datasets="$data_mu_2022preEE${bkgs}:"
        datasets="${datasets}$data_mu_2022postEE${bkgs}:"
        datasets="${datasets}$data_mu_2022preBPix${bkgs}:"
        datasets="${datasets}$data_mu_2022postBPix${bkgs}"
        processes='dy_tt_m50,dy_ll_m50,vv,tt,st,wj,h_ggf_htt_sm_prod_sm,data'
        workflow='htcondor'
     ;;
"22post_and23pre_mt")
        configs="run3_2022_postEE_mutau,run3_2023_preBPix_mutau"
        datasets="$data_mu_2022postEE$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$:$data_mu_2023preBPix$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$"
        processes='dy_tt_m50,dy_ll_m50,vv,tt,st,wj,data'
        workflow='htcondor'
     ;;
##############################
####### 2022preEE ############
##############################
    "22pre_mt_lim")
        config="run3_2022_preEE_mutau_limited"	
        datasets='WtoLNu_madgraphMLM'
        processes='wj'
        workflow='local'
    ;;
    "22pre_mt")
        config="run3_2022_preEE_mutau"
        data=$data_mu_2022preEE
        datasets=$data$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$signal
        processes='dy_tt_m50,dy_ll_m50,vv,tt,st,wj,h_ggf_htt_sm_prod_sm,data'
        workflow='htcondor'
     ;;
###############################
######## 2022postEE ###########
###############################    
    "22post_mt_lim")
        config="run3_2022_postEE_mutau_limited"	
        datasets='WtoLNu_madgraphMLM' 
        processes='wj' 
        workflow='local'
    ;;
    "22post_mt")
        config="run3_2022_postEE_mutau"
        data=$data_mu_2022postEE
        signal='h_ggf_htt_sm_prod_sm_filtered,'
      	datasets=$data$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$signal
        processes='dy_tt_m50,dy_ll_m50,vv,tt,st,wj,h_ggf_htt_sm_prod_sm,data'
	    workflow='local'
    ;;
##############################
####### 2023preBPix ############
##############################
    "23pre_mt_lim")
        config="run3_2023_preBPix_mutau_limited"	
        datasets='WtoLNu_madgraphMLM' 
        processes='wj' 
        workflow='local'
    ;;
    "23pre_mt")
        config="run3_2023_preBPix_mutau"
        data=$data_mu_2023preBPix
        signal='h_ggf_htt_sm_prod_sm_filtered,'
      	datasets=$data$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$signal 
        processes='dy_tt_m50,dy_ll_m50,vv,tt,st,wj,h_ggf_htt_sm_prod_sm,data'
	    workflow='local'
    ;;
#     ;;
# ##############################
# ####### 2023postBPix ############
# ##############################
   "23post_mt_lim")
        config="run3_2023_postBPix_mutau_limited"
        datasets='WtoLNu_madgraphMLM'
        processes='wj'
        workflow='local'
    ;;

    "23post_mt")
        config="run3_2023_postBPix_mutau"
        data=$data_mu_2023postBPix
        datasets=$data$bkg_wj$bkg_dy$bkg_vv$bkg_top$bkg_ttbar$signal
        processes='dy_tt_m50,dy_ll_m50,vv,tt,st,wj,h_ggf_htt_sm_prod_sm,data'
        workflow='htcondor'
     ;;

 
    *)
    echo "Unknown run argument!"
    exit

esac
}
