#!/bin/bash

set_common_vars() {

version="desy_dev"

variables_emu_list=( 
  "emu_mt_tot" "emu_mt_emu" "D_zeta" "D_zeta_check"
  "emu_mt_e" "emu_mt_mu" "N_jets_pT_20_eta_4_7_Tight"
  "leading_jet_eta" "subleading_jet_eta" "leading_jet_phi"
  "subleading_jet_phi" "N_b_jets" 
  "leading_jet_pt" "subleading_jet_pt" "dijet_delta_eta" "mjj" 
  "leading_b_jet_eta" "subleading_b_jet_eta" "leading_b_jet_phi"
  "subleading_b_jet_phi" "leading_b_jet_pt" "subleading_b_jet_pt" 
  "di_b_jet_delta_eta" "mb_jb_j" 
  "emu_lep0_pt" "emu_lep0_eta" "emu_lep0_phi"
  "emu_lep0_ip_sig" "emu_lep1_pt" "emu_lep1_eta" "emu_lep1_phi" "emu_lep1_ip_sig"
  "emu_mvis" "emu_delta_r" "emu_pt" "puppi_met_pt" "puppi_met_phi"
  "pt_H" "hcand_emu_fastMTT_mass" 
)
#   "bdt_raw_score_ggh_M100" "bdt_raw_score_bbh_M100" "bdt_raw_score_dy_M100" 
#   "bdt_raw_score_tt_M100" "bdt_raw_score_wj_M100" "bdt_raw_score_st_M100"
#   "bdt_raw_score_mb_M100"
# mssm_ggF_signal_list=(
#   "h_ggf_htt_60" "h_ggf_htt_65" "h_ggf_htt_70" "h_ggf_htt_75"
#   "h_ggf_htt_80" "h_ggf_htt_85" "h_ggf_htt_90" "h_ggf_htt_95"
#   "h_ggf_htt_100" "h_ggf_htt_105" "h_ggf_htt_110" "h_ggf_htt_115"
#   "h_ggf_htt_120" "h_ggf_htt_125" "h_ggf_htt_130" "h_ggf_htt_135"
#   "h_ggf_htt_140" "h_ggf_htt_160" "h_ggf_htt_180" "h_ggf_htt_200"
#   "h_ggf_htt_250" "h_ggf_htt_300" "h_ggf_htt_350" "h_ggf_htt_400"
#   "h_ggf_htt_450" "h_ggf_htt_500" "h_ggf_htt_600" "h_ggf_htt_700"
#   "h_ggf_htt_800" "h_ggf_htt_900" "h_ggf_htt_1000" "h_ggf_htt_1100"
#   "h_ggf_htt_1200" "h_ggf_htt_1400" "h_ggf_htt_1600" "h_ggf_htt_1800"
#   "h_ggf_htt_2000" "h_ggf_htt_2300" "h_ggf_htt_2600" "h_ggf_htt_2900"
#   "h_ggf_htt_3200" "h_ggf_htt_3500"
# )
# mssm_bbh_signal_list=(
#   "bbh_htt_60" "bbh_htt_65" "bbh_htt_70" "bbh_htt_75" "bbh_htt_80" "bbh_htt_85"
#   "bbh_htt_90" "bbh_htt_95" "bbh_htt_100" "bbh_htt_105" "bbh_htt_110" "bbh_htt_115"
#   "bbh_htt_120" "bbh_htt_125" "bbh_htt_130" "bbh_htt_135" "bbh_htt_140" "bbh_htt_160"
#   "bbh_htt_180" "bbh_htt_200" "bbh_htt_250" "bbh_htt_300" "bbh_htt_350" "bbh_htt_400"
#   "bbh_htt_450" "bbh_htt_500" "bbh_htt_600" "bbh_htt_700" "bbh_htt_800" "bbh_htt_900"
#   "bbh_htt_1000" "bbh_htt_1100" "bbh_htt_1200" "bbh_htt_1400" "bbh_htt_1600"
#   "bbh_htt_1800" "bbh_htt_2000" "bbh_htt_2300" "bbh_htt_2600" "bbh_htt_2900"
#   "bbh_htt_3200" "bbh_htt_3500"
# )

variables_emu=$(IFS=,; echo "${variables_emu_list[*]}")


data_egamma_2022preEE='data_egamma_C,data_egamma_D,'
data_muoneg_2022preEE='data_muoneg_C,data_muoneg_D,'
data_mu_2022preEE='data_mu_C,data_mu_D,data_singlemu_C,'

data_egamma_2022postEE='data_egamma_E,data_egamma_F,data_egamma_G,'
data_muoneg_2022postEE='data_muoneg_E,data_muoneg_F,data_muoneg_G,'
data_mu_2022postEE='data_mu_E,data_mu_F,data_mu_G,'

data_egamma_2023preBPix='data_egamma_Cv123,data_egamma_Cv4,'
data_muoneg_2023preBPix='data_muoneg_Cv123,data_muoneg_Cv4,'
data_mu_2023preBPix='data_mu_Cv123,data_mu_Cv4,'

data_egamma_2023postBPix='data_egamma_D,'
data_muoneg_2023postBPix='data_muoneg_D,'
data_mu_2023postBPix='data_mu_D,'

bkg_dy='DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,'
bkg_wj='WtoLNu_madgraphMLM,WtoLNu_1J_madgraphMLM,WtoLNu_2J_madgraphMLM,WtoLNu_3J_madgraphMLM,WtoLNu_4J_madgraphMLM,'
bkg_vv='WW,WZ,ZZ,'
bkg_vvv='WWW_4F,WWZ_4F,WZZ,ZZZ,'
bkg_vh_htt='WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,ZHto2Tau_UncorrelatedDecay_UnFiltered,'
bkg_higgs='GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,VBFHto2Tau_UncorrelatedDecay_UnFiltered,'
bkg_top='TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,'
bkg_ttbar='TTto2L2Nu,TTto4Q,TTtoLNu2Q,'
signal='bbh_htt_100,h_ggf_htt_100'
signal_bbh='bbh_htt_60,bbh_htt_65,bbh_htt_70,bbh_htt_75,bbh_htt_80,bbh_htt_85,bbh_htt_90,bbh_htt_95,bbh_htt_100,bbh_htt_105,bbh_htt_110,bbh_htt_115,bbh_htt_120,bbh_htt_125,bbh_htt_130,bbh_htt_135,bbh_htt_140,bbh_htt_160,bbh_htt_180,bbh_htt_200,bbh_htt_250,bbh_htt_300,bbh_htt_350,bbh_htt_400,bbh_htt_450,bbh_htt_500,bbh_htt_600,bbh_htt_700,bbh_htt_800,bbh_htt_900,bbh_htt_1000,bbh_htt_1100,bbh_htt_1200,bbh_htt_1400,bbh_htt_1600,bbh_htt_1800,bbh_htt_2000,bbh_htt_2300,bbh_htt_2600,bbh_htt_2900,bbh_htt_3200,bbh_htt_3500,'
signal_ggf='h_ggf_htt_60,h_ggf_htt_65,h_ggf_htt_70,h_ggf_htt_75,h_ggf_htt_80,h_ggf_htt_85,h_ggf_htt_90,h_ggf_htt_95,h_ggf_htt_100,h_ggf_htt_105,h_ggf_htt_110,h_ggf_htt_115,h_ggf_htt_120,h_ggf_htt_125,h_ggf_htt_130,h_ggf_htt_135,h_ggf_htt_140,h_ggf_htt_160,h_ggf_htt_180,h_ggf_htt_200,h_ggf_htt_250,h_ggf_htt_300,h_ggf_htt_350,h_ggf_htt_400,h_ggf_htt_450,h_ggf_htt_500,h_ggf_htt_600,h_ggf_htt_700,h_ggf_htt_800,h_ggf_htt_900,h_ggf_htt_1000,h_ggf_htt_1100,h_ggf_htt_1200,h_ggf_htt_1400,h_ggf_htt_1600,h_ggf_htt_1800,h_ggf_htt_2000,h_ggf_htt_2300,h_ggf_htt_2600,h_ggf_htt_2900,h_ggf_htt_3200,h_ggf_htt_3500'

case $1 in

##############################
####### Combined datasets ############
##############################
"22and23_emu")
        config="run3_2022_preEE_emu,run3_2022_postEE_emu,run3_2023_preBPix_emu,run3_2023_postBPix_emu"
        bkgs="$bkg_dy$bkg_wj$bkg_vv$bkg_vvv$bkg_vh_htt$bkg_higgs$bkg_top$bkg_ttbar$signal,"
        datasets="$data_egamma_2022preEE$data_mu_2022preEE${bkgs}:"
        datasets="${datasets}$data_egamma_2022postEE$data_mu_2022postEE${bkgs}:"
        datasets="${datasets}$data_egamma_2023preBPix$data_mu_2023preBPix${bkgs}:"
        datasets="${datasets}$data_egamma_2023postBPix$data_mu_2023postBPix${bkgs}"
        categories='cat_emu_sr'
        processes="data,dy_lep,dy_tt,h_ggf_htt_sm_prod_sm,st,tt,h_vbf_htt_sm,vh_htt,wj,vv,vvv,$signal"
        variables=$variables_emu
        workflow='htcondor'
     ;;
"22_emu")
        config="run3_2022_preEE_emu"
        bkgs="$bkg_dy$bkg_wj$bkg_vv$bkg_vvv$bkg_vh_htt$bkg_higgs$bkg_top$bkg_ttbar$signal"
        datasets="$data_egamma_2022preEE$data_mu_2022preEE${bkgs}"
        categories='cat_emu_sr,cat_emu_sr__nj0,cat_emu_sr__nj1'
        processes="data,dy_lep,dy_tt,h_ggf_htt_sm_prod_sm,st,tt,h_vbf_htt_sm,vh_htt,wj,vv,vvv,$signal,"
        variables=$variables_emu
        workflow='htcondor'
     ;;
"22EE_emu")
        config="run3_2022_postEE_emu"
        bkgs="$bkg_dy$bkg_wj$bkg_vv$bkg_vvv$bkg_vh_htt$bkg_higgs$bkg_top$bkg_ttbar$signal"
        datasets="$data_egamma_2022postEE$data_mu_2022postEE${bkgs}"
        categories='cat_emu_sr,cat_emu_sr__nj0,cat_emu_sr__nj1'
        processes="data,dy_lep,dy_tt,h_ggf_htt_sm_prod_sm,st,tt,h_vbf_htt_sm,vh_htt,wj,vv,vvv,$signal,"
        variables=$variables_emu
        workflow='htcondor'
     ;;
"23_emu")
        config="run3_2023_preBPix_emu"
        bkgs="$bkg_dy$bkg_wj$bkg_vv$bkg_vvv$bkg_vh_htt$bkg_higgs$bkg_top$bkg_ttbar$signal"
        datasets="$data_egamma_2023preBPix$data_mu_2023preBPix${bkgs}"
        categories='cat_emu_sr' #,cat_emu_sr__nj0,cat_emu_sr__nj1'
        processes="data,dy_lep,dy_tt,h_ggf_htt_sm_prod_sm,st,tt,h_vbf_htt_sm,vh_htt,wj,vv,vvv,$signal,"
        variables=$variables_emu
        workflow='htcondor'
     ;;
"23_emu_shift")
        config="run3_2023_preBPix_emu"
        datasets="DYto2Tau_MLL_50_0J_amcatnloFXFX"
        categories='cat_emu_sr' #,cat_emu_sr__nj0,cat_emu_sr__nj1'
        processes="dy_lep"
        variables="D_zeta"
        workflow='local'
     ;;
"23EE_emu")
        config="run3_2023_postBPix_emu"
        bkgs="$bkg_dy$bkg_wj$bkg_vv$bkg_vvv$bkg_vh_htt$bkg_higgs$bkg_top$bkg_ttbar$signal"
        datasets="$data_egamma_2023postBPix$data_mu_2023postBPix${bkgs}"
        categories='cat_emu_sr' #,cat_emu_sr__nj0,cat_emu_sr__nj1'
        processes="data,dy_lep,dy_tt,h_ggf_htt_sm_prod_sm,st,tt,h_vbf_htt_sm,vh_htt,wj,vv,vvv,$signal,"
        variables=$variables_emu
        workflow='htcondor'
     ;;
"22preEE_lim")
        config="run3_2022_preEE_emu_limited"
        bkgs="$bkg_dy$bkg_wj$bkg_vv$bkg_vvv$bkg_vh_htt$bkg_higgs$bkg_top$bkg_ttbar$signal"
        datasets="$data_egamma_2022preEE$data_mu_2022preEE${bkgs}:"
        categories='cat_emu_sr'
        processes="data,dy_lep,dy_tt,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,wj,vv,vvv,$signal,"
        variables=$variables_emu
        workflow='local'
     ;;
"22and23_emu_sig")
        config="run3_2022_preEE_emu,run3_2022_postEE_emu,run3_2023_preBPix_emu,run3_2023_postBPix_emu"
        datasets='bbh_htt_3500,h_ggf_htt_60'
        categories='cat_emu_sr'
        processes='bbh_htt_3500,h_ggf_htt_60'
        variables=$variables_emu
        workflow='htcondor'
     ;;
#########################
####### 2022preEE #######
#########################        
    "run3_2022preEE_emu_jec_lim")
      config="run3_2022_preEE_emu_limited"
      processes='data,tt_dl,tt_fh,tt_sl'
      datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,TTtoLNu2Q,TTto2L2Nu,TTto4Q'
	  categories='cat_emu_sr' #$categories_emu
	  variables='leading_b_jet_eta,subleading_b_jet_eta,leading_b_jet_phi,subleading_b_jet_phi,leading_b_jet_pt,subleading_b_jet_pt,di_b_jet_delta_eta,mb_jb_j,di_b_jet_delta_eta,mb_jb_j'
	  workflow='local'
    ;;

    "run3_2022preEE_emu_DC")
      config="run3_2022_preEE_emu"
      processes='data,dy_lep,dy_tt,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,wj,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
      datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_madgraphMLM,WtoLNu_1J_madgraphMLM,WtoLNu_2J_madgraphMLM,WtoLNu_3J_madgraphMLM,WtoLNu_4J_madgraphMLM,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	    categories='cat_emu_sr__bdt_ggh_M100,cat_emu_sr__bdt_ggh_M100' #$categories_emu
	    variables='bdt_raw_score_ggh_M100,bdt_raw_score_bbh_M100'
	    workflow='local'
    ;;
    "run3_2022preEE_emu_signals")
      config="run3_2022_preEE_emu"
      processes=$mssm_ggF_signal,$mssm_bbh_signal
      datasets=$mssm_ggF_signal,$mssm_bbh_signal
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2022postEE_emu_signals")
      config="run3_2022_postEE_emu"
      processes=$mssm_ggF_signal,$mssm_bbh_signal
      datasets=$mssm_ggF_signal,$mssm_bbh_signal
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2023preBPix_emu_signals")
      config="run3_2023_preBPix_emu"
      processes=$mssm_ggF_signal,$mssm_bbh_signal
      datasets=$mssm_ggF_signal,$mssm_bbh_signal
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2023postBPix_emu_signals")
      config="run3_2023_postBPix_emu"
      processes=$mssm_ggF_signal,$mssm_bbh_signal
      datasets=$mssm_ggF_signal,$mssm_bbh_signal
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2022preEE_emu_D_zeta_check")
      config="run3_2022_preEE_emu"
      processes='data,dy_lep,dy_tt,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,wj,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
      datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_madgraphMLM,WtoLNu_1J_madgraphMLM,WtoLNu_2J_madgraphMLM,WtoLNu_3J_madgraphMLM,WtoLNu_4J_madgraphMLM,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	    categories=$categories_emu
	    variables="D_zeta_check"
	    workflow='htcondor'
    ;;

    #DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX
    "run3_2022preEE_emu_SM_Higgs")
        config="run3_2022_preEE_emu"
        processes='h_ggf_htt,h_vbf_htt,'
        datasets='VBFHto2Tau_UncorrelatedDecay_UnFiltered,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,'
	      categories=$categories_emu
	      variables=$variables_emu
	      workflow='local'
    ;;
    "run3_2022preEE_emu_FF")
        config="run3_2022_preEE_emu"
        data=$data_egamma_2022preEE$data_mu_2022preEE
        bkg_ewk=$bkg_ewk
        bkg_single_top=$bkg_single_top
        bkg_ttbar=$bkg_ttbar
        datasets=$data$bkg_ewk$bkg_single_top$bkg_ttbar
        processes='dy_lep,vv,tt,st,wj,data'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2022preEE_emu_sig")
        config="run3_2022_preEE_emu"
        data=$data_egamma_2022preEE$data_mu_2022preEE
        bkg_ewk=$bkg_ewk
        bkg_single_top=$bkg_single_top
        bkg_ttbar=$bkg_ttbar
        datasets=$mssm_signal
        processes='h_ggf_htt_masses'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;


################################
###### 2022postEE_limited ######
################################  
    "run3_2022postEE_emu_lim")
        config="run3_2022_postEE_emu_limited"
        processes='tt_dl,dy_tt_m50_0j'
        datasets='TTto2L2Nu,DYto2Tau_MLL_50_0J_amcatnloFXFX'
	     categories='cat_emu_sr' #__nj0,cat_emu_sr__nj1'
	     variables="D_zeta"
	     workflow='local'
    ;;

#########################
####### 2022postEE ######
#########################
    
    "run3_2022postEE_emu_jec")
      config="run3_2022_postEE_emu_limited"
      processes='data,tt_dl,tt_sl'
      datasets='data_egamma_E,data_egamma_F,data_egamma_G,data_mu_E,data_mu_F,data_mu_G,TTtoLNu2Q,TTto2L2Nu'
	    categories='cat_emu_sr' #$categories_emu
	    variables='N_b_jets'
	    workflow='local'
    ;;

    "run3_2022postEE_emu")
      config="run3_2022_postEE_emu"
      processes='data,dy_lep,dy_tt,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,wj,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
      datasets='data_egamma_E,data_egamma_F,data_egamma_G,data_mu_E,data_mu_F,data_mu_G,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_madgraphMLM,WtoLNu_1J_madgraphMLM,WtoLNu_2J_madgraphMLM,WtoLNu_3J_madgraphMLM,WtoLNu_4J_madgraphMLM,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	  categories='cat_emu_sr__bdt_ggh_M100,cat_emu_sr__bdt_ggh_M100' #$categories_emu
	  variables='bdt_raw_score_ggh_M100,bdt_raw_score_bbh_M100'		
      categories=$categories_emu	    
      workflow='local'
    ;;

    "run3_2022postEE_emu_2D")
        config="run3_2022_postEE_emu"
        data=$data_egamma_2022postEE$data_mu_2022postEE
        bkg_ewk=$bkg_ewk
        bkg_single_top=$bkg_single_top
        bkg_ttbar=$bkg_ttbar
        datasets=$bkg_ewk$bkg_single_top$bkg_ttbar
        processes='dy_lep,vv,tt,st,wj'
	      categories=$categories_emu
	      variables='leading_jet_eta-leading_jet_phi'
	      workflow='htcondor'
    ;;

###################################
####### 2023preBPix_limited #######
###################################
    "run3_2023preBPix_emu_lim")
        config="run3_2023_preBPix_emu_limited"
        processes='data,tt_dl'
        datasets='data_egamma_Cv123,TTto2L2Nu'
	     categories='cat_emu_sr'
	     variables="D_zeta"
	     workflow='local'
    ;;

##############################
####### 2023preBPix ##########
##############################
    "run3_2023preBPix_emu")
        config="run3_2023_preBPix_emu"
        processes='data,dy_lep,dy_tt,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,wj,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
        datasets='data_egamma_Cv123,data_egamma_Cv4,data_mu_Cv4,data_mu_Cv123,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_madgraphMLM,WtoLNu_1J_madgraphMLM,WtoLNu_2J_madgraphMLM,WtoLNu_3J_madgraphMLM,WtoLNu_4J_madgraphMLM,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	      categories='cat_emu_sr__bdt_ggh_M100,cat_emu_sr__bdt_ggh_M100' #$categories_emu
	      variables='bdt_raw_score_ggh_M100,bdt_raw_score_bbh_M100'
	      workflow='htcondor'
    ;;


####################################
####### 2023postBPix_limited #######
####################################
    "run3_2023postBPix_emu_lim")
        config="run3_2023_postBPix_emu_limited"
        processes='data,tt_dl'
        datasets='data_egamma_D,TTto2L2Nu'
	     categories='cat_emu_sr'
	     variables="D_zeta"
	     workflow='local'
    ;;

############################
####### 2023postBPix #######
############################
    "run3_2023postBPix_emu")
        config="run3_2023_postBPix_emu"
        processes='data,dy_lep,dy_tt,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,wj,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
        datasets='data_egamma_D,data_mu_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_madgraphMLM,WtoLNu_1J_madgraphMLM,WtoLNu_2J_madgraphMLM,WtoLNu_3J_madgraphMLM,WtoLNu_4J_madgraphMLM,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	      categories='cat_emu_sr__bdt_ggh_M100,cat_emu_sr__bdt_ggh_M100' #$categories_emu
	      variables='bdt_raw_score_ggh_M100,bdt_raw_score_bbh_M100'	
        # categories=$categories_emu
	      # variables=$variables_emu
	     workflow='htcondor'
    ;;

    "run3_2023postBPix_emu_FF")
        config="run3_2023_postBPix_emu"
        data=$data_egamma_2023postBPix$data_mu_2023postBPix
        bkg_ewk=$bkg_ewk
        bkg_single_top=$bkg_single_top
        bkg_ttbar=$bkg_ttbar
        datasets=$data$bkg_ewk$bkg_single_top$bkg_ttbar
        processes='dy_lep,vv,tt,st,wj,data'
        categories=$categories_emu
        variables=$variables_emu
        workflow='htcondor'
    ;;

    *)
    echo "Unknown run argument!"
    exit

esac
}