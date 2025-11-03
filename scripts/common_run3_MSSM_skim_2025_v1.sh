#!/bin/bash

set_common_vars() {

version="desy_dev"

data_egamma_2022preEE_list=( "data_egamma_C" "data_egamma_D" )
data_muoneg_2022preEE_list=( "data_muoneg_C" "data_muoneg_D" )
data_mu_2022preEE_list=( "data_mu_C" "data_mu_D" "data_singlemu_C" )
data_egamma_2022postEE_list=( "data_egamma_E" "data_egamma_F" "data_egamma_G" )
data_muoneg_2022postEE_list=( "data_muoneg_E" "data_muoneg_F" "data_muoneg_G" )
data_mu_2022postEE_list=( "data_mu_E" "data_mu_F" "data_mu_G" )
data_egamma_2023preBPix_list=( "data_egamma_Cv123" "data_egamma_Cv4" )
data_egamma_2023postBPix_list=( "data_egamma_2023_D" )
data_muoneg_2023preBPix_list=( "data_muoneg_Cv123" "data_muoneg_Cv4" )
data_muoneg_2023postBPix_list=( "data_muoneg_2023_D" )
data_mu_2023preBPix_list=( "data_mu_Cv4" "data_mu_Cv123" )
data_mu_2023postBPix_list=( "data_mu_2023_D" )
categories_emu_list=( "cat_emu_sr" ) # "cat_emu_sr__bdt_sig_M100") # "cat_emu_sr__bdt_dy_M100" "cat_emu_sr__bdt_tt_M100" "cat_emu_sr__bdt_wj_M100" )
categories_emu_qcd_list=( "cat_emu_ar_qcd" "cat_emu_dr_num_qcd" "cat_emu_dr_den_qcd" )
variables_emu_list=( 
  "emu_mt_tot" "emu_mt_emu" "D_zeta"
  "emu_mt_e" "emu_mt_mu" "N_jets_pT_20_eta_4_7_Tight"
  "leading_jet_eta" "subleading_jet_eta" "leading_jet_phi"
  "subleading_jet_phi" "N_b_jets"
  "leading_jet_pt" "subleading_jet_pt" "delta_eta_jj" "mjj" 
  "emu_lep0_pt" "emu_lep0_eta" "emu_lep0_phi"
  "emu_lep0_ip_sig" "emu_lep1_pt" "emu_lep1_eta" "emu_lep1_phi" "emu_lep1_ip_sig"
  "emu_mvis" "emu_delta_r" "emu_pt" "puppi_met_pt" "puppi_met_phi"
  "bdt_raw_score_sig_M100" "bdt_raw_score_dy_M100" "bdt_raw_score_tt_M100" "bdt_raw_score_wj_M100" #"pt_H" "hcand_emu_fastMTT_mass" 
)
bkg_SM_higgs_list=(
  "GluGluHTo2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay"
  "VBFHToTauTau_UncorrelatedDecay_UnFiltered"
)
bkg_ewk_list=(
  "WW" "WZ" "ZZ" "WWW_4F" "WWZ_4F" "WZZ" "ZZZ"
  "DYto2Tau_MLL_50_0J_amcatnloFXFX"
  "DYto2Tau_MLL_50_1J_amcatnloFXFX"
  "DYto2Tau_MLL_50_2J_amcatnloFXFX"
  # "DYto2L_M_10to50_amcatnloFXFX"
  # "DYto2L_M_50_0J_amcatnloFXFX"
  # "DYto2L_M_50_1J_amcatnloFXFX"
  # "DYto2L_M_50_2J_amcatnloFXFX"
  "DYto2L_M_50_amcatnloFXFX"
  "WtoLNu_amcatnloFXFX"    
)
bkg_vh_htt_list=(
  "WminusHToTauTau_UncorrelatedDecay_UnFiltered"
  "WplusHToTauTau_UncorrelatedDecay_UnFiltered"
  "ZHToTauTau_UncorrelatedDecay_UnFiltered"
)
bkg_single_top_list=(
  "TbarWplusto4Q"
  "TWminusto4Q"
  "TbarWplusto2L2Nu" 
  "TbarWplustoLNu2Q"
  "TWminusto2L2Nu" 
  "TWminustoLNu2Q" 
)
bkg_ttbar_list=(
  "TTtoLNu2"
  "TTto2L2Nu"
  "TTto4Q"
)
mssm_ggF_signal_list=(
  "h_ggf_htt_60" "h_ggf_htt_65" "h_ggf_htt_70" "h_ggf_htt_75"
  "h_ggf_htt_80" "h_ggf_htt_85" "h_ggf_htt_90" "h_ggf_htt_95"
  "h_ggf_htt_100" "h_ggf_htt_105" "h_ggf_htt_110" "h_ggf_htt_115"
  "h_ggf_htt_120" "h_ggf_htt_125" "h_ggf_htt_130" "h_ggf_htt_135"
  "h_ggf_htt_140" "h_ggf_htt_160" "h_ggf_htt_180" "h_ggf_htt_200"
  "h_ggf_htt_250" "h_ggf_htt_300" "h_ggf_htt_350" "h_ggf_htt_400"
  "h_ggf_htt_450" "h_ggf_htt_500" "h_ggf_htt_600" "h_ggf_htt_700"
  "h_ggf_htt_800" "h_ggf_htt_900" "h_ggf_htt_1000" "h_ggf_htt_1100"
  "h_ggf_htt_1200" "h_ggf_htt_1400" "h_ggf_htt_1600" "h_ggf_htt_1800"
  "h_ggf_htt_2000" "h_ggf_htt_2300" "h_ggf_htt_2600" "h_ggf_htt_2900"
  "h_ggf_htt_3200" "h_ggf_htt_3500"
)
mssm_bbh_signal_list=(
  "bbh_htt_60" "bbh_htt_65" "bbh_htt_70" "bbh_htt_75" "bbh_htt_80" "bbh_htt_85"
  "bbh_htt_90" "bbh_htt_95" "bbh_htt_100" "bbh_htt_105" "bbh_htt_110" "bbh_htt_115"
  "bbh_htt_120" "bbh_htt_125" "bbh_htt_130" "bbh_htt_135" "bbh_htt_140" "bbh_htt_160"
  "bbh_htt_180" "bbh_htt_200" "bbh_htt_250" "bbh_htt_300" "bbh_htt_350" "bbh_htt_400"
  "bbh_htt_450" "bbh_htt_500" "bbh_htt_600" "bbh_htt_700" "bbh_htt_800" "bbh_htt_900"
  "bbh_htt_1000" "bbh_htt_1100" "bbh_htt_1200" "bbh_htt_1400" "bbh_htt_1600"
  "bbh_htt_1800" "bbh_htt_2000" "bbh_htt_2300" "bbh_htt_2600" "bbh_htt_2900"
  "bbh_htt_3200" "bbh_htt_3500"
)
# --- Build CSV strings ---
data_egamma_2022preEE=$(IFS=,; echo "${data_egamma_2022preEE_list[*]}")
data_muoneg_2022preEE=$(IFS=,; echo "${data_muoneg_2022preEE_list[*]}")
data_mu_2022preEE=$(IFS=,; echo "${data_mu_2022preEE_list[*]}")
data_egamma_2022postEE=$(IFS=,; echo "${data_egamma_2022postEE_list[*]}")
data_muoneg_2022postEE=$(IFS=,; echo "${data_muoneg_2022postEE_list[*]}")
data_mu_2022postEE=$(IFS=,; echo "${data_mu_2022postEE_list[*]}")
data_egamma_2023preBPix=$(IFS=,; echo "${data_egamma_2023preBPix_list[*]}")
data_egamma_2023postBPix=$(IFS=,; echo "${data_egamma_2023postBPix_list[*]}")
data_muoneg_2023preBPix=$(IFS=,; echo "${data_muoneg_2023preBPix_list[*]}")
data_muoneg_2023postBPix=$(IFS=,; echo "${data_muoneg_2023postBPix_list[*]}")
data_mu_2023preBPix=$(IFS=,; echo "${data_mu_2023preBPix_list[*]}")
data_mu_2023postBPix=$(IFS=,; echo "${data_mu_2023postBPix_list[*]}")
categories_emu=$(IFS=,; echo "${categories_emu_list[*]}")
categories_emu_qcd=$(IFS=,; echo "${categories_emu_qcd_list[*]}")
variables_emu=$(IFS=,; echo "${variables_emu_list[*]}")
bkg_SM_higgs=$(IFS=,; echo "${bkg_SM_higgs_list[*]}")
bkg_ewk=$(IFS=,; echo "${bkg_ewk_list[*]}")
bkg_vh_htt=$(IFS=,; echo "${bkg_vh_htt_list[*]}")
bkg_single_top=$(IFS=,; echo "${bkg_single_top_list[*]}")
bkg_ttbar=$(IFS=,; echo "${bkg_ttbar_list[*]}")
mssm_ggF_signal=$(IFS=,; echo "${mssm_ggF_signal_list[*]}")
mssm_bbh_signal=$(IFS=,; echo "${mssm_bbh_signal_list[*]}")

case $1 in
################################
####### 2022preEE_limited ######
################################
    "run3_2022preEE_emu_DY_split")
      config="run3_2022_preEE_emu"
      datasets='DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX' 
      processes='dy_ll,dy_tautau_nj' #'dy_z2ee,dy_z2mumu,dy_z2tautau'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2022postEE_emu_DY_split")
      config="run3_2022_postEE_emu_limited"
      datasets='DYto2L_M_50_amcatnloFXFX' 
      processes='dy_m50toinf' #'dy_z2ee,dy_z2mumu,dy_z2tautau'
	    categories='cat_emu_sr'
	    variables='D_zeta'
	    workflow='local'
    ;;
    "run3_2023preBPix_emu_DY_split")
      config="run3_2023_preBPix_emu_limited"
      datasets='DYto2L_M_50_amcatnloFXFX' 
      processes='dy_m50toinf' #'dy_z2ee,dy_z2mumu,dy_z2tautau'
	    categories='cat_emu_sr'
	    variables='D_zeta'
	    workflow='local'
    ;;
    "run3_2023postBPix_emu_DY_split")
      config="run3_2023_postBPix_emu_limited"
      datasets='DYto2L_M_50_amcatnloFXFX' 
      processes='dy_m50toinf' #'dy_z2ee,dy_z2mumu,dy_z2tautau'
	    categories='cat_emu_sr'
	    variables='D_zeta'
	    workflow='local'
    ;;
    "run3_2022preEE_emu_DY_incl_ALL")
      config="run3_2022_preEE_emu"
      processes='data,dy_m50toinf,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
      datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_50_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2022preEE_emu_03_test")
      config="run3_2022_preEE_emu_limited"
      processes='data,dy,dy_tautau_nj,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
      datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	  categories=$categories_emu
	  variables=$variables_emu
	  workflow='local'
    ;;
    "run3_2022preEE_emu_BSM_sig")
      config="run3_2022_preEE_emu"
      datasets="bbh_htt_60"
      processes="bbh_htt_60"
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='local'
    ;;

    "run3_2022postEE_emu_BSM_sig")
      config="run3_2022_postEE_emu"
      datasets="bbh_htt_60"
      processes="bbh_htt_60"
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='local'
    ;;
    "run3_2023preBPix_emu_BSM_sig")
      config="run3_2023_preBPix_emu"
      datasets="bbh_htt_60"
      processes="bbh_htt_60"
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='local'
    ;;
    "run3_2023postBPix_emu_BSM_sig")
      config="run3_2023_postBPix_emu"
      datasets="bbh_htt_60"
      processes="bbh_htt_60"
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='local'
    ;;

#########################
####### 2022preEE #######
#########################        

    "run3_2022preEE_emu")
      config="run3_2022_preEE_emu"
      processes='data,dy_ll,dy_tautau_nj,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
      datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2022preEE_emu_no_dy_tautau_1_2j")
        config="run3_2022_preEE_emu"
        processes='data,dy_ll,dy_tautau_m50toinf_0j,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
        datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2022preEE_emu_no_10_50")
        config="run3_2022_preEE_emu"
        processes='data,dy_ll,dy_tautau_nj,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
        datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2022preEE_emu_tt")
        config="run3_2022_preEE_emu"
        processes='tt,'
        datasets='TTto2L2Nu,TTto4Q,TTtoLNu2Q,'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
    "run3_2022preEE_emu_nodata")
        config="run3_2022_preEE_emu"
        processes='dy_ll,dy_tautau_nj,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
        datasets='DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;

    "run3_2022preEE_emu_no_LL")
        config="run3_2022_preEE_emu"
        processes='data,dy_tautau_nj,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
        categories=$categories_emu
        datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
        variables=$variables_emu
        workflow='local'
    ;;
    "run3_2022preEE_emu_no_tautau")
        config="run3_2022_preEE_emu"
        processes='data,dy,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
        categories=$categories_emu
        datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
        variables=$variables_emu
        workflow='local'
    ;;
    "run3_2022preEE_emu_BSM_ggF")
        config="run3_2022_preEE_emu"
        processes='h_ggf_htt_100,'
        categories=$categories_emu
        datasets='h_ggf_htt_100,'
        variables=$variables_emu
        workflow='local'
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
    "run3_2022preEE_emu_Z")
        config="run3_2022_preEE_emu"
        processes='dy_m10to50,dy_m50toinf_0j,dy_m50toinf_1j,dy_m50toinf_2j,dy_m50toinf,dy_tautau_m50toinf_0j,dy_tautau_m50toinf_1j,dy_tautau_m50toinf_2j'
        datasets='DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX'
	      categories=$categories_emu
	      variables=$variables_emu
	      workflow='htcondor'
    ;;
    "run3_2022preEE_emu_Z_Veto")
        config="run3_2022_preEE_emu"
        processes='dy_m50toinf_0j,dy_m50toinf_1j,dy_m50toinf_2j'
        datasets=',DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,'
	      categories=$categories_emu
	      variables=$variables_emu
	      workflow='htcondor'
    ;;
    "run3_2022preEE_emu_Z_incl")
        config="run3_2022_preEE_emu"
        processes='dy_m10to50,dy_m50toinf'
        datasets='DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX'
	      categories=$categories_emu
	      variables=$variables_emu
	      workflow='htcondor'
    ;;
    "run3_2022preEE_emu_0J")
        config="run3_2022_preEE_emu"
        processes='dy_m50toinf_0j,dy_tautau_m50toinf_0j'
        datasets='DYto2L_M_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,'
	      categories=$categories_emu
	      variables=$variables_emu
	      workflow='htcondor'
    ;;
    "run3_2022preEE_emu_1J")
        config="run3_2022_preEE_emu"
        processes='dy_m50toinf_1j,dy_tautau_m50toinf_1j'
        datasets='DYto2L_M_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,'
	      categories=$categories_emu
	      variables=$variables_emu
	      workflow='htcondor'
    ;;
    "run3_2022preEE_emu_2J")
        config="run3_2022_preEE_emu"
        processes='dy_m50toinf_2j,dy_tautau_m50toinf_2j'
        datasets='DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,'
	      categories=$categories_emu
	      variables=$variables_emu
	      workflow='htcondor'
    ;;
    "run3_2022preEE_emu_no_sticht")
        config="run3_2022_preEE_emu_limited"
        processes='dy_m10to50,dy_m50toinf,'
        datasets='DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,'
	      categories=$categories_emu
	      variables=$variables_emu
	      workflow='local'
    ;;
      "run3_2022preEE_emu_not_full")
        config="run3_2022_preEE_emu"
        processes='data,dy,h_ggf_htt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
        datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	      categories=$categories_emu
	      variables=$variables_emu
	      workflow='htcondor'
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
        datasets='WtoLNu_amcatnloFXFX'
        processes='w_lnu'
	    categories='cat_emu_sr__bdt_sig_M60'
	    variables='bdt_raw_score_sig_M60'
	    workflow='htcondor'
    ;;
    "run3_2022postEE_emu_dy_tautau_lim")
      config="run3_2022_postEE_emu_limited"
      datasets='DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX'
      processes='dy_tautau_nj'
	    categories='cat_emu_sr__bdt_sig_M60'
	    variables='bdt_raw_score_sig_M60'
	    workflow='local'
    ;;


#########################
####### 2022postEE ######
#########################

    "run3_2022postEE_emu")
        config="run3_2022_postEE_emu"
        processes='data,dy_ll,dy_tautau_nj,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,bbh_htt_100,h_ggf_htt_100,'
        datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,bbh_htt_100,h_ggf_htt_100,'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
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
        data=$data_egamma_2023preBPix,$data_mu_2023preBPix
        bkg_SM_higgs=$bkg_SM_higgs
        bkg_vh_htt=$bkg_vh_htt
        bkg_ewk=$bkg_ewk
        bkg_single_top=$bkg_single_top
        bkg_ttbar=$bkg_ttbar
        mssm_signal=$mssm_ggF_signal,$mssm_bbh_signal
        datasets=$data,$bkg_SM_higgs,$bkg_vh_htt,$bkg_ewk,$bkg_single_top,$bkg_ttbar,$mssm_signal
        processes='data,SM_higgs,vh_htt,dy_ll,dy_tautau,w_lnu,vv,vvv,st,tt,h_ggf_htt_masses,bbh_htt_masses'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;
##############################
####### 2023preBPix ##########
##############################
    "run3_2023preBPix_emu")
        config="run3_2023_preBPix_emu"
        processes='data,dy_ll,dy_tautau_nj,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,'
        datasets='data_egamma_C,data_egamma_D,data_singlemu_C,data_mu_C,data_mu_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;

####################################
####### 2023postBPix_limited #######
####################################

    "run3_2023postBPix_emu_lim")
        config="run3_2023_postBPix_emu_limited"
        data=$data_egamma_2023postBPix,$data_mu_2023postBPix
        bkg_SM_higgs=$bkg_SM_higgs
        bkg_vh_htt=$bkg_vh_htt
        bkg_ewk=$bkg_ewk
        bkg_single_top=$bkg_single_top
        bkg_ttbar=$bkg_ttbar
        mssm_signal=$mssm_ggF_signal,$mssm_bbh_signal
        datasets=$data,$bkg_SM_higgs,$bkg_vh_htt,$bkg_ewk,$bkg_single_top,$bkg_ttbar,$mssm_signal
        processes='data,SM_higgs,vh_htt,dy_ll,dy_tautau,w_lnu,vv,vvv,st,tt,h_ggf_htt_masses,bbh_htt_masses'
        categories=$categories_emu
        variables=$variables_emu
        workflow='htcondor'
    ;;
    "run3_2023postBPix_emu")
        config="run3_2023_postBPix_emu"
        processes='data,dy_ll,dy_tautau_nj,h_ggf_htt,st,tt,h_vbf_htt,vh_htt,w_lnu,vv,vvv,'
        datasets='data_egamma_2023_D,data_mu_2023_D,DYto2L_M_10to50_amcatnloFXFX,DYto2L_M_50_amcatnloFXFX,DYto2L_M_50_0J_amcatnloFXFX,DYto2L_M_50_1J_amcatnloFXFX,DYto2L_M_50_2J_amcatnloFXFX,DYto2Tau_MLL_50_0J_amcatnloFXFX,DYto2Tau_MLL_50_1J_amcatnloFXFX,DYto2Tau_MLL_50_2J_amcatnloFXFX,GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay,TbarWplusto4Q,TWminusto4Q,TbarWplusto2L2Nu,TbarWplustoLNu2Q,TWminusto2L2Nu,TWminustoLNu2Q,TTto2L2Nu,TTto4Q,TTtoLNu2Q,VBFHto2Tau_UncorrelatedDecay_UnFiltered,WminusHto2Tau_UncorrelatedDecay_UnFiltered,WplusHto2Tau_UncorrelatedDecay_UnFiltered,WtoLNu_amcatnloFXFX,WW,WZ,ZZ,WWW_4F,WWZ_4F,WZZ,ZZZ,ZHto2Tau_UncorrelatedDecay_UnFiltered,'
	    categories=$categories_emu
	    variables=$variables_emu
	    workflow='htcondor'
    ;;

############################
####### 2023postBPix #######
############################

    # "run3_2023postBPix_emu")
    #     config="run3_2023_postBPix_emu"
    #     data=$data_egamma_2023postBPix,$data_mu_2023postBPix
    #     bkg_SM_higgs=$bkg_SM_higgs
    #     bkg_vh_htt=$bkg_vh_htt
    #     bkg_ewk=$bkg_ewk
    #     bkg_single_top=$bkg_single_top
    #     bkg_ttbar=$bkg_ttbar
    #     mssm_signal=$mssm_ggF_signal,$mssm_bbh_signal
    #     datasets=$data,$bkg_SM_higgs,$bkg_vh_htt,$bkg_ewk,$bkg_single_top,$bkg_ttbar,$mssm_signal
    #     processes='data,SM_higgs,vh_htt,dy_ll,dy_tautau,w_lnu,vv,vvv,st,tt,h_ggf_htt_masses,bbh_htt_masses'
    #     categories=$categories_emu
    #     variables=$variables_emu
    #     workflow='htcondor'
    # ;;
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
