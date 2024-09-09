#!/bin/bash

set_common_vars() {

version="no_mt_cut"
case $1 in
    "run3ts_lim")
        config="run3_2022_preEE_tau_spinner_limited"
        datasets='h_ggf_htt_filtered' #,h_ggf_htt_unfiltered,zh_htt_unfiltered'
        processes='h_ggf_htt,data_mu' #,zh_htt'
        workflow='local'
    ;;
    "run3ts")
        config="run3_2022_preEE_tau_spinner"
        datasets='h_ggf_htt_filtered' #,h_ggf_htt_unfiltered,zh_htt_unfiltered'
        processes='h_ggf_htt' #,zh_htt'
        workflow='local'
    ;;
    "run3_lim")
        config="run3_2022_preEE_limited"
        datasets='dy_incl' #,h_ggf_htt_unfiltered,zh_htt_unfiltered'
        processes='dy_lep' #,zh_htt'
        workflow='local'
    ;;
    "run3")
        config="run3_2022_preEE"
        data='data_mu_C,data_mu_D,data_tau_C,data_tau_D,data_e_C,data_e_D,'
        bkg_ewk='wj_incl,ww,wz,zz,dy_incl,'
        bkg_top='st_twchannel_t_fh,st_twchannel_t_sl,st_twchannel_t_dl,st_twchannel_tbar_sl,st_twchannel_tbar_dl,'
        bkg_ttbar='tt_sl,tt_dl,tt_fh'
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes='dy_lep,vv,tt,st,wj,data'
        workflow='htcondor'
    ;;
    # "run3")
    #     config="run3_2022_preEE_nano_tau_v12"
    #     #Datasets to use
    #     data='data_mu_c,data_mu_d,data_mu_e,'
    #     bkg_ewk='wj_incl,ww,wz,zz,wj, dy_incl,'
    #     bkg_top='st_t_bbarq,st_tbar_bq,'`
    #     `'st_t_wminus_to_lnu2q,st_t_wminus_to_2l2nu,'`
    #     `'st_tbar_wplus_to_lnu2q,st_tbar_wplus_to_2l2nu,'
    #     bkg_ttbar='tt_sl,tt_dl,tt_fh'
    #     datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
    #     processes="dy_z2mumu,dy_z2tautau,vv,tt,st,wj,data"
    # ;;
    *)
    echo "Unknown run argument!"
    exit
esac
}
