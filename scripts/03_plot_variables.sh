#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes 'dy_lep,vv,tt,st,wj,data'
        #--datasets $datasets
        --version $version
        --categories tautau_signal_reg #'mutau_signal_reg,mutau_signal_reg_no_mt'
        --cf.CalibrateEvents-workflow $workflow
        --cf.SelectEvents-workflow $workflow
        --cf.ReduceEvents-workflow $workflow
        --cf.MergeReducedEvents-workflow local
        --variables 
        #'phi_cp_mu_pi,phi_cp_mu_rho,'`
        #`'phi_cp_mu_pi_reg1,phi_cp_mu_pi_reg2,phi_cp_mu_rho_reg1,phi_cp_mu_rho_reg2,'`
        #`'hcand_leg1_ip_sig,hcand_leg2_ip_sig,alpha_mu_rho,alpha_mu_pi,'`
        #`'phi_cp_mu_pi_2bin,phi_cp_mu_pi_reg1_2bin,phi_cp_mu_pi_reg2_2bin,'`
        #`'phi_cp_mu_rho_2bin,phi_cp_mu_rho_reg1_2bin,phi_cp_mu_rho_reg2_2bin,'`
        #`'hcand_leg1_ip_x,hcand_leg2_ip_x,hcand_leg1_ip_y,hcand_leg2_ip_y,hcand_leg1_ip_z,hcand_leg2_ip_z'
        'phi_cp_rho_rho,phi_cp_rho_rho_reg1,phi_cp_rho_rho_reg2,'`
        `'phi_cp_rho_rho_2bin,phi_cp_rho_rho_reg1_2bin,phi_cp_rho_rho_reg2_2bin'

        #'puppi_met_phi,puppi_met_pt,mT,hcand_mass,hcand_leg1_pt,hcand_leg2_pt,hcand_leg1_eta,hcand_leg2_phi,hcand_leg1_phi,hcand_leg2_eta,hcand_leg1_mass,hcand_leg2_mass,hcand_leg2_decayModePNet,hcand_leg2_decayMode'
        #'mc_weight,pu_weight,tau_weight,muon_weight'
        #'mT,hcand_mass,hcand_leg1_pt,hcand_leg2_pt,hcand_leg1_eta,hcand_leg2_phi,hcand_leg1_phi,hcand_leg2_eta,hcand_leg1_mass,hcand_leg2_mass,hcand_leg2_decayModePNet,hcand_leg2_decayMode,mc_weight,pu_weight,tau_weight,muon_weight'
        
        
        --general-settings "cms-label=pw"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"