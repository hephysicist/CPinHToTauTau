#!/bin/bash
source ./common_run3.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes  $processes
        --datasets $datasets
        
        --categories  $categories
        
        
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
        
        --cf.MergeReducedEvents-workflow $workflow
        --cf.MergeReducedEvents-version $version
        
        --cf.MergeSelectionStats-version $version
        --cf.ProvideReducedEvents-version $version
        --version ff_exp_and_p2_tau_eta2p3
        --variables $variables
        --file-types pdf,png
       # Currently there are three methods that are implemented as hist hooks:
       # 1. ff_method: general fake factor method that requires ff weights to be present at the events tree
       # 2. ff_method_dr_closure_test: Calclulate fake contribution and apply it to the dr_num regions for the closure tests
       # 3. good_old_abcd: estimates QCD contribution by taking events from same sign region and transfer factors from inv. lep iso
         
        --hist-hooks good_old_abcd #ff_method_dr_closure_test
        --general-settings "cms-label=pw,yscale=log"
        --process-settings "h_ggf_htt,unstack,scale=stack,color=#0000FF"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"
