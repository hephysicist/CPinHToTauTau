#!/bin/bash
source ./cp_samples_2025_v1.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
prod_version=apr_28
args=(
        --configs $configs
        #--categories 'cat_mutau_sr__fake_incl__hig__cat0__tau2a1_3pr,'`
        #`'cat_mutau_sr__fake_incl__hig__cat1__tau2a1_3pr,'`
        #`'cat_mutau_sr__fake_incl__hig__cat2__tau2a1_3pr,'`
        #`'cat_mutau_sr__fake_incl__gtau,cat_mutau_sr__fake_incl__fake'


        --cf.CalibrateEvents-workflow $workflow
        --cf.CalibrateEvents-version $version

        --cf.SelectEvents-workflow $workflow
        --cf.SelectEvents-version $version
        
        --cf.ReduceEvents-workflow $workflow
        --cf.ReduceEvents-version $version
        --cf.ReduceEvents-tasks-per-job 1 

        --cf.MergeReducedEvents-workflow $workflow
        --cf.MergeReducedEvents-version $version

        
        --cf.MergeSelectionStats-version $version
        --cf.ProvideReducedEvents-version $version

        --cf.ProduceColumns-workflow $workflow
        --cf.ProduceColumns-version $prod_version
        --cf.ProduceColumns-tasks-per-job 3

        --cf.CreateHistograms-workflow $workflow
        --cf.CreateHistograms-version $prod_version
        --cf.CreateHistograms-tasks-per-job 5
        --cf.MergeHistograms-workflow $workflow
        --cf.MergeHistograms-version  $prod_version
        --cf.MergeHistograms-tasks-per-job 5
        --cf.MergeShiftedHistograms-version $prod_version
        --cf.MergeShiftedHistograms-tasks-per-job 1
        --cf.MergeShiftedHistograms-workflow $workflow
        --version $prod_version
       
        --inference-model hcp_model
        --hist-hooks ff_method
        "${@:2}"
    )
echo law run cf.CreateDatacards "${args[@]}"
law run cf.CreateDatacards "${args[@]}"
