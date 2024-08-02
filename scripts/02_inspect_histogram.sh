#!/bin/bash
args=(
    /afs/cern.ch/user/s/stzakhar/work/higgs_cp/data/cf_store/analysis_higgs_cp/cf.MergeHistograms/run3_2022_postEE_nano_tau_v12_limited/data_mu_f/nominal/calib__example/sel__default/prod__example/0/hist__muon_eta.pickle
)

cf_inspect "${args[@]}"

