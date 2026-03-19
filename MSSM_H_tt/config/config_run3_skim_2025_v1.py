
# coding: utf-8

"""
Configuration of the MSSM analysis.
"""

import functools
import yaml
import law
import order as od
import os
from scinum import Number

from columnflow.util import DotDict, maybe_import, dev_sandbox
from columnflow.config_util import (
    get_root_processes_from_campaign, 
    add_category, add_shift_aliases,
    verify_config_processes,get_shifts_from_sources
)

ak = maybe_import("awkward")

def add_run3(ana: od.Analysis,
             campaign: od.Campaign,
             channel               = None,
             config_name           = None,
             config_id             = None,
             limit_dataset_files   = None,) -> od.Config :

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)
    
    # create a config by passing the campaign, so id and name will be identical
    cfg = ana.add_config(campaign,
                        name  = config_name,
                        id    = config_id)

    # gather campaign data
    run = campaign.x.run
    cfg.x.year = campaign.x.year
    cfg.x.tag = campaign.x.tag
    year = cfg.x.year

    # validations
    if campaign.x.year == 2022:
        assert campaign.x.tag in ["preEE", "postEE"]
    if campaign.x.year == 2023:
        assert campaign.x.tag in ["preBPix", "postBPix"]

    # gather campaign data
    year = campaign.x.year
    year2 = year % 100

    def tag_caster(campaign: od.Campaign) -> str:
        #Helper function to cast campaign tags to the tags used in POG groups for the scale factors
        year = campaign.x.year
        tag = campaign.x.tag
        out_tag = ""
        e_sf_tag = ""
        e_scale_corrector = ""
        e_smearing_corrector = "" 
        if year in [2017,2018]  : out_tag = "UL"
        elif tag == "preEE"     : 
            out_tag = "Summer22"
            e_sf_tag = "2022Re-recoBCD"
            e_scale_corrector = "EGMScale_EleEtaR9_2022preEE" #"2022Re-recoBCD_ScaleJSON"
            e_smearing_corrector = "EGMSmearAndSyst_EleEtaR9_2022preEE" #"2022Re-recoBCD_SmearingJSON"
        elif tag == "postEE"    : 
            out_tag = "Summer22EE"
            e_sf_tag = "2022Re-recoE+PromptFG"
            e_scale_corrector = "EGMScale_EleEtaR9_2022postEE" #"2022Re-recoE+PromptFG_ScaleJSON"
            e_smearing_corrector = "EGMSmearAndSyst_EleEtaR9_2022postEE"#2022Re-recoE+PromptFG_SmearingJSON"
        elif tag == "preBPix"   : 
            out_tag = "Summer23"
            e_sf_tag = "2023PromptC"
            e_scale_corrector = "EGMScale_EleEtaR9_2023preBPIX" #"2022Re-recoE+PromptFG_ScaleJSON"
            e_smearing_corrector = "EGMSmearAndSyst_EleEtaR9_2023preBPIX" #2022Re-recoE+PromptFG_SmearingJSON"
        elif tag == "postBPix"  : 
            out_tag = "Summer23BPix"
            e_sf_tag = "2023PromptD"
            e_scale_corrector = "EGMScale_EleEtaR9_2023postBPIX" #"2022Re-recoE+PromptFG_ScaleJSON"
            e_smearing_corrector = "EGMSmearAndSyst_EleEtaR9_2023postBPIX" #2022Re-recoE+PromptFG_SmearingJSON"
        elif tag == "preVFP"    : out_tag = "preVFP_UL"
        elif tag == "postVFP"   : out_tag = "postVFP_UL"    
        return out_tag, e_sf_tag, e_scale_corrector, e_smearing_corrector, tag
        
    tag, e_sf_tag, e_scale_corrector, e_smearing_corrector, tau_tag = tag_caster(campaign)

    # Map Run 3 eras to ACD metadata release tags used in /cat/metadata
    acd_tag_map = {
        (2022, "preEE"):   "Run3-22CDSep23-Summer22-NanoAODv12",
        (2022, "postEE"):  "Run3-22EFGSep23-Summer22EE-NanoAODv12",
        (2023, "preBPix"): "Run3-23CSep23-Summer23-NanoAODv12",
        (2023, "postBPix"): "Run3-23DSep23-Summer23BPix-NanoAODv12",
    }

    try:
        json_acd_tag = acd_tag_map[(year, campaign.x.tag)]
    except KeyError as e:
        raise RuntimeError(f"No json_acd_tag mapping for year={year}, tag={campaign.x.tag}") from e

    # (optional) expose on cfg for downstream visibility/debugging
    cfg.x.json_acd_tag = json_acd_tag

    # add processes we are interested in
    
    process_names = [
        "data", 
        "data_mu",
        "data_tau",
        "data_e",
        "data_egamma",
        "data_muoneg",
        "data_singlemu",
        # DY->ll
        "dy_lep",
        "dy_lep_m10to50",
        "dy_ll_m50_0j",
        "dy_ll_m50_1j",
        "dy_ll_m50_2j",
        "dy_ll_m50",  
        # DY->tautau
        "dy_tt_m50",
        "dy_tt_m50_0j",
        "dy_tt_m50_1j",
        "dy_tt_m50_2j",
        # single top
        "st",
        #single top t-channel        
        "st_tchannel_tbar",
        "st_tchannel_t",
        #single top s-channel   
        "st_schannel_t_lep",
        "st_schannel_tbar_lep",
        # single top tW channel
        "st_twchannel_tbar_fh",
        "st_twchannel_t_fh",
        "st_twchannel_tbar_dl",
        "st_twchannel_tbar_sl",
        "st_twchannel_t_dl",
        "st_twchannel_t_sl",
        # tt
        "tt",
        "tt_dl",
        "tt_fh",
        "tt_sl",
        # SM H->tautau 
        "h_ggf_htt_sm_prod_sm",
        "h_vbf_htt_sm",
        # vh_htt
        "vh_htt",
        "zh_htt_flat",
        "wph_htt_flat",
        "wmh_htt_flat",
        # W + jets
        "wj",
        "wj_1j",
        "wj_2j",
        "wj_3j",
        "wj_4j",
        # vv 
        "vv",
        "ww",
        "wz",
        "zz",
        # vvv
        "vvv",
        "www",
        "wwz",
        "zzz"
        ]

    from MSSM_H_tt.config.mass_points import read_bdt_masses
    MASS_POINTS = read_bdt_masses()

    for mass in MASS_POINTS:
        process_names.append(f"h_ggf_htt_{mass}")
        process_names.append(f"bbh_htt_{mass}")

    for process_name in process_names:
        # add the process
        proc = cfg.add_process(procs.get(process_name))
        #for signal datasets create special tag
        if process_name.startswith("h_"):
            proc.add_tag("signal")
           
    # add datasets we need to study
    dataset_names_2022preEE = [
        #data
        "data_egamma_C",
        "data_egamma_D",
        "data_muoneg_C",
        "data_muoneg_D",
        "data_singlemu_C",
        "data_mu_C",
        "data_mu_D",
        "data_tau_C",
        "data_tau_D",
        # DY->ll
        "DYto2L_M_10to50_amcatnloFXFX",
        "DYto2L_M_50_0J_amcatnloFXFX",
        "DYto2L_M_50_1J_amcatnloFXFX",
        "DYto2L_M_50_2J_amcatnloFXFX",
        "DYto2L_M_50_amcatnloFXFX",
        # DY->tautau
        "DYto2Tau_MLL_50_0J_amcatnloFXFX",
        "DYto2Tau_MLL_50_1J_amcatnloFXFX",
        "DYto2Tau_MLL_50_2J_amcatnloFXFX",
        # SM ggH->tautau
        "GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay",
        # Single top
        "TbarWplusto4Q",
        "TWminusto4Q",
        "TbarWplusto2L2Nu", 
        "TbarWplustoLNu2Q",
        "TWminusto2L2Nu", 
        "TWminustoLNu2Q", 
        # tt
        "TTto2L2Nu",
        "TTto4Q",
        "TTtoLNu2Q",
        # SM VBFH->tautau
        "VBFHto2Tau_UncorrelatedDecay_UnFiltered",
        # WH->tautau
        "WminusHto2Tau_UncorrelatedDecay_UnFiltered",
        "WplusHto2Tau_UncorrelatedDecay_UnFiltered",
        # W+jets                        
        # "WtoLNu_amcatnloFXFX",                  
        "WtoLNu_1J_madgraphMLM",
        "WtoLNu_2J_madgraphMLM",
        "WtoLNu_3J_madgraphMLM",
        "WtoLNu_4J_madgraphMLM",
        "WtoLNu_madgraphMLM",
        # Diboson
        "WW",
        "WZ",
        "ZZ",
        # Triboson
        "WWW_4F",
        "WWZ_4F",
        "WZZ",
        "ZZZ",
        # ZH->tautau
        "ZHto2Tau_UncorrelatedDecay_UnFiltered",
        ]


    dataset_names_2022postEE = [
        #data
        "data_egamma_E",
        "data_egamma_F",
        "data_egamma_G",
        "data_muoneg_E",
        "data_muoneg_F",
        "data_muoneg_G",
        "data_mu_E",
        "data_mu_F",
        "data_mu_G",
        "data_tau_E",
        "data_tau_F",
        "data_tau_G",
        # DY->ll
        "DYto2L_M_10to50_amcatnloFXFX",
        "DYto2L_M_50_0J_amcatnloFXFX",
        "DYto2L_M_50_1J_amcatnloFXFX",
        "DYto2L_M_50_2J_amcatnloFXFX",
        "DYto2L_M_50_amcatnloFXFX",
        # DY->tautau
        "DYto2Tau_MLL_50_0J_amcatnloFXFX",
        "DYto2Tau_MLL_50_1J_amcatnloFXFX",
        "DYto2Tau_MLL_50_2J_amcatnloFXFX",
        # SM ggH->tautau
        "GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay",
        # Single top
        "TbarWplusto4Q",
        "TWminusto4Q",
        "TbarWplusto2L2Nu", 
        "TbarWplustoLNu2Q",
        "TWminusto2L2Nu", 
        "TWminustoLNu2Q", 
        # tt
        "TTto2L2Nu",
        "TTto4Q",
        "TTtoLNu2Q",
        # SM VBFH->tautau
        "VBFHto2Tau_UncorrelatedDecay_UnFiltered",
        # WH->tautau
        "WminusHto2Tau_UncorrelatedDecay_UnFiltered",
        "WplusHto2Tau_UncorrelatedDecay_UnFiltered",
        # W+jets                          
        # "WtoLNu_amcatnloFXFX",                  
        "WtoLNu_1J_madgraphMLM",
        "WtoLNu_2J_madgraphMLM",
        "WtoLNu_3J_madgraphMLM",
        "WtoLNu_4J_madgraphMLM",
        "WtoLNu_madgraphMLM",
        # Diboson
        "WW",
        "WZ",
        "ZZ",
        # Triboson
        "WWW_4F",
        "WWZ_4F",
        "WZZ",
        "ZZZ",
        # ZH->tautau
        "ZHto2Tau_UncorrelatedDecay_UnFiltered",
        ]
    
    dataset_names_2023preBPix = [
        #data
        "data_egamma_Cv123",
        "data_egamma_Cv4",
        "data_muoneg_Cv123",
        "data_muoneg_Cv4",
        "data_mu_Cv123",
        "data_mu_Cv4",
        # DY->ll
        "DYto2L_M_10to50_amcatnloFXFX",
        "DYto2L_M_50_0J_amcatnloFXFX",
        "DYto2L_M_50_1J_amcatnloFXFX",
        "DYto2L_M_50_2J_amcatnloFXFX",
        "DYto2L_M_50_amcatnloFXFX",
        # DY->tautau
        "DYto2Tau_MLL_50_0J_amcatnloFXFX",
        "DYto2Tau_MLL_50_1J_amcatnloFXFX",
        "DYto2Tau_MLL_50_2J_amcatnloFXFX",
        # SM ggH->tautau
        "GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay",
        # Single top
        "TbarWplusto4Q",
        "TWminusto4Q",
        "TbarWplusto2L2Nu", 
        "TbarWplustoLNu2Q",
        "TWminusto2L2Nu", 
        "TWminustoLNu2Q", 
        # tt
        "TTto2L2Nu",
        "TTto4Q",
        "TTtoLNu2Q",
        # SM VBFH->tautau
        "VBFHto2Tau_UncorrelatedDecay_UnFiltered",
        # WH->tautau
        "WminusHto2Tau_UncorrelatedDecay_UnFiltered",
        "WplusHto2Tau_UncorrelatedDecay_UnFiltered",
        # W+jets                          
        # "WtoLNu_amcatnloFXFX",                  
        "WtoLNu_1J_madgraphMLM",
        "WtoLNu_2J_madgraphMLM",
        "WtoLNu_3J_madgraphMLM",
        "WtoLNu_4J_madgraphMLM",
        "WtoLNu_madgraphMLM",
        # Diboson
        "WW",
        "WZ",
        "ZZ",
        # Triboson
        "WWW_4F",
        "WWZ_4F",
        "WZZ",
        "ZZZ",
        # ZH->tautau
        "ZHto2Tau_UncorrelatedDecay_UnFiltered",
        ]
    
    dataset_names_2023postBPix = [
        #data
        "data_egamma_D",
        "data_muoneg_D",
        "data_mu_D",
        # DY->ll
        "DYto2L_M_10to50_amcatnloFXFX",
        "DYto2L_M_50_0J_amcatnloFXFX",
        "DYto2L_M_50_1J_amcatnloFXFX",
        "DYto2L_M_50_2J_amcatnloFXFX",
        "DYto2L_M_50_amcatnloFXFX",
        # DY->tautau
        "DYto2Tau_MLL_50_0J_amcatnloFXFX",
        "DYto2Tau_MLL_50_1J_amcatnloFXFX",
        "DYto2Tau_MLL_50_2J_amcatnloFXFX",
        # SM ggH->tautau
        "GluGluHto2Tau_UncorrelatedDecay_SM_UnFiltered_ProdAndDecay",
        # Single top
        "TbarWplusto4Q",
        "TWminusto4Q",
        "TbarWplusto2L2Nu", 
        "TbarWplustoLNu2Q",
        "TWminusto2L2Nu", 
        "TWminustoLNu2Q", 
        # tt
        "TTto2L2Nu",
        "TTto4Q",
        "TTtoLNu2Q",
        # SM VBFH->tautau
        "VBFHto2Tau_UncorrelatedDecay_UnFiltered",
        # WH->tautau
        "WminusHto2Tau_UncorrelatedDecay_UnFiltered",
        "WplusHto2Tau_UncorrelatedDecay_UnFiltered",
        # W+jets                           
        # "WtoLNu_amcatnloFXFX",                  
        "WtoLNu_1J_madgraphMLM",
        "WtoLNu_2J_madgraphMLM",
        "WtoLNu_3J_madgraphMLM",
        "WtoLNu_4J_madgraphMLM",
        "WtoLNu_madgraphMLM",
        # Diboson
        "WW",
        "WZ",
        "ZZ",
        # Triboson
        "WWW_4F",
        "WWZ_4F",
        "WZZ",
        "ZZZ",
        # ZH->tautau
        "ZHto2Tau_UncorrelatedDecay_UnFiltered",
        ]
    for mass in MASS_POINTS:
        dataset_names_2022preEE.append(f"h_ggf_htt_{mass}")
        dataset_names_2022postEE.append(f"h_ggf_htt_{mass}")
        dataset_names_2023preBPix.append(f"h_ggf_htt_{mass}")
        dataset_names_2023postBPix.append(f"h_ggf_htt_{mass}")
        dataset_names_2022preEE.append(f"bbh_htt_{mass}")
        dataset_names_2022postEE.append(f"bbh_htt_{mass}")
        dataset_names_2023preBPix.append(f"bbh_htt_{mass}")
        dataset_names_2023postBPix.append(f"bbh_htt_{mass}")

    dataset_era = {
        "Summer22": dataset_names_2022preEE,
        "Summer22EE" : dataset_names_2022postEE,
        "Summer23" : dataset_names_2023preBPix,
        "Summer23BPix" : dataset_names_2023postBPix
    }
    dataset_names = dataset_era[tag]

    for dataset_name in dataset_names:
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))
        if dataset_name.startswith("h_ggf_htt_"):
            dataset.add_tag("signal")
            dataset.add_tag("is_mc")
        if dataset_name.startswith("bbh_htt_"):
            dataset.add_tag("signal")
            dataset.add_tag("is_mc")
        if dataset.name.startswith("TTto"):
            dataset.add_tag({"has_top", "ttbar", "tt"})    
        # for testing purposes, limit the number of files to 1
        for info in dataset.info.values():
            if limit_dataset_files:
                info.n_files = min(info.n_files, limit_dataset_files) #<<< REMOVE THIS FOR THE FULL DATASET

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

    #Adding the triggers 

    from MSSM_H_tt.config.triggers import add_triggers_run3
    if year>=2022: add_triggers_run3(cfg)
    
    #Adding the met filters 
    from MSSM_H_tt.config.met_filters import add_met_filters
    add_met_filters(cfg)

    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    cfg.x.default_calibrator = "main"
    cfg.x.default_selector = "main"
    cfg.x.default_reducer = "cf_default"
    cfg.x.default_producer = "main"
    cfg.x.default_weight_producer = "main"
    cfg.x.default_hist_producer = "httcp_hist_producer"
    cfg.x.default_ml_model = None
    cfg.x.default_inference_model = "example"
    cfg.x.default_categories = ("incl",)
    cfg.x.default_variables = ("event","channel_id")

    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {
        "data" : [
            "data_mu", 
            "data_tau",
            "data_e",
            "data_egamma",
            "data_muoneg",
            "data_singlemu"
            ],
        "vv": [
            "ww",
            "wz",
            "zz",
            ],
        "vvv": [
            "www",
            "wwz",
            "zzz"
            ],
        "tt": [
            "tt_dl",
            "tt_fh",
            "tt_sl",
            ],
        "st": [
            "TbarWplusto4Q",
            "TWminusto4Q",
            "TbarWplusto2L2Nu", 
            "TbarWplustoLNu2Q",
            "TWminusto2L2Nu", 
            "TWminustoLNu2Q", 
            ],
        "SM_higgs": [
            "h_ggf_htt_sm_prod_sm",
            "h_vbf_htt_sm",
            ],
        "vh_htt": [
            "zh_htt_flat",
            "wph_htt_flat",
            "wmh_htt_flat"
            ],
        "wj": [         
            "wj",
            "wj_1j",
            "wj_2j",
            "wj_3j",
            "wj_4j",     
            ],
        "dy_tt_m50": [
            "dy_tt_m50_0j",
            "dy_tt_m50_1j",
            "dy_tt_m50_2j",
            ],
        "dy_lep": [
            "dy_lep_m10to50",
            "dy_ll_m50_0j",
            "dy_ll_m50_1j",
            "dy_ll_m50_2j",
            "dy_ll_m50",  
            ],
        "h_ggf_htt_masses":["h_ggf_htt_60","h_ggf_htt_65","h_ggf_htt_70",
                            "h_ggf_htt_75","h_ggf_htt_80","h_ggf_htt_85",
                            "h_ggf_htt_90","h_ggf_htt_95","h_ggf_htt_100",
                            "h_ggf_htt_105","h_ggf_htt_110","h_ggf_htt_115",
                            "h_ggf_htt_120","h_ggf_htt_125","h_ggf_htt_130",
                            "h_ggf_htt_135","h_ggf_htt_140","h_ggf_htt_160",
                            "h_ggf_htt_180","h_ggf_htt_200","h_ggf_htt_250",
                            "h_ggf_htt_300","h_ggf_htt_350","h_ggf_htt_400",
                            "h_ggf_htt_450","h_ggf_htt_500","h_ggf_htt_600",
                            "h_ggf_htt_700","h_ggf_htt_800","h_ggf_htt_900",
                            "h_ggf_htt_1000","h_ggf_htt_1100","h_ggf_htt_1200",
                            "h_ggf_htt_1400","h_ggf_htt_1600","h_ggf_htt_1800",
                            "h_ggf_htt_2000","h_ggf_htt_2300","h_ggf_htt_2600",
                            "h_ggf_htt_2900","h_ggf_htt_3200","h_ggf_htt_3500"
                            ],
        "bbh_htt_masses": ["bbh_htt_60","bbh_htt_65","bbh_htt_70",
                            "bbh_htt_75","bbh_htt_80","bbh_htt_85",
                            "bbh_htt_90","bbh_htt_95","bbh_htt_100",
                            "bbh_htt_105","bbh_htt_110","bbh_htt_115",
                            "bbh_htt_120","bbh_htt_125","bbh_htt_130",
                            "bbh_htt_135","bbh_htt_140","bbh_htt_160",
                            "bbh_htt_180","bbh_htt_200","bbh_htt_250",
                            "bbh_htt_300","bbh_htt_350","bbh_htt_400",
                            "bbh_htt_450","bbh_htt_500","bbh_htt_600",
                            "bbh_htt_700","bbh_htt_800","bbh_htt_900",
                            "bbh_htt_1000","bbh_htt_1100","bbh_htt_1200",
                            "bbh_htt_1400","bbh_htt_1600","bbh_htt_1800",
                            "bbh_htt_2000","bbh_htt_2300","bbh_htt_2600",
                            "bbh_htt_2900","bbh_htt_3200","bbh_htt_3500"
                            ],
        }

    # dataset groups for conveniently looping over certain datasets
    # (used in wrapper_factory and during plotting)
    cfg.x.dataset_groups = {}

    # category groups for conveniently looping over certain categories
    # (used during plotting)
    cfg.x.category_groups = {}

    # variable groups for conveniently looping over certain variables
    # (used during plotting)
    cfg.x.variable_groups = {}


    # selector step groups for conveniently looping over certain steps
    # (used in cutflow tasks)
    cfg.x.selector_step_groups = {
        "default": ["json", "met_filter", "dl_res_veto", "trigger", "lepton", "jet"],
    }
    #  cfg.x.selector_step_labels = {"Initial":0, 
    #                                "Trigger": , "Muon"}
     
    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False
    
    # # jec configuration
    # # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=201
    # jerc_postfix = ""
    # if year == 2016 and campaign.x.vfp == "post":
    #     jerc_postfix = "APV"
    # elif year == 2022 and campaign.x.tag == "postEE":
    #     jerc_postfix = "EE"
    
    # elif year == 2023 and campaign.x.tag == "postBPix":
    #     jerc_postfix = "BPix"
    # jet_type = "AK4PFPuppi"
    # if year < 2022:
    #     jerc_campaign = f"Summer19UL{year2}{jerc_postfix}"
    #     jet_type = "AK4PFchs"
    # elif year==2022:
    #     jerc_campaign = f"Summer{year2}{jerc_postfix}_22Sep2023"
    #     jet_type = "AK4PFPuppi"
    # elif year==2023:
    #     jerc_campaign = f"Summer{year2}{jerc_postfix}Prompt23"
        
    # cfg.x.jec = DotDict.wrap({
    #     "campaign": jerc_campaign,
    #      "version": {2016: "V7", 2017: "V5", 2018: "V5", 2022: "V2", 2023:"V1"}[year],
    #     "jet_type": jet_type,
    #     "levels_DATA": ["L1L2L3Res"], #"L2Relative", "L2L3Residual", "L3Absolute", "L1L2L3Res" 
    #     "levels_MC": ["L1L2L3Res"], 
    #     "levels_for_type1_met": ["L1L2L3Res"], 
    #     "uncertainty_sources": [
    #         # "AbsoluteStat",
    #         # "AbsoluteScale",
    #         # "AbsoluteSample",
    #         # "AbsoluteFlavMap",
    #         # "AbsoluteMPFBias",
    #         # "Fragmentation",
    #         # "SinglePionECAL",
    #         # "SinglePionHCAL",
    #         # "FlavorQCD",
    #         # "TimePtEta",
    #         # "RelativeJEREC1",
    #         # "RelativeJEREC2",
    #         # "RelativeJERHF",
    #         # "RelativePtBB",
    #         # "RelativePtEC1",
    #         # "RelativePtEC2",
    #         # "RelativePtHF",
    #         # "RelativeBal",
    #         # "RelativeSample",
    #         # "RelativeFSR",
    #         # "RelativeStatFSR",
    #         # "RelativeStatEC",
    #         # "RelativeStatHF",
    #         # "PileUpDataMC",
    #         # "PileUpPtRef",
    #         # "PileUpPtBB",
    #         # "PileUpPtEC1",
    #         # "PileUpPtEC2",
    #         # "PileUpPtHF",
    #         # "PileUpMuZero",
    #         # "PileUpEnvelope",
    #         # "SubTotalPileUp",
    #         # "SubTotalRelative",
    #         # "SubTotalPt",
    #         # "SubTotalScale",
    #         # "SubTotalAbsolute",
    #         # "SubTotalMC",
    #         "Total",
    #         # "TotalNoFlavor",
    #         # "TotalNoTime",
    #         # "TotalNoFlavorNoTime",
    #         # "FlavorZJet",
    #         # "FlavorPhotonJet",
    #         # "FlavorPureGluon",
    #         # "FlavorPureQuark",
    #         # "FlavorPureCharm",
    #         # "FlavorPureBottom",
    #         # # "TimeRunA",
    #         # # "TimeRunB",
    #         # # "TimeRunC",
    #         # # "TimeRunD",
    #         # "CorrelationGroupMPFInSitu",
    #         # "CorrelationGroupIntercalibration",
    #         # "CorrelationGroupbJES",
    #         # "CorrelationGroupFlavor",
    #         # "CorrelationGroupUncorrelated",
    #     ],
    # })

    ################################################################################################
    # jet settings
    # TODO: keep a single table somewhere that configures all settings: btag correlation, year
    #       dependence, usage in calibrator, etc
    ################################################################################################
    # common jec/jer settings configuration
    if run == 3:
        # https://cms-jerc.web.cern.ch/Recommendations/#2022
        jerc_postfix = {
            (2022, ""): "_22Sep2023",
            (2022, "EE"): "_22Sep2023",
            (2023, ""): "Prompt23",
            (2023, "BPix"): "Prompt23",
            (2024, ""): "Prompt24",
        }[(year, campaign.x.postfix)]
        jec_campaign = f"Summer{year2}{campaign.x.postfix}{jerc_postfix}"
        jec_version = {
            (2022, ""): "V3",
            (2022, "EE"): "V3",
            (2023, ""): "V2",
            (2023, "BPix"): "V3",
            (2024, ""): "V2",
        }[(year, campaign.x.postfix)]
        jer_campaign = f"Summer{year2}{campaign.x.postfix}{jerc_postfix}"
        if year == 2024:
            jer_campaign = "Summer23BPixPrompt23"  # https://cms-jerc.web.cern.ch/Recommendations/#2024_1
        # special "Run" fragment in 2023 jer campaign
        if year == 2023:
            jer_campaign += f"_Run{'Cv1234' if campaign.x.tag =='preBPix' else 'D'}"
        if year == 2024:
            jer_campaign += "_RunD"
        jer_version = "JR" + {2022: "V1", 2023: "V1", 2024: "V1"}[year]
        jet_type = "AK4PFPuppi"
    else:
        assert False

    # full list of jec sources in a fixed order that is used to assign consistent ids across configs
    # (please add new sources at the bottom to preserve the order of existing ones)
    # the boolean flag decides whether to use them in the JEC config and if shifts should be created for them
    # https://cms-jerc.web.cern.ch/Recommendations/#uncertainites-and-correlations
    jec_source_era = f"{year}{campaign.x.postfix}"
    all_jec_sources = {
        "AbsoluteFlavMap": False,
        "AbsoluteMPFBias": False,
        "AbsoluteSample": False,
        "AbsoluteScale": False,
        "AbsoluteStat": False,
        "FlavorPhotonJet": False,
        "FlavorPureBottom": False,
        "FlavorPureCharm": False,
        "FlavorPureGluon": False,
        "FlavorPureQuark": False,
        "FlavorQCD": False,
        "FlavorZJet": False,
        "Fragmentation": False,
        "PileUpDataMC": False,
        "PileUpEnvelope": False,
        "PileUpMuZero": False,
        "PileUpPtBB": False,
        "PileUpPtEC1": False,
        "PileUpPtEC2": False,
        "PileUpPtHF": False,
        "PileUpPtRef": False,
        "RelativeBal": False,
        "RelativeFSR": False,
        "RelativeJEREC1": False,
        "RelativeJEREC2": False,
        "RelativeJERHF": False,
        "RelativePtBB": False,
        "RelativePtEC1": False,
        "RelativePtEC2": False,
        "RelativePtHF": False,
        "RelativeSample": False,
        "RelativeStatEC": False,
        "RelativeStatFSR": False,
        "RelativeStatHF": False,
        "SinglePionECAL": False,
        "SinglePionHCAL": False,
        "SubTotalAbsolute": False,
        "SubTotalMC": False,
        "SubTotalPileUp": False,
        "SubTotalPt": False,
        "SubTotalRelative": False,
        "SubTotalScale": False,
        "TimePtEta": False,
        "Total": True,
        "TotalNoFlavor": False,
        "TotalNoFlavorNoTime": False,
        "TotalNoTime": False,
        "CorrelationGroupFlavor": False,
        "CorrelationGroupIntercalibration": False,
        "CorrelationGroupMPFInSitu": False,
        "CorrelationGroupUncorrelated": False,
        "CorrelationGroupbJES": False,
        "Regrouped_Absolute": False,
        f"Regrouped_Absolute_{jec_source_era}": False,
        "Regrouped_BBEC1": False,
        f"Regrouped_BBEC1_{jec_source_era}": False,
        "Regrouped_EC2": False,
        f"Regrouped_EC2_{jec_source_era}": False,
        "Regrouped_FlavorQCD": False,
        "Regrouped_HF": False,
        f"Regrouped_HF_{jec_source_era}": False,
        "Regrouped_RelativeBal": False,
        f"Regrouped_RelativeSample_{jec_source_era}": False,
        "Regrouped_Total": False,
    }

    cfg.x.jec = DotDict.wrap({
        "Jet": {
            "campaign": jec_campaign,
            "version": jec_version,
            "data_per_era": year == 2022,  # 2022 JEC has the era in the correction set name
            "jet_type": jet_type,
            "levels": ["L1L2L3Res"], #["L1FastJet", "L2Relative", "L2L3Residual", "L3Absolute"],
            "levels_for_type1_met": ["L1L2L3Res"], # ["L1FastJet"],
            "uncertainty_sources": [src for src, flag in all_jec_sources.items() if flag],
        },
    })

    # JER
    cfg.x.jer = DotDict.wrap({
        "Jet": {
            "campaign": jer_campaign,
            "version": jer_version,
            "jet_type": jet_type,
        },
    })

    # updated jet id
    from columnflow.production.cms.jet import JetIdConfig
    cfg.x.jet_id = JetIdConfig(corrections={"AK4PUPPI_Tight": 2, "AK4PUPPI_TightLeptonVeto": 3})
    # trigger sf corrector
    cfg.x.jet_trigger_corrector = "jetlegSFs"
    

    # JEC uncertainty sources propagated to btag scale factors
    # (names derived from contents in BTV correctionlib file)
    cfg.x.btag_sf_jec_sources = [
        "",  # total
        "Absolute",
        "AbsoluteMPFBias",
        "AbsoluteScale",
        "AbsoluteStat",
        f"Absolute_{year}",
        "BBEC1",
        f"BBEC1_{year}",
        "EC2",
        f"EC2_{year}",
        "FlavorQCD",
        "Fragmentation",
        "HF",
        f"HF_{year}",
        "PileUpDataMC",
        "PileUpPtBB",
        "PileUpPtEC1",
        "PileUpPtEC2",
        "PileUpPtHF",
        "PileUpPtRef",
        "RelativeBal",
        "RelativeFSR",
        "RelativeJEREC1",
        "RelativeJEREC2",
        "RelativeJERHF",
        "RelativePtBB",
        "RelativePtEC1",
        "RelativePtEC2",
        "RelativePtHF",
        "RelativeSample",
        f"RelativeSample_{year}",
        "RelativeStatEC",
        "RelativeStatFSR",
        "RelativeStatHF",
        "SinglePionECAL",
        "SinglePionHCAL",
        "TimePtEta",
    ]
    ################################################################################################
    # met settings
    ################################################################################################


    if run == 3:
        cfg.x.met_name = "PuppiMET"
        cfg.x.raw_met_name = "RawPuppiMET"
        
    # ##################################
    # # Parameters fot top pT reweight #
    # ##################################
    # https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_the
    cfg.x.top_pt_reweighting_params = {
            "a": 0.0615,
            "a_up": 0.0615 * 1.5,
            "a_down": 0.0615 * 0.5,
            "b": -0.0005,
            "b_up": -0.0005 * 1.5,
            "b_down": -0.0005 * 0.5,
        }
    # top pt reweighting
    # https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting?rev=31
    # theory-based method preferred
    # from MSSM_H_tt.production.top_pt_weight import TopPtWeightFromTheoryConfig
    # cfg.x.top_pt_weight = TopPtWeightFromTheoryConfig(params={
    #     "a": 0.103,
    #     "b": -0.0118,
    #     "c": -0.000134,
    #     "d": 0.973,
    # })
    # data-based method preferred
    # from columnflow.production.cms.top_pt_weight import TopPtWeightFromDataConfig
    # cfg.x.top_pt_weight = TopPtWeightFromDataConfig(
    #     params={
    #         "a": 0.0615,
    #         "a_up": 0.0615 * 1.5,
    #         "a_down": 0.0615 * 0.5,
    #         "b": -0.0005,
    #         "b_up": -0.0005 * 1.5,
    #         "b_down": -0.0005 * 0.5,
    #     },
    #     pt_max=500.0,
    # )


   ################################################################################################
    # luminosity and normalization
    ################################################################################################

    # lumi values in 1/pb (= 1000/fb)
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3?rev=36
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis
    # difference pre-post VFP: https://cds.cern.ch/record/2854610/files/DP2023_006.pdf
    # Lumis for Run3 within the Twiki are outdated as stated here:
    # https://cms-talk.web.cern.ch/t/luminosity-in-run2023c/116859/2
    # Run3 Lumis can be calculated with brilcalc tool https://twiki.cern.ch/twiki/bin/view/CMS/BrilcalcQuickStart?rev=15
    # CClub computed this already: https://gitlab.cern.ch/cclubbtautau/AnalysisCore/-/issues/49

    if year == 2022 and campaign.x.tag == "preEE":
        cfg.x.luminosity = Number(7_980.4541, {
            "lumi_13p6TeV_2022": 0.014j,
            "lumi_13p6TeV_1": 0.0138j,
        })
    elif year == 2022 and campaign.x.tag == "postEE":
        cfg.x.luminosity = Number(26_671.6097, {
            "lumi_13p6TeV_2022": 0.014j,
            "lumi_13p6TeV_1": 0.0138j,
        })
    elif year == 2023 and campaign.x.tag == "preBPix":
        cfg.x.luminosity = Number(18_062.6591, {
            "lumi_13p6TeV_2023": 0.013j,
            "lumi_13p6TeV_1": 0.0017j,
            "lumi_13p6TeV_2": 0.0127j,
        })
    elif year == 2023 and campaign.x.tag == "postBPix":
        cfg.x.luminosity = Number(9_693.1301, {
            "lumi_13p6TeV_2023": 0.013j,
            "lumi_13p6TeV_1": 0.0017j,
            "lumi_13p6TeV_2": 0.0127j,
        })
    elif year == 2024:
        cfg.x.luminosity = Number(104_675.143180, {
            "lumi_13p6TeV_2024": 0.016j,
            "lumi_13p6TeV_1": 0.0020j,
            "lumi_13p6TeV_2": 0.0068j,
            "lumi_13p6TeV_3": 0.0144j,
        })
    else:
        assert False

    # minimum bias cross section in mb (milli) for creating PU weights, values from
    # https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData?rev=52#Recommended_cross_section
    cfg.x.minbias_xs = Number(69.2, 0.046j)

   
    # names of muon correction sets and working points
    # (used in the muon producer)   
  
    cfg.x.deep_tau = DotDict.wrap({
        "tagger": "DeepTau2018v2p5",
        "vs_e"          : {"mutau": "VVLoose",
                           "etau": "Tight",
                           "tautau": "VVLoose"},        
        "vs_mu"         : {"mutau": "Tight",
                           "etau": "VLoose",
                           "tautau": "VLoose"},
        "vs_jet"        : {"mutau": "Medium",
                           "etau": "Medium",
                           "tautau": "Medium"},
        "vs_e_jet_wps"  : {"VVVLoose"   : 1,
                           "VVLoose"    : 2,
                           "VLoose"     : 3,
                           "Loose"      : 4,
                           "Medium"     : 5,
                           "Tight"      : 6,
                           "VTight"     : 7,
                           "VVTight"    : 8},
        "vs_mu_wps"     : {"VLoose" : 1,
                           "Loose"  : 2,
                           "Medium" : 3,
                           "Tight"  : 4}
        })
################################################################################################
# b tagging
################################################################################################
    # name of the btag_sf correction set and jec uncertainties to propagate through
    cfg.x.btag_sf = ("deepJet_shape", cfg.x.btag_sf_jec_sources)

    cfg.x.btag_working_points = DotDict.wrap(
            {   2022 : {
                    "preEE":{
                        "deepjet" : { #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
                            "loose"  : 0.0583,
                            "medium" : 0.3086,
                            "tight"  : 0.7183,
                            "xtight" : 0.8111,
                            "xxtight": 0.9512,
                        },
                        "particleNet": {
                            "loose"  : 0.047,
                            "medium" : 0.245, 
                            "tight"  : 0.6734, 
                            "xtight" : 0.7862, 
                            "xxtight": 0.961, 
                        },
                        "robustParticleTransformer": {
                            "loose"  : 0.0849, 
                            "medium" : 0.4319, 
                            "tight"  : 0.8482, 
                            "xtight" : 0.9151, 
                            "xxtight": 0.9874,
                        },                    
                    },
                    "postEE":{
                        "deepjet" : { #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
                            "loose"  : 0.0614,
                            "medium" : 0.3196,
                            "tight"  : 0.73,
                            "xtight" : 0.8184,
                            "xxtight": 0.9542, 
                        },
                        "particleNet": {
                            "loose"  : 0.0499,
                            "medium" : 0.2605, 
                            "tight"  : 0.6915, 
                            "xtight" : 0.8033, 
                            "xxtight": 0.9664, 
                        },
                        "robustParticleTransformer": {
                            "loose"  : 0.0897, 
                            "medium" : 0.451, 
                            "tight"  : 0.8604, 
                            "xtight" : 0.9234, 
                            "xxtight": 0.9893,
                        },             
                        
                    },
                },
                2023 : {
                    "preBPix":{
                        "deepjet" : { #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
                            "loose"  : 0.0479,
                            "medium" : 0.2431,
                            "tight"  : 0.6553,
                            "xtight" : 0.7667, 
                            "xxtight": 0.9459, 
                        },
                        "particleNet": {
                            "loose"  : 0.0358, 
                            "medium" : 0.1917,
                            "tight"  : 0.6172, 
                            "xtight" : 0.7515, 
                            "xxtight": 0.9659, 
                        },
                        "robustParticleTransformer": {
                            "loose"  : 0.0681,
                            "medium" : 0.3487, 
                            "tight"  : 0.7969, 
                            "xtight" : 0.8882,
                            "xxtight": 0.9883,
                        },             
                    },
                    "postBPix":{
                        "deepjet" : { #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
                            "loose"  : 0.048,
                            "medium" : 0.2435,
                            "tight"  : 0.6563,
                            "xtight" : 0.7671,
                            "xxtight": 0.9483,
                        },
                        "particleNet": {
                            "loose"  : 0.0359,
                            "medium" : 0.1919,
                            "tight"  : 0.6133,
                            "xtight" : 0.7544,
                            "xxtight": 0.9688,
                        },
                        "robustParticleTransformer": {
                            "loose"  : 0.0683,
                            "medium" : 0.3494,
                            "tight"  : 0.7994,
                            "xtight" : 0.8877,
                            "xxtight": 0.9883,
                        },             
                    },
                }
            }
        )
    # JEC uncertainty sources propagated to btag scale factors
    # (names derived from contents in BTV correctionlib file)
    cfg.x.btag_sf_jec_sources = [
        "",  # same as "Total"
        "Absolute",
        "AbsoluteMPFBias",
        "AbsoluteScale",
        "AbsoluteStat",
        f"Absolute_{year}",
        "BBEC1",
        f"BBEC1_{year}",
        "EC2",
        f"EC2_{year}",
        "FlavorQCD",
        "Fragmentation",
        "HF",
        f"HF_{year}",
        "PileUpDataMC",
        "PileUpPtBB",
        "PileUpPtEC1",
        "PileUpPtEC2",
        "PileUpPtHF",
        "PileUpPtRef",
        "RelativeBal",
        "RelativeFSR",
        "RelativeJEREC1",
        "RelativeJEREC2",
        "RelativeJERHF",
        "RelativePtBB",
        "RelativePtEC1",
        "RelativePtEC2",
        "RelativePtHF",
        "RelativeSample",
        f"RelativeSample_{year}",
        "RelativeStatEC",
        "RelativeStatFSR",
        "RelativeStatHF",
        "SinglePionECAL",
        "SinglePionHCAL",
        "TimePtEta",
    ]

    from columnflow.production.cms.btag import BTagSFConfig
    cfg.x.btag_sf_deepjet = BTagSFConfig(
        correction_set="deepJet_shape",
        jec_sources=cfg.x.btag_sf_jec_sources,
        discriminator="btagDeepFlavB",
    )
    if year >= 2022:
        cfg.x.btag_sf_pnet = BTagSFConfig(
            correction_set="particleNet_shape",
            jec_sources=cfg.x.btag_sf_jec_sources,
            discriminator="btagPNetB",
        )   
    ################################################################################################
    # json file paths
    ################################################################################################       
     
    jsonpog_dir = "/eos/user/a/anigamov/htt_corrections_mirror/jsonpog-integration_latest/POG/"
    jsonpog_tau_dir = "/eos/user/a/anigamov/htt_corrections_mirror/jsonpog-integration_tau_latest/POG"
    corr_dir = "/eos/user/a/anigamov/htt_corrections_mirror/"
    # stiching_eos_path = "/eos/user/j/jmalvaso/SWAN_projects/higgs_MSSM/"
    #CMS Analysis Corrections Documentation: https://cms-analysis-corrections.docs.cern.ch/
    json_acd_path="/cvmfs/cms-griddata.cern.ch/cat/metadata/" 
    sz_path = "/afs/cern.ch/user/d/dmroy/public/ZpT_RecCorr_V5"
    golden_ls = { 
        2022 : "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json", 
        2023 : "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json"
    }    

    cfg.x.external_files = DotDict.wrap({
        "lumi": {
            "golden": (golden_ls[year], "v1"),
            "normtag": ("/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json", "v1"), #/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags
        }, 
        "pu_sf"                    : (f"{json_acd_path}LUM/{json_acd_tag}/latest/puWeights.json.gz", "v2"),
        "muon_correction"          : f"{json_acd_path}MUO/{json_acd_tag}/latest/muon_Z.json.gz",
        "HLT_mu_eff"               : f"{corr_dir}hleprare/TriggerScaleFactors/{cfg.x.year}{campaign.x.tag}/MuHlt_abseta_pt_wEff.json",
        "electron_scaling_smearing": f"{json_acd_path}EGM/{json_acd_tag}/latest/electronSS_EtDependent.json.gz",
        "electron_idiso"           : f"{json_acd_path}EGM/{json_acd_tag}/latest/electron.json.gz",
        "electron_trigger"         : f"{json_acd_path}EGM/{json_acd_tag}/latest/electronHlt.json.gz",
        "tau_correction"           : f"{json_acd_path}TAU/{json_acd_tag}/latest/tau.json.gz", #tau_DeepTau2018v2p5_{cfg.x.year}_{tau_tag}
        "zpt_weight"              : (f"{corr_dir}dy_ptll/DY_pTll_weights_{cfg.x.year}{campaign.x.tag}.json.gz","v2"),
        #"jet_jerc"                 : (f"{jsonpog_dir}JME/{cfg.x.year}_{tag}/jet_jerc.json.gz", "v2"),
        "jet_jerc"                 : (f"{json_acd_path}JME/{json_acd_tag}/latest/jet_jerc.json.gz", "v2"),
        "jet_veto_map"             : (f"{json_acd_path}JME/{json_acd_tag}/latest/jetvetomaps.json.gz", "v2"),
        "btag_sf_corr"             : (f"{json_acd_path}BTV/{json_acd_tag}/latest/btagging.json.gz", "v2"),
        "met_recoil"               : (f"{corr_dir}dy_ptll/Recoil_corrections_{cfg.x.year}{campaign.x.tag}.json.gz", "v2"),
    })
    #/eos/user/a/anigamov/htt_corrections_mirror/dy_ptll'
    # --------------------------------------------------------------------------------------------- #
    # electron settings
    # names of electron correction sets and working points
    # (used in the electron_sf producer)
    # --------------------------------------------------------------------------------------------- #
    cfg.x.electron_sf = DotDict.wrap({
        "ID": {"corrector": "Electron-ID-SF",
               "year": e_sf_tag,
               "wp":"wp90noiso"},
        "scale": {"corrector": e_scale_corrector},
        "smearing": {"corrector": e_smearing_corrector},
        "trig": {"corrector": "Electron-HLT-SF",
                 "year": e_sf_tag,
                 "wp": "HLT_SF_Ele30_TightID"},
        "xtrig": {"corrector": "Electron-HLT-SF",
                  "year": e_sf_tag,
                  "wp": "HLT_SF_Ele24_TightID"},
        "MC_eff": {"corrector": "Electron-HLT-McEff",
                  "year": e_sf_tag,
                  "wp": "HLT_SF_Ele30_TightID"},
        "Data_eff": {"corrector": "Electron-HLT-DataEff",
                     "year": e_sf_tag,
                     "wp": "HLT_SF_Ele30_TightID"},
    })
    
    # --------------------------------------------------------------------------------------------- #
    # muon settings
    # names of muon correction sets and working points
    # (used in the muon producer)
    # --------------------------------------------------------------------------------------------- #
    cfg.x.muon_id_wp = 'Medium'
    cfg.x.muon_sf = DotDict.wrap({ 
                                  
       'ID': {'corrector': f"NUM_{cfg.x.muon_id_wp}ID_DEN_TrackerMuons",},
        
        'iso': {'corrector': f"NUM_TightPFIso_DEN_{cfg.x.muon_id_wp}ID",},
        
        'trig': {'corrector': f"NUM_IsoMu24_DEN_CutBasedId{cfg.x.muon_id_wp}_and_PFIsoMedium",},
                 #TODO Worging points for trig and iso does not match!
                
        "trig_mc_eff": {"corrector": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_MCeff",
                 "year": f"{year}_{tag}"},
        
        "trig_data_eff": {"corrector": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_DATAeff",
                 "year": f"{year}_{tag}"},
        
        'xtrig': {'corrector': f"NUM_IsoMu20_DEN_CutBasedId{cfg.x.muon_id_wp}_and_PFIsoMedium",},
        
        "MC_eff_mutau": {"corrector": "NUM_IsoMu20_DEN_CutBasedIdTight_and_PFIsoTight_MCeff"},
        
        "Data_eff_mutau": {"corrector": "NUM_IsoMu20_DEN_CutBasedIdTight_and_PFIsoTight_DATAeff"},
    })
    
    # --------------------------------------------------------------------------------------------- #
    # stitching settings
    # names of tau correction sets and working points
    # (used in the stitching weights producer)
    # --------------------------------------------------------------------------------------------- #
    cfg.x.stitching_settings = DotDict.wrap({
        "corrector": "stitch_weight",   # must match the name in your correction file
        "era": f"{campaign.x.year}{campaign.x.tag}",     
    })
    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    
    from MSSM_H_tt.config.variables import keep_columns
    keep_columns(cfg)
 
    cfg.add_shift(name="nominal", id=0)

    cfg.add_shift(name="tau_weight_down", id=1, type="shape")
    cfg.add_shift(name="tau_weight_up", id=2, type="shape")
    add_shift_aliases(cfg, "tau_weight", {"tau_weight": "tau_weight_{direction}"})
    
    cfg.add_shift(name="muon_weight_down", id=3, type="shape")
    cfg.add_shift(name="muon_weight_up", id=4, type="shape")
    add_shift_aliases(cfg, "muon_weight", {"muon_weight": "muon_weight_{direction}"})

    cfg.add_shift(name="electron_weight_down", id=8, type="shape")
    cfg.add_shift(name="electron_weight_up", id=9, type="shape")
    add_shift_aliases(cfg, "electron_weight", {"electron_weight": "electron_weight_{direction}"})
    
    cfg.add_shift(name="top_pt_weight_down", id=10, type="shape")
    cfg.add_shift(name="top_pt_weight_up", id=11, type="shape")
    add_shift_aliases(cfg, "top_pt_weight", {"top_pt_weight": "top_pt_weight_{direction}"})
    
    cfg.add_shift(name="Trigger_SF_weight_down", id=12, type="shape")
    cfg.add_shift(name="Trigger_SF_weight_up", id=13, type="shape")
    add_shift_aliases(cfg, "Trigger_SF_weight", {"Trigger_SF_weight": "Trigger_SF_weight_{direction}"})
    
    cfg.add_shift(name="zpt_weight_down", id=14, type="shape")
    cfg.add_shift(name="zpt_weight_up", id=15, type="shape")
    add_shift_aliases(cfg, "zpt_weight", {"zpt_weight": "zpt_weight_{direction}"})
    
    cfg.add_shift(name="pu_weight_down", id=16, type="shape")
    cfg.add_shift(name="pu_weight_up", id=17, type="shape")
    add_shift_aliases(cfg,"pu_weight",{"pu_weight": "pu_weight_{direction}"})
    # add column aliases for shift jec
    for i, (jec_source, flag) in enumerate(all_jec_sources.items()):
        if not flag:
            continue
        cfg.add_shift(
            name=f"jec_{jec_source}_up",
            id=5000 + 2 * i,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        cfg.add_shift(
            name=f"jec_{jec_source}_down",
            id=5001 + 2 * i,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        add_shift_aliases(
            cfg,
            f"jec_{jec_source}",
            {
                "Jet.pt": "Jet.pt_{name}",
                "Jet.mass": "Jet.mass_{name}",
                f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
                f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
            },
        )
        # # TODO: check the JEC de/correlation across years and the interplay with btag weights
        # if ("" if jec_source == "Total" else jec_source) in cfg.x.btag_sf_jec_sources:
        #     add_shift_aliases(
        #         cfg,
        #         f"jec_{jec_source}",
        #         {
        #             "normalized_btag_weight_pnet": "normalized_btag_weight_pnet_{name}",
        #             "normalized_njet_btag_weight_pnet": "normalized_njet_btag_weight_pnet_{name}",
        #         },
        #     )

    cfg.add_shift(name="jer_up", id=6000, type="shape", tags={"jer"})
    cfg.add_shift(name="jer_down", id=6001, type="shape", tags={"jer"})
    add_shift_aliases(
        cfg,
        "jer",
        {
            "Jet.pt": "Jet.pt_{name}",
            "Jet.mass": "Jet.mass_{name}",
            f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
            f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
        },
    )

    # cfg.add_shift(name="btag_weight_SF_up", id=15, type="shape")
    # cfg.add_shift(name="btag_weight_SF_down", id=16, type="shape")
    # add_shift_aliases(cfg, "btag_weight_SF", {"btag_weight_SF_nom": "btag_weight_SF_nom_{direction}"})
    
    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    get_shifts = functools.partial(get_shifts_from_sources, cfg)   
    
    cfg.x.event_weights = DotDict({
        "normalization_weight": [],
        # "filter_weight": [],
        # "mc_weight":[],
        "tau_weight": get_shifts("tau_weight"),
        "pu_weight": get_shifts("pu_weight"),
        "zpt_weight":get_shifts("zpt_weight"),
        "muon_weight": get_shifts("muon_weight"),
        "electron_weight": get_shifts("electron_weight"), 
        "top_pt_weight" : get_shifts("top_pt_weight"),     
        "pu_weight": get_shifts("pu_weight"),
        "btag_weight_SF_nom": [],
        "Trigger_SF_weight": get_shifts("Trigger_SF_weight"),
        "stitching_weight": [],
    })

    # thisdir = os.path.dirname(os.path.abspath(__file__))
    
    # with open(os.path.join(thisdir, "jec_sources.yaml"), "r") as f:
    #     all_jec_sources = yaml.load(f, yaml.Loader)["names"]

    # for jec_source in cfg.x.jec["uncertainty_sources"]:
    #     idx = all_jec_sources.index(jec_source)
    #     cfg.add_shift(
    #         name=f"jec_{jec_source}_up",
    #         id=5000 + 2 * idx,
    #         type="shape",
    #         tags={"jec"},
    #         aux={"jec_source": jec_source},
    #     )
    #     cfg.add_shift(
    #         name=f"jec_{jec_source}_down",
    #         id=5001 + 2 * idx,
    #         type="shape",
    #         tags={"jec"},
    #         aux={"jec_source": jec_source},
    #     )
    #     add_shift_aliases(
    #         cfg,
    #         f"jec_{jec_source}",
    #         {"Jet.pt": "Jet.pt_{name}", "Jet.mass": "Jet.mass_{name}"},
    #     )

    #     if jec_source in ["Total", *cfg.x.btag_sf_jec_sources]:
    #         # when jec_source is a known btag SF source, add aliases for btag weight column
    #         add_shift_aliases(
    #             cfg,
    #             f"jec_{jec_source}",
    #             {
    #                 "btag_weight": f"btag_weight_jec_{jec_source}_" + "{direction}",
    #                 "normalized_btag_weight": f"normalized_btag_weight_jec_{jec_source}_" + "{direction}",
    #                 "normalized_njet_btag_weight": f"normalized_njet_btag_weight_jec_{jec_source}_" + "{direction}",
    #             },
    #         )

    # cfg.add_shift(name="jer_up", id=6000, type="shape", tags={"jer"})
    # cfg.add_shift(name="jer_down", id=6001, type="shape", tags={"jer"})
    # add_shift_aliases(cfg, "jer", {"Jet.pt": "Jet.pt_{name}", "Jet.mass": "Jet.mass_{name}"})
  
    

    
    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    def set_version(cls, inst, params):
        # per default, use the version set on the command line
        version = inst.version 
        return version if version else "dev"
            
        
    cfg.x.versions = {
        "cf.CalibrateEvents"    : set_version,
        "cf.SelectEvents"       : set_version,
        "cf.MergeSelectionStats": set_version,
        "cf.MergeSelectionMasks": set_version,
        "cf.ReduceEvents"       : set_version,
        "cf.MergeReductionStats": set_version,
        "cf.MergeReducedEvents" : set_version,
    }
    # channels
    # processing only one channel at once, electron calibration is channel dependent!
    channel_id = {"etau"  : 1,
                  "mutau" : 2,
                  "emu"   : 3,
                  "tautau": 4}
    cfg.add_channel(name=channel,   id=channel_id[channel])
    
    cfg.x.ch_objects = DotDict.wrap({
        "etau"   : {"lep0" : "Electron",
                    "lep1" : "Tau"    },
        "mutau"  : {"lep0" : "Muon",
                    "lep1" : "Tau"    },
        "emu"    : {"lep0" : "Electron",
                    "lep1" : "Muon"   },
        "tautau" : {"lep0" : "Tau",
                    "lep1" : "Tau"    },
    })

    cfg.x.met_recoil = DotDict.wrap({
        "datasets" : {},
    })
    cfg.x.met_recoil["datasets"].update({
                                # DY -> ll
                                "DYto2L_M_50_0J_amcatnloFXFX": "NLO",
                                "DYto2L_M_50_1J_amcatnloFXFX": "NLO",
                                "DYto2L_M_50_2J_amcatnloFXFX": "NLO",
                                "DYto2L_M_50_amcatnloFXFX": "NLO",
                                " DYto2L_M_10to50_amcatnloFXFX":"NLO",
                                # DY -> tautau (POWHEG + aMC@NLO)
                                "DYto2Tau_MLL_50_0J_amcatnloFXFX": "NLO",
                                "DYto2Tau_MLL_50_1J_amcatnloFXFX": "NLO",
                                "DYto2Tau_MLL_50_2J_amcatnloFXFX": "NLO",
                                # W + jets (MG5 MLM ~ LO)
                                "WtoLNu_1J_madgraphMLM": "LO",
                                "WtoLNu_2J_madgraphMLM": "LO",
                                "WtoLNu_3J_madgraphMLM": "LO",
                                "WtoLNu_4J_madgraphMLM": "LO",
                                "WtoLNu_madgraphMLM": "LO",
                            })
    stitch_samples = [
            "DYto2L_M_50_amcatnloFXFX",
            "DYto2L_M_50_0J_amcatnloFXFX",
            "DYto2L_M_50_1J_amcatnloFXFX",
            "DYto2L_M_50_2J_amcatnloFXFX",
            "WtoLNu_1J_madgraphMLM",
            "WtoLNu_2J_madgraphMLM",
            "WtoLNu_3J_madgraphMLM",
            "WtoLNu_4J_madgraphMLM",
            "WtoLNu_madgraphMLM",
        ]

    cfg.x.stitch_samples = stitch_samples
    
    bugged_DYto2Tau_samples = [
            "DYto2L_M_10to50_amcatnloFXFX",
            "DYto2L_M_50_amcatnloFXFX",
            "DYto2L_M_50_0J_amcatnloFXFX",
            "DYto2L_M_50_1J_amcatnloFXFX",
            "DYto2L_M_50_2J_amcatnloFXFX",
        ]

    cfg.x.bugged_DYto2Tau_samples = bugged_DYto2Tau_samples
    
    # cfg.x.fake_factor_method = DotDict.wrap({
    # "axes": {"delta_r" : {
    #             "var_route": [f"hcand_{channel}","delta_r"],
    #             "ax_str"  : "Variable([0.3,3,3.5,4,6], name='delta_r', label='Delta R', underflow=False, overflow=False)",
    #             },
    #          "N_jets"  : {
    #             "var_route" : ["n_jets"],
    #             "ax_str"  : "Integer(0, 3, name='N_jets', label='Number of jets', underflow=False, overflow=False)",
    #             },
    #          "N_b_jets": {
    #             "var_route" : ["N_b_jets"],
    #             "ax_str"   : "Integer(0, 2, name='N_b_jets', label='Number of b jets',underflow=False, overflow=False)",
    #             },
    # },
    # "columns" : ["ff_weight_qcd"],
    # "shifts"  : ["up", "nominal", "down"]
    # })   
    
    # shift groups for conveniently looping over certain shifts
    # (used during plotting)
    cfg.x.shift_groups = {
        "jec": [
            shift_inst.name for shift_inst in cfg.shifts
            if shift_inst.has_tag(("jec", "jer"))
        ],
    }
    
    if cfg.campaign.x("custom").get("creator") == "desy":  
        def get_dataset_lfns(dataset_inst: od.Dataset, shift_inst: od.Shift, dataset_key: str) -> list[str]:
            # destructure dataset_key into parts and create the lfn base directory
            print(f"Creating custom get_dataset_lfns for {config_name}")   
            try:
               basepath = cfg.campaign.x("custom").get("location")
            except:
                print("Did not find any basebath in the campaigns")
                basepath = "" 
            lfn_base = law.wlcg.WLCGDirectoryTarget(
                f"{basepath}{dataset_key}",
                fs="wlcg_fs_eos",
            )
            print(f"lfn basedir:{lfn_base}")
            # loop though files and interpret paths as lfns
            return [
                lfn_base.child(basename, type="f").path
                for basename in lfn_base.listdir(pattern="*.root")
            ]
        # define the lfn retrieval function
        cfg.x.get_dataset_lfns = get_dataset_lfns
        # define a custom sandbox
        cfg.x.get_dataset_lfns_sandbox = dev_sandbox("bash::$CF_BASE/sandboxes/cf.sh")
        # define custom remote fs"s to look at
        cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: "wlcg_fs_eos"
        
    # add categories using the "add_category" tool which adds auto-generated ids
    from MSSM_H_tt.config.categories import add_categories
    add_categories(cfg,channel=channel)
        
    from MSSM_H_tt.config.variables import add_variables
    add_variables(cfg)
    
    from data_driven.hist_hooks import add_hist_hooks
    add_hist_hooks(ana)
        