
# coding: utf-8

"""
Configuration of the higgs_cp analysis.
"""

import functools
import itertools
import yaml
import law
import order as od
import os
from scinum import Number

from columnflow.util import maybe_import, dev_sandbox
from law.util import DotDict
from columnflow.config_util import (
    get_root_processes_from_campaign, 
    add_category, add_shift_aliases,
    verify_config_processes,get_shifts_from_sources
)
from columnflow.calibration.cms.tau import TECConfig
from columnflow.calibration.cms.egamma import EGammaCorrectionConfig

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

    # def tag_caster(campaign: od.Campaign) -> str:
    #     #Helper function to cast campaign tags to the tags used in POG groups for the scale factors
    #     year = campaign.x.year
    #     tag = campaign.x.tag
    #     out_tag = ''
    #     e_sf_tag = ''
    #     e_scale_corrector = ''
    #     e_smearing_corrector = '' 
    #     if year in [2017,2018]  : out_tag = 'UL'
    #     elif tag == "preEE"     : 
    #         out_tag = "Summer22"
    #         e_sf_tag = "2022Re-recoBCD"
    #         e_scale_corrector = "2022Re-recoBCD_ScaleJSON"
    #         e_smearing_corrector = "2022Re-recoBCD_SmearingJSON"
    #     elif tag == "postEE"    : 
    #         out_tag = "Summer22EE"
    #         e_sf_tag = "2022Re-recoE+PromptFG"
    #         e_scale_corrector = "2022Re-recoE+PromptFG_ScaleJSON"
    #         e_smearing_corrector = "2022Re-recoE+PromptFG_SmearingJSON"
    #     elif tag == "preBPix"   : 
    #         out_tag = "Summer23"
    #         e_sf_tag = "2023PromptC"
    #         e_scale_corrector = "2022Re-recoE+PromptFG_ScaleJSON"
    #         e_smearing_corrector = "2022Re-recoE+PromptFG_SmearingJSON"
    #     elif tag == "postBPix"  : 
    #         out_tag = "Summer23BPix"
    #         e_sf_tag = "2023PromptD"
    #         e_scale_corrector = "2022Re-recoE+PromptFG_ScaleJSON"
    #         e_smearing_corrector = "2022Re-recoE+PromptFG_SmearingJSON"
    #     elif tag == "preVFP"    : out_tag = "preVFP_UL"
    #     elif tag == "postVFP"   : out_tag = "postVFP_UL"    
    #     return out_tag, e_sf_tag, e_scale_corrector, e_smearing_corrector, tag
        
    # tag, e_sf_tag, e_scale_corrector, e_smearing_corrector, tau_tag = tag_caster(campaign)
    
    tag_dict = {
        "preEE"     : {
            "short_tag"             : "",
            "long_tag"              : "preEE",
            "pog_tag"               : "Summer22",
            "e_sf_tag"              : "2022Re-recoBCD",
            "e_scale_corrector"     : "2022Re-recoBCD_ScaleJSON",
            "e_smearing_corrector"  : "2022Re-recoBCD_SmearingJSON",
            "jerc_postfix"          : "",
            "cat_tag"               : "Run3-22CDSep23-Summer22-NanoAODv12",
            },
        "postEE"    : {
            "short_tag"             : "EE",
            "long_tag"              : "postEE",
            "pog_tag"               : "Summer22EE",
            "e_sf_tag"              : "2022Re-recoE+PromptFG",
            "e_scale_corrector"     : "2022Re-recoE+PromptFG_ScaleJSON",
            "e_smearing_corrector"  : "2022Re-recoE+PromptFG_SmearingJSON",
            "jerc_postfix"          : "EE",
            "cat_tag"               : "Run3-22EFGSep23-Summer22EE-NanoAODv12",
            },
        "preBPix"   : {
            "short_tag"             : "",
            "long_tag"              : "preBPix",
            "pog_tag"               : "Summer23",
            "e_sf_tag"              : "2023PromptC",
            "e_scale_corrector"     : "2022Re-recoE+PromptFG_ScaleJSON",
            "e_smearing_corrector"  : "2022Re-recoE+PromptFG_SmearingJSON",
            "jerc_postfix"          : "",
            "cat_tag"               : "Run3-23CSep23-Summer23-NanoAODv12",
            },
        "postBPix"  : {
            "short_tag"             : "BPix",
            "long_tag"              : "postBPix",
            "pog_tag"               : "Summer23BPix",
            "e_sf_tag"              : "2023PromptD",
            "e_scale_corrector"     : "2022Re-recoE+PromptFG_ScaleJSON",
            "e_smearing_corrector"  : "2022Re-recoE+PromptFG_SmearingJSON",
            "jerc_postfix"          : "BPix",
            "cat_tag"               : "Run3-23DSep23-Summer23BPix-NanoAODv12",
            },
        ""  : { #Default values in case tag is empty
            "short_tag"             : "",
            "long_tag"              : "",
            "pog_tag"               : "",
            "e_sf_tag"              : "",
            "e_scale_corrector"     : "",
            "e_smearing_corrector"  : "",
            "jerc_postfix"          : "",
            },
    }
    
    tags = tag_dict[campaign.x.tag]
    tag = tags['long_tag']
    
    # add processes we are interested in
    process_names = [
        "data", 
        "data_mu",
        "data_tau",
        "data_e",
        "data_singlemu",
        #Drell-Yan
        "dy_lep",
        #"dy_z2ee",
        #"dy_z2mumu",
        #"dy_z2tautau",
        "dy_ll_m10to50",
        "dy_ll_m50",
        "dy_tt_m50",
        # "dy_lep_m10to50",
        #W + jets
        "wj",
        "wj_1j",
        "wj_2j",
        "wj_3j",
        "wj_4j",
        #diboson + single top
        "vvt",
        #diboson
        "vv", #diboson inclusive
        "ww",
        "wz",
        "zz",
        #ttbar
        "tt",#ttbar inclusive
        "tt_sl",
        "tt_dl",
        "tt_fh",
        #single top
        "st",
        #single top t-channel        
        "st_tchannel_tbar",
        "st_tchannel_t",
        #single top s-channel   
        "st_schannel_t_lep",
        "st_schannel_tbar_lep",
        # single top tW channel
        "st_twchannel_t_fh","st_twchannel_t_sl","st_twchannel_t_dl",
        "st_twchannel_tbar_sl","st_twchannel_tbar_dl","st_twchannel_tbar_fh",
        #ggF signal
        "h_ggf_htt_sm",
        "h_ggf_htt_sm_prod_sm","h_ggf_htt_sm_prod_mm","h_ggf_htt_sm_prod_cpo",
        "h_ggf_htt_mm",
        "h_ggf_htt_mm_prod_sm","h_ggf_htt_mm_prod_mm","h_ggf_htt_mm_prod_cpo",
        "h_ggf_htt_cpo",
        "h_ggf_htt_cpo_prod_sm","h_ggf_htt_cpo_prod_mm","h_ggf_htt_cpo_prod_cpo",
        "h_ggf_htt_flat",
        "h_ggf_htt_flat_prod_sm","h_ggf_htt_flat_prod_mm","h_ggf_htt_flat_prod_cpo", 
        #VBF signal
        "h_vbf_htt_cpo","h_vbf_htt_sm","h_vbf_htt_mm","h_vbf_htt_flat",
        #VH signal
        "zh_htt_cpo","zh_htt_sm","zh_htt_mm","zh_htt_flat",
        "wh_htt_cpo","wh_htt_sm","wh_htt_mm","wh_htt_flat",
        "wph_htt_cpo","wph_htt_sm","wph_htt_mm","wph_htt_flat",
        "wmh_htt_cpo","wmh_htt_sm","wmh_htt_mm","wmh_htt_flat",
    ]
    for process_name in process_names:
        # add the process
        #print(f'importing {process_name} for {campaign.name}')
        proc = cfg.add_process(procs.get(process_name))
        #for signal datasets create special tag
        if process_name.startswith("h_"):
            proc.add_tag("signal")
           
    # add datasets we need to study
    from httcp.config.datasets import add_datasets_2025_skim_v2,add_datasets_2024_skim_v1
    datasets = add_datasets_2025_skim_v2()
    #Simultaneously select the dataset list depending on the year,campaign tag, and iterate over the constituents
    for dataset_name in datasets[f"{year}{tag}"]:
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))
        if dataset_name.startswith("h_") or dataset_name.startswith("zh_") or dataset_name.startswith("wh_"):
            dataset.add_tag("signal")   
        if dataset.name.startswith("tt_"):
            dataset.add_tag({"has_top", "ttbar", "tt"})    
        # for testing purposes, limit the number of files to 1
        for info in dataset.info.values():
            if limit_dataset_files:
                info.n_files = min(info.n_files, limit_dataset_files) #<<< REMOVE THIS FOR THE FULL DATASET
    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)
    
    cfg.x.stitch_samples = cfg.x.datasets2apply_tau_veto = [
        'DYto2L_M_50_amcatnloFXFX',
        'DYto2L_M_50_0J_amcatnloFXFX',
        'DYto2L_M_50_1J_amcatnloFXFX',
        'DYto2L_M_50_2J_amcatnloFXFX',
        "WtoLNu_1J_madgraphMLM",
        "WtoLNu_2J_madgraphMLM",
        "WtoLNu_3J_madgraphMLM",
        "WtoLNu_4J_madgraphMLM",
        "WtoLNu_madgraphMLM",
    ]

    #Adding the triggers 
    from httcp.config.triggers import add_triggers_run3
    add_triggers_run3(cfg)

    from httcp.config.met_filters import add_met_filters
    add_met_filters(cfg)

    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    cfg.x.default_calibrator = "main"
    cfg.x.default_selector = "main"
    cfg.x.default_reducer = "cf_default"
    cfg.x.default_producer = "main"
    cfg.x.default_hist_producer = "httcp_hist_producer"
    cfg.x.default_ml_model = None
    cfg.x.default_inference_model = "hcp_model"
    cfg.x.default_categories = ("incl",)
    cfg.x.default_variables = ("event")

    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {
        "data" : ["data_mu", "data_tau","data_e"],
        "vv"   : ["ww", "wz", "zz"],
        "tt"   : ["tt_sl","tt_dl","tt_fh"],
        "st"   : ["st_tchannel_tbar","st_tchannel_t","st_schannel_tbar_lep","st_schannel_t_lep",
               "st_twchannel_t_fh","st_twchannel_t_sl","st_twchannel_t_dl",
               "st_twchannel_tbar_sl","st_twchannel_tbar_dl","st_twchannel_tbar_fh",],
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

    # shift groups for conveniently looping over certain shifts
    # (used during plotting)
    cfg.x.shift_groups = {}

    # selector step groups for conveniently looping over certain steps
    # (used in cutflow tasks)
    cfg.x.selector_step_groups = {
        "default": ["json", "met_filter", "dl_res_veto", "trigger", "lepton", "jet"],
    }
     
    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False

    jet_type = "AK4PFPuppi"
    if year < 2022:
        jerc_campaign = f"Summer19UL{year2}{tags['jerc_postfix']}"
        jet_type = "AK4PFchs"
    elif year==2022:
        jerc_campaign = f"Summer{year2}{tags['jerc_postfix']}_22Sep2023"
    elif year==2023:
        jerc_campaign = f"Summer{year2}{tags['jerc_postfix']}Prompt23"
        
    cfg.x.jec = DotDict.wrap({
        "campaign": jerc_campaign,
        "version": {2016: "V7", 2017: "V5", 2018: "V5", 2022: "V2", 2023:"V1"}[year],
        "jet_type": jet_type,
        "levels": ["L1L2L3Res"], #"L2Relative", "L2L3Residual", "L3Absolute", "L1L2L3Res" 
        "levels_DATA": ["L1L2L3Res"], #"L2Relative", "L2L3Residual", "L3Absolute", "L1L2L3Res" 
        "levels_MC": ["L1L2L3Res"], 
        "levels_for_type1_met": ["L1L2L3Res"], 
        "uncertainty_sources": [
            # "AbsoluteStat",
            # "AbsoluteScale",
            # "AbsoluteSample",
            # "AbsoluteFlavMap",
            # "AbsoluteMPFBias",
            # "Fragmentation",
            # "SinglePionECAL",
            # "SinglePionHCAL",
            # "FlavorQCD",
            # "TimePtEta",
            # "RelativeJEREC1",
            # "RelativeJEREC2",
            # "RelativeJERHF",
            # "RelativePtBB",
            # "RelativePtEC1",
            # "RelativePtEC2",
            # "RelativePtHF",
            # "RelativeBal",
            # "RelativeSample",
            # "RelativeFSR",
            # "RelativeStatFSR",
            # "RelativeStatEC",
            # "RelativeStatHF",
            # "PileUpDataMC",
            # "PileUpPtRef",
            # "PileUpPtBB",
            # "PileUpPtEC1",
            # "PileUpPtEC2",
            # "PileUpPtHF",
            # "PileUpMuZero",
            # "PileUpEnvelope",
            # "SubTotalPileUp",
            # "SubTotalRelative",
            # "SubTotalPt",
            # "SubTotalScale",
            # "SubTotalAbsolute",
            # "SubTotalMC",
            "Total",
            # "TotalNoFlavor",
            # "TotalNoTime",
            # "TotalNoFlavorNoTime",
            # "FlavorZJet",
            # "FlavorPhotonJet",
            # "FlavorPureGluon",
            # "FlavorPureQuark",
            # "FlavorPureCharm",
            # "FlavorPureBottom",
            # # "TimeRunA",
            # # "TimeRunB",
            # # "TimeRunC",
            # # "TimeRunD",
            # "CorrelationGroupMPFInSitu",
            # "CorrelationGroupIntercalibration",
            # "CorrelationGroupbJES",
            # "CorrelationGroupFlavor",
            # "CorrelationGroupUncorrelated",
        ],
    })
    
    ###############################################################################################
    # met settings
    ################################################################################################

    # name of the MET phi correction set
    # (used in the met_phi calibrator)
    cfg.x.met_name = 'PuppiMET'
    #cfg.x.met_phi_correction_set = r"{variable}_metphicorr_pfmet_{data_source}"
    
    ###############################################################################################
    # JER settings
    ################################################################################################
    
    # # JER
    # # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=107
    # # TODO: get jerc working for Run3
    # cfg.x.jer = DotDict.wrap({
    #     "campaign": jerc_campaign,
    #     "version": {2016: "JRV3", 2017: "JRV2", 2018: "JRV2", 2022: "JRV1"}[year],
    #     "jet_type": jet_type,
    # })

    # # JEC uncertainty sources propagated to btag scale factors
    # # (names derived from contents in BTV correctionlib file)
    # cfg.x.btag_sf_jec_sources = [
    #     "",  # total
    #     "Absolute",
    #     "AbsoluteMPFBias",
    #     "AbsoluteScale",
    #     "AbsoluteStat",
    #     f"Absolute_{year}",
    #     "BBEC1",
    #     f"BBEC1_{year}",
    #     "EC2",
    #     f"EC2_{year}",
    #     "FlavorQCD",
    #     "Fragmentation",
    #     "HF",
    #     f"HF_{year}",
    #     "PileUpDataMC",
    #     "PileUpPtBB",
    #     "PileUpPtEC1",
    #     "PileUpPtEC2",
    #     "PileUpPtHF",
    #     "PileUpPtRef",
    #     "RelativeBal",
    #     "RelativeFSR",
    #     "RelativeJEREC1",
    #     "RelativeJEREC2",
    #     "RelativeJERHF",
    #     "RelativePtBB",
    #     "RelativePtEC1",
    #     "RelativePtEC2",
    #     "RelativePtHF",
    #     "RelativeSample",
    #     f"RelativeSample_{year}",
    #     "RelativeStatEC",
    #     "RelativeStatFSR",
    #     "RelativeStatHF",
    #     "SinglePionECAL",
    #     "SinglePionHCAL",
    #     "TimePtEta",
    # ]
    
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

    ################################
    # luminosity and normalization #
    ################################

    # lumi values in inverse pb
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2?rev=2#Combination_and_correlations
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis
    # difference pre-post VFP: https://cds.cern.ch/record/2854610/files/DP2023_006.pdf
    
    
    lumi_dict = {
        "2022preEE"     : Number(7_980.4,  {"lumi_13p6TeV_correlated": 0.014j,}),
        "2022postEE"    : Number(26_671.7, {"lumi_13p6TeV_correlated": 0.014j,}),
        "2023preBPix"   : Number(17_794,   {"lumi_13p6TeV_correlated": 0.0j,}),
        "2023postBPix"  : Number(9_451,    {"lumi_13p6TeV_correlated": 0.0j,}),
        "2024"          : Number(109_080,  {"lumi_13p6TeV_correlated": 0.0j,}),
    }
    cfg.x.luminosity = lumi_dict[f"{year}{tag}"]
    
    cfg.x.deep_tau = DotDict.wrap({
        "tagger": "DeepTau2018v2p5",
        "vs_e"          : {"mutau": "VVLoose",
                           "etau": "Tight",
                           "tautau": "Tight"},        
        "vs_mu"         : {"mutau": "Tight",
                           "etau": "Tight",
                           "tautau": "Tight"},
        "vs_jet"        : {"mutau": "VTight",
                           "etau": "Medium",
                           "tautau": "Medium"},
        "vs_e_jet_wps"  : {'VVVLoose'   : 1,
                           'VVLoose'    : 2,
                           'VLoose'     : 3,
                           'Loose'      : 4,
                           'Medium'     : 5,
                           'Tight'      : 6,
                           'VTight'     : 7,
                           'VVTight'    : 8},
        "vs_mu_wps"     : {'VLoose' : 1,
                           'Loose'  : 2,
                           'Medium' : 3,
                           'Tight'  : 4}
        })
    cfg.x.btag_working_points = DotDict.wrap(
        {
            2016 : {
                "deepjet": { #TODO: make a link to this numbers
                    "loose": 0.0532,
                    "medium": 0.3040,
                    "tight": 0.7476,
                },
                "deepcsv": {
                    "loose": 0.1355,
                    "medium": 0.4506,
                    "tight": 0.7738,
                },
            },
            2022 : {
                "preEE":{
                    "deepjet" : { #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
                        "loose": 0.0583,
                        "medium": 0.3086,
                        "tight": 0.7183,
                    },
                },
                "postEE":{
                    "deepjet" : { #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
                        "loose": 0.0614,
                        "medium": 0.3196,
                        "tight": 0.73,
                    },
                },
            },
            2023 : {
                "preBPix":{
                    "deepjet" : { #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
                        "loose": 0.0479,
                        "medium": 0.2431,
                        "tight": 0.6553,
                    },
                },
                "postBPix":{
                    "deepjet" : { #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
                        "loose": 0.048,
                        "medium": 0.2435,
                        "tight": 0.6563,
                    },
                },

            }

                    
                
        },
    )
            
    cfg.x.tec = TECConfig(
            tagger=cfg.x.deep_tau.tagger,
            corrector_kwargs={"wp": getattr(cfg.x.deep_tau.vs_jet, channel), "wp_VSe": getattr(cfg.x.deep_tau.vs_e, channel)},
            )
    
    cfg.x.eec = EGammaCorrectionConfig(
                correction_set=f"EGMSmearAndSyst_ElePTsplit_{str(year)}{campaign.x.tag}",
                value_type="scale",
                uncertainty_type="escale",
                compound=True,
            )
    
    cfg.x.eer = EGammaCorrectionConfig(
                correction_set=f"EGMSmearAndSyst_ElePTsplit_{str(year)}{campaign.x.tag}",
                value_type="smear",
                uncertainty_type="esmear",
    )
    ##########################
    ###### mT cut value ######
    ##########################
    
    cfg.x.mt_cut_value = 65
    
    ##########################
    ##########################
    ##########################
    
    jsonpog_dir = "/eos/user/a/anigamov/htt_corrections_mirror/jsonpog-integration_latest/POG/"
    jsonpog_tau_dir = "//eos/user/a/anigamov/htt_corrections_mirror/jsonpog-integration_tau_latest/POG/"
    corr_dir = "/eos/user/a/anigamov/htt_corrections_mirror/"
    ml_dir = "/eos/user/s/stzakhar/TauTheDifference/Training/models/"
    golden_ls = { 
        2022 : "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json", 
        2023 : "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json"
    }  
    
    pog_tag = tags['pog_tag']
    short_tag = tags['short_tag']
    cfg.x.long_tag = long_tag = tags['long_tag']
    ch_short = channel.replace('mu','m').replace('tau','t')
    
    jsons_2024_v1 = DotDict.wrap({
        "lumi": {
            "golden": (golden_ls[year], "v1"),
            "normtag": ("/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json", "v1"), #/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags
        },

        "pu_sf"                         : (f"{jsonpog_dir}LUM/{year}_{pog_tag}/puWeights.json.gz", "v1"),
        "muon_correction"               : f"{jsonpog_dir}MUO/{year}_{pog_tag}/muon_Z.json.gz",
        "cross_mutau_mu_leg"            : f"{corr_dir}hleprare/TriggerScaleFactors/{year}{tag}/CrossMuTauHlt_MuLeg_v1.json",
        "HLT_mu_eff"                    : f"{corr_dir}hleprare/TriggerScaleFactors/{year}{tag}/MuHlt_abseta_pt_wEff.json",
        "electron_scaling_smearing"     : f"{jsonpog_dir}EGM/{year}_{pog_tag}/electronSS.json.gz",
        "electron_idiso"                : f"{jsonpog_dir}EGM/{year}_{pog_tag}/electron.json.gz",
        "electron_trigger"              : f"{jsonpog_dir}EGM/{year}_{pog_tag}/electronHlt.json.gz",
        "tau_correction"                : f"{jsonpog_tau_dir}TAU/{year}_{tag}/tau_DeepTau2018v2p5_{year}_{tag}.json.gz",
        "zpt_weight"                    : f"{corr_dir}/dy_ptll/DY_pTll_recoil_corrections_{year}{long_tag}.json.gz",
        "jet_jerc"                      : (f"{jsonpog_dir}JME/{year}_{pog_tag}/jet_jerc.json.gz", "v2"),
        "jet_veto_map"                  : (f"{jsonpog_dir}JME/{year}_{pog_tag}/jetvetomaps.json.gz", "v2"),
        "fake_factors"                  : (f"{corr_dir}fake_factors_{channel}_22and23_mt{cfg.x.mt_cut_value}_4bins.json", "v2"),
        "met_recoil"                    : (f"{corr_dir}hleprare/RecoilCorrlib/Recoil_corrections_{cfg.x.year}{tag}_v2.json.gz", "v2"),
        #"met_phi_corr": (f"{jsonpog_dir}JME/{cfg.x.year}{pog_tag}/met{cfg.x.year}.json.gz", "v2"), #FIXME: there is no json present in the jsonpog-integration for this year, I retrieve the json frm: https://cms-talk.web.cern.ch/t/2022-met-xy-corrections/53414/2 but it seems corrupted
        "ip_corr"                       : f"{corr_dir}ip_correction/ip_correction_Run3_{year}{short_tag}.json",
        "ml_model_even"                 : f"{ml_dir}{ch_short}/EVEN/model_{ch_short}_EVEN.json",
        "ml_model_odd"                  : f"{ml_dir}{ch_short}/ODD/model_{ch_short}_ODD.json",
        "filter_eff"                    : f"{corr_dir}/filter_eff/2024_v2/Run3_{year}{short_tag}.yaml",
    })
    
    
    cat_path = "/cvmfs/cms-griddata.cern.ch/cat/metadata"
    cat_tag = tags['cat_tag']
    
    jsons_2025_v2 = DotDict.wrap({
        "lumi": {
            "golden": (golden_ls[year], "v1"),
            "normtag": ("/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json", "v1"), #/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags
        },
        "pu_sf"                         : (f"{jsonpog_dir}LUM/{year}_{pog_tag}/puWeights.json.gz", "v1"),
        "muon_correction"               : f"{jsonpog_dir}MUO/{year}_{pog_tag}/muon_Z.json.gz",
        "cross_mutau_mu_leg"            : f"{corr_dir}hleprare/TriggerScaleFactors/{year}{tag}/CrossMuTauHlt_MuLeg_v1.json",
        "electron_scaling_smearing"     : f"{jsonpog_dir}EGM/{year}_{pog_tag}/electronSS.json.gz",
        "electron_idiso"                : f"{jsonpog_dir}EGM/{year}_{pog_tag}/electron.json.gz",
        "electron_trigger"              : f"{jsonpog_dir}EGM/{year}_{pog_tag}/electronHlt.json.gz",
        "tau_correction"                : f"{cat_path}/TAU/{cat_tag}/2025-10-01/tau_DeepTau2018v2p5_{year}_{tag}.json.gz",
        "tes"                           : (f"{corr_dir}measured_by_ic/tes/tau_es_dm_DeepTau2018v2p5_{year}_{tag}.json.gz", "v1"),
        "tau_sf"                        : f"{corr_dir}measured_by_ic/tau_sf/tau_sf_pt-dm_DeepTau2018v2p5VSjet_{year}_{tag}.json.gz",
        "tau_trigger_sf"                : f"{corr_dir}measured_by_ic/tau_trigger_sf/tau_trigger_DeepTau2018v2p5_{year}_{tag}.json.gz",
        "zpt_weight"                    : f"{corr_dir}dy_ptll/DY_pTll_weights_{year}{tag}.json.gz",
        "met_recoil"                    : f"{corr_dir}dy_ptll/DY_pTll_recoil_corrections_{year}{tag}.json.gz",
        "jet_jerc"                      : (f"{jsonpog_dir}JME/{year}_{pog_tag}/jet_jerc.json.gz", "v2"),
        "jet_veto_map"                  : (f"{jsonpog_dir}JME/{year}_{pog_tag}/jetvetomaps.json.gz", "v2"),
        "fake_factors"                  : (f"{corr_dir}fake_factors_v1_2025_{channel}_22and23_mt{cfg.x.mt_cut_value}_dr_iso0p05_mt70.json", "v2"),
        "ip_corr"                       : f"{corr_dir}ip_correction/ip_correction_Run3_{year}{short_tag}.json",
        "ml_model_even"                 : f"{corr_dir}signal_classifier/model_EVEN.json",
        "ml_model_odd"                  : f"{corr_dir}signal_classifier/model_ODD.json",
        "filter_eff"                    : f"{corr_dir}filter_eff/2025_v1/Run3_{year}{short_tag}.yaml",
        "stitching"                     : f"{corr_dir}stitching_weights.json",

    })
    
    
    
    cfg.x.external_files = jsons_2025_v2
    #--------------------------------------------------------------------------------------------- #
    # electron settings
    # names of electron correction sets and working points
    # (used in the electron_sf producer)
    # --------------------------------------------------------------------------------------------- #
    cfg.x.electron_sf = DotDict.wrap({
        'ID': {'corrector': "Electron-ID-SF",
               'year': tags['e_sf_tag'],
               'wp':"wp80iso"},
        'scale': {'corrector': tags['e_scale_corrector']},
        'smearing': {'corrector': tags['e_smearing_corrector']},
        'trig': {'corrector': "Electron-HLT-SF",
                 'year': tags['e_sf_tag'],
                 'wp': "HLT_SF_Ele30_TightID"},
        'xtrig': {'corrector': "Electron-HLT-SF",
                  'year': tags['e_sf_tag'],
                  'wp': "HLT_SF_Ele24_TightID"}
    })
    
    # --------------------------------------------------------------------------------------------- #
    # muon settings
    # names of muon correction sets and working points
    # (used in the muon producer)
    # --------------------------------------------------------------------------------------------- #
    cfg.x.muon_id_wp = 'Medium'
    cfg.x.muon_sf = DotDict.wrap({ 
                                  
        'ID': {'corrector': f"NUM_{cfg.x.muon_id_wp}ID_DEN_TrackerMuons",
               },
        
        'iso': {'corrector': f"NUM_TightPFIso_DEN_{cfg.x.muon_id_wp}ID",
                },
        
        'trig': {'corrector': f"NUM_IsoMu24_DEN_CutBasedId{cfg.x.muon_id_wp}_and_PFIsoMedium", #TODO Worging points for trig and iso does not match!
                 },
        
        'xtrig': {'corrector': f"NUM_IsoMu20_DEN_CutBasedId{cfg.x.muon_id_wp}_and_PFIsoMedium", 
                  },
    })
    
    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    
    ##########
    # shifts #
    ##########
    from httcp.config.variables import keep_columns
    keep_columns(cfg)
 
    cfg.add_shift(name="nominal", id=0)

   
    for i, (match, dm) in enumerate(itertools.product(["jet", "e"], [0, 1, 2, 10, 11])):
        cfg.add_shift(name=f"tec_{match}_dm{dm}_up", id=20 + 2 * i, type="shape", tags={"tec"})
        cfg.add_shift(name=f"tec_{match}_dm{dm}_down", id=21 + 2 * i, type="shape", tags={"tec"})
        add_shift_aliases(
            cfg,
            f"tec_{match}_dm{dm}",
            {
                "Tau.pt": "Tau.pt_{name}",
                "Tau.mass": "Tau.mass_{name}",
                f"{cfg.x.met_name}.pt": f"{cfg.x.met_name}.pt_{{name}}",
                f"{cfg.x.met_name}.phi": f"{cfg.x.met_name}.phi_{{name}}",
            },
        )
        
    cfg.add_shift(name="tau_up", id=1, type="shape")
    cfg.add_shift(name="tau_down", id=2, type="shape")
    add_shift_aliases(cfg, "tau", {"tau_weight": "tau_weight_{direction}"})
    
    cfg.add_shift(name="mu_up", id=3, type="shape")
    cfg.add_shift(name="mu_down", id=4, type="shape")
    add_shift_aliases(cfg, "mu", {"muon_weight": "muon_weight_{direction}"})
    
    # cfg.add_shift(name="ts_up", id=5, type="shape") #cp-even
    # cfg.add_shift(name="ts_down", id=7, type="shape") #cp-odd
    # add_shift_aliases(cfg, "ts", {"tauspinner_weight": "tauspinner_weight_{direction}"})
    
    # cfg.add_shift(name="electron_up", id=8, type="shape")
    # cfg.add_shift(name="electron_down", id=9, type="shape")
    # add_shift_aliases(cfg, "electron", {"electron_weight": "electron_weight_{direction}"})
    
    cfg.add_shift(name="top_pt_up", id=10, type="shape")
    cfg.add_shift(name="top_pt_down", id=11, type="shape")
    add_shift_aliases(cfg, "top_pt", {"top_pt_weight": "top_pt_weight_{direction}"})
    
    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    get_shifts = functools.partial(get_shifts_from_sources, cfg)   
 
    cfg.x.event_weights = DotDict({
        "normalization_weight": [],
        "filter_weight": [],
        #"mc_weight":[],
        "tau_weight_nom": get_shifts("tau"),
        "pu_weight": [],
        "tauspinner_weight": [],
        "zpt_weight":[],
        "muon_weight_nom": get_shifts("mu"),
        #"electron_weight_nom": get_shifts("electron"), 
        "top_pt_weight" : [],       
        "trigger_weight_mutau_nom": [],
        "stitching_weight": [],
    })
    for dataset in cfg.datasets:
        if dataset.has_tag("ttbar"):
            dataset.x.event_weights = {"top_pt_weight": get_shifts("top_pt")} 
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
        return version if version else 'dev'
            
        
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
    channel_id = {'etau': 1,
                  'mutau': 2,
                  'tautau': 4}
    cfg.add_channel(name=channel,   id=channel_id[channel])
    
    cfg.x.ch_objects = DotDict.wrap({
        'etau'   : {'lep0' : 'Electron',
                    'lep1' : 'Tau'},
        'mutau'  : {'lep0' : 'Muon',
                    'lep1' : 'Tau'},
        'tautau' : {'lep0' : 'Tau',
                    'lep1' : 'Tau'},
    })
    cfg.x.fake_factor_method = DotDict.wrap({
    "axes": {'tau_pt': {
                'var_route' : [f'hcand_{channel}', 'lep1', 'pt'],
                'ax_str'    : 'Variable([20,30,40,60,200], name="tau_pt", label="Tau pt", underflow=False, overflow=False)',
                },
             'tau_dm_pnet': {
                'var_route' : [f'hcand_{channel}', 'lep1', 'decayModePNet'],
                'ax_str'    : 'IntCategory([0,1,2,10,11], name="tau_dm_pnet", label="Tau PNet decayMode")',
             },
             "n_jets": {
                'var_route' : ['n_jets'],
                'ax_str'    : 'Integer(0, 3, name="n_jets", label="Number of jets",underflow=False, overflow=False)',
            },
    },
    "columns" : ['ff_weight_wj','ff_weight_qcd'],
    "shifts"  : ["up", "nom", "down"]
    })
    
    cfg.x.dy_ptll_corrs = DotDict.wrap({
        'datasets' : {
            "DYto2L_M_10to50_amcatnloFXFX": "NLO",
            "DYto2L_M_50_0J_amcatnloFXFX": "NLO",
            "DYto2L_M_50_1J_amcatnloFXFX": "NLO",
            "DYto2L_M_50_2J_amcatnloFXFX": "NLO",
            "DYto2L_M_50_amcatnloFXFX": "NLO",
            # DY->tautau
            "DYto2Tau_MLL_50_0J_amcatnloFXFX": "NLO",
            "DYto2Tau_MLL_50_1J_amcatnloFXFX": "NLO",
            "DYto2Tau_MLL_50_2J_amcatnloFXFX": "NLO",
            
            "WtoLNu_amcatnloFXFX": "NLO",
            },
       
       
        # 'datasets' : {'dy_lep_madgraph'  : "LO",
        #               'wj_incl_madgraph' : "LO"},
    })
    
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
        # define custom remote fs's to look at
        cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: "wlcg_fs_eos"
        
    # add categories using the "add_category" tool which adds auto-generated ids
    from httcp.config.categories import add_categories
    add_categories(cfg,channel=channel)
        
    from httcp.config.variables import add_variables
    add_variables(cfg)
    
    from data_driven.hist_hooks import add_hist_hooks
    add_hist_hooks(ana)
    
    from httcp.config.styles import setup_plot_styles
    setup_plot_styles(cfg)
