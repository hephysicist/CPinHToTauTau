# coding: utf-8

"""
Configuration of the higgs_cp analysis.
"""

import functools

import law
import order as od
from scinum import Number

from columnflow.util import DotDict, maybe_import, dev_sandbox
from columnflow.config_util import (
    get_root_processes_from_campaign, get_shifts_from_sources,
    add_shift_aliases, get_shifts_from_sources,
    verify_config_processes,
)

ak = maybe_import("awkward")


def add_run2_2018 (ana: od.Analysis,
                      campaign: od.Campaign,
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
    
    # add processes we are interested in
    process_names = [
        "data",
        "h_ggf_htt",
        #Drell-Yan
        #"dy_lep",
        #"dy_z2mumu",
        #"dy_z2tautau",
        #W + jets
        #"wj",
        #diboson
        # "vv", #diboson inclusive
        # "ww",
        # "wz",
        # "zz",
        # #ttbar 
        # "tt", # ttbar inclusive
        # "tt_sl",
        # "tt_dl",
        # "tt_fh",
        # #single top
        # "st",
        # "st_tchannel_t",
        # "st_twchannel",
    ]
    for process_name in process_names:
        # add the process
        proc = cfg.add_process(procs.get(process_name))
        if proc.is_mc:
            if proc.name == "dy_lep": proc.color1 = (223,102,72)
            if proc.name == "h_ggf_htt": proc.color1 = (51,53,204)
            if proc.name == "wj": proc.color1 = (201,89,84)
            if proc.name == "tt_sl": proc.color1 = (153,153,204)
            if proc.name == "tt_dl": proc.color1 = (184,184,227)
            if proc.name == "tt_fh": proc.color1 = (87,87,141)
            if proc.name == "ww" : proc.color1 = (102,204,102)
            if proc.name == "wz" : proc.color1 = (49,157,49)
            if proc.name == "zz" : proc.color1 = (120,214,120)
            
            if proc.name == "vv" : proc.color1 = (102,204,102)

        # configuration of colors, labels, etc. can happen here
        

    # add datasets we need to study
    dataset_names = [
        "data_mu_a",
        "h_ggf_htt"
        #"data_mu_c",
        #"data_mu_d",
        #"data_mu_e",
        #"data_mu_f",
        #Drell-Yan
        #"dy_incl",
        #W+jets
        #"wj_incl",
        # #Diboson
        # "ww",
        # "wz",
        # "zz",
        # #ttbar
        # "tt_sl",
        # "tt_dl",
        # "tt_fh",
        # #single top t-channel
        # "st_t_bbarq",
        # "st_tbar_bq",
        # # single top tW channel
        # "st_t_wminus_to_lnu2q",
        # "st_t_wminus_to_2l2nu",
        # "st_tbar_wplus_to_lnu2q",
        # "st_tbar_wplus_to_2l2nu",
        # "signal"
        ]
    
    for dataset_name in dataset_names:
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))

        # for testing purposes, limit the number of files to 1
        for info in dataset.info.values():
            if limit_dataset_files:
                info.n_files = min(info.n_files, limit_dataset_files) #<<< REMOVE THIS FOR THE FULL DATASET

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

  
    #from httcp.config.triggers import add_triggers_run2_UL2017
    #add_triggers_run2_UL2017(cfg)
    from httcp.config.triggers import add_triggers_2018
    add_triggers_2018(cfg)
    
    from httcp.config.met_filters import add_met_filters
    add_met_filters(cfg)

    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    cfg.x.default_calibrator = "main"
    cfg.x.default_selector = "main"
    cfg.x.default_producer = "main"
    cfg.x.default_ml_model = None
    cfg.x.default_inference_model = "example"
    cfg.x.default_categories = ("incl",)
    #cfg.x.default_variables = ("n_jet", "jet1_pt")
    cfg.x.default_variables = ("event","channel_id")

    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {
        "diboson": ["ww", "wz", "zz"],
        "tt" : ["tt_sl","tt_dl","tt_fh"]
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
    #  cfg.x.selector_step_labels = {"Initial":0, 
    #                                "Trigger": , "Muon"}
     
    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False

    # lumi values in inverse pb
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2?rev=2#Combination_and_correlations
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#DATA_AN2
    #Only F and G eras
    # lumi values in inverse pb
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2?rev=2#Combination_and_correlations
    # difference pre-post VFP: https://cds.cern.ch/record/2854610/files/DP2023_006.pdf
    if  cfg.x.year == 2016:
        if campaign.x.vfp == "pre":
            cfg.x.luminosity = Number(19_500, {
                "lumi_13TeV_2016": 0.01j,
                "lumi_13TeV_correlated": 0.006j,
            })
        if campaign.x.vfp == "post":
            cfg.x.luminosity = Number(16_800, {
                "lumi_13TeV_2016": 0.01j,
                "lumi_13TeV_correlated": 0.006j,
            })
    elif  cfg.x.year == 2017:
        cfg.x.luminosity = Number(41_480, {
            "lumi_13TeV_2017": 0.02j,
            "lumi_13TeV_1718": 0.006j,
            "lumi_13TeV_correlated": 0.009j,
        })
    else:  # 2018
        cfg.x.luminosity = Number(59_830, {
            "lumi_13TeV_2017": 0.015j,
            "lumi_13TeV_1718": 0.002j,
            "lumi_13TeV_correlated": 0.02j,
        })
    # minimum bias cross section in mb (milli) for creating PU weights, values from
    # https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData?rev=45#Recommended_cross_section
    cfg.x.minbias_xs = Number(69.2, 0.046j)


    

    # names of muon correction sets and working points
    # (used in the muon producer)
    #cfg.x.muon_sf_names = ("NUM_TightRelIso_DEN_TightIDandIPCut", f"{cfg.x.year}")

   
  
    cfg.x.deep_tau = DotDict.wrap({
        "tagger": "DeepTau2017v2p1", # Placeholder!!!! WPs are for "DeepTau2018v2p5",
        "vs_e"          : "VVLoose",
        "vs_mu"         : "Tight",
        "vs_jet"        : "Medium",
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
            2017 : {
                "deepcsv": { #https://btv-wiki.docs.cern.ch/ScaleFactors/UL2017/
                    "loose": 0.1355,
                    "medium": 0.4506,
                    "tight": 0.7738,
                },
                "deepjet": {
                    "loose": 0.0532,
                    "medium": 0.3040,
                    "tight": 0.7476,
                },
            },
            2018 : {
                "deepcsv": { #https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
                    "loose": 0.1208,
                    "medium": 0.4168 ,
                    "tight": 0.7665,
                },
                "deepjet": {
                    "loose": 0.0490,
                    "medium": 0.2783,
                    "tight": 0.7100,
                },
            },                             
        },
    )
    corr_dir = "/afs/desy.de/user/s/stzakhar/nfs/CPinHToTauTau/httcp/corrections"
    
    def tag_caster(campaign: od.Campaign) -> str:
        #Helper function to cast campaign tags to the tags used in POG groups for the scale factors
        year = campaign.x.year
        tag = campaign.x.tag
        out_tag = ''
        if year in [2017,2018]  : out_tag = '_UL'
        elif tag == "preEE"     : out_tag = "_Summer22"
        elif tag == "postEE"    : out_tag = "_Summer22EE"
        elif tag == "preBpix"   : out_tag = "_Summer23"
        elif tag == "postBpix"  : out_tag = "_Summer23BPix"
        elif tag == "preVFP"    : out_tag = "preVFP_UL"
        elif tag == "postVFP"   : out_tag = "postVFP_UL"
        return out_tag
    
    tag = tag_caster(campaign)
    cfg.x.external_files = DotDict.wrap({
        # lumi files
        "lumi": {
            "golden": (f"{corr_dir}/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt", "v1"),  # https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
            "normtag": (f"{corr_dir}/normtag_PHYSICS.json", "v1"), #/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags
        },
        "pu_sf": (f"{corr_dir}/jsonpog-integration/POG/LUM/{cfg.x.year}{tag}/puWeights.json.gz", "v1"),
        "muon_correction" : f"{corr_dir}/jsonpog-integration/POG/MUO/{cfg.x.year}{tag}/muon_Z.json.gz",
        "tau_correction"  : f"{corr_dir}/jsonpog-integration/POG/TAU/{cfg.x.year}{tag}/tau.json.gz",
    })

    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    
    from httcp.config.variables import keep_columns
    keep_columns(cfg)
    
    # register shifts
    cfg.add_shift(name="nominal", id=0)

    cfg.add_shift(name="tau_up", id=1, type="shape")
    cfg.add_shift(name="tau_down", id=2, type="shape")
    add_shift_aliases(cfg, "tau", {"tau_weight": "tau_weight_{direction}"})
    
    cfg.add_shift(name="mu_up", id=3, type="shape")
    cfg.add_shift(name="mu_down", id=4, type="shape")
    add_shift_aliases(cfg, "mu", {"muon_weight": "muon_weight_{direction}"})
    
    
    
    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight"  : [],
        "pu_weight"             : [],
        "muon_weight_nom"       : [],
        "tau_weight_nom"        : []
    })
    
    cfg.x.default_weight_producer = "all_weights"

    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    def set_version(cls, inst, params):
        # per default, use the version set on the command line
        version = inst.version 
        return version if version else 'dev1'
            
        
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
    cfg.add_channel(name="etau",   id=1)
    cfg.add_channel(name="mutau",  id=2)
    cfg.add_channel(name="tautau", id=4)
    
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
    add_categories(cfg)
        
    from httcp.config.variables import add_variables
    add_variables(cfg)