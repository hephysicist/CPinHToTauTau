"""
Configuration of the HCPToTauTau analysis.
"""

import os
import functools

import law
import order as od
from scinum import Number
from typing import Optional

from columnflow.util import DotDict, maybe_import, dev_sandbox
from columnflow.columnar_util import EMPTY_FLOAT
from columnflow.config_util import (
    get_root_processes_from_campaign, add_shift_aliases, get_shifts_from_sources,
    verify_config_processes,
)

from httcp.util import get_dataset_lfns

ak = maybe_import("awkward")


def add_config(
        analysis: od.Analysis,
        campaign: od.Campaign,
        config_name: Optional[str] = None,
        config_id: Optional[int] = None,
        limit_dataset_files: Optional[int] = None,
) -> od.Config:
    print(f"ADDING CONFIG ---> {config_name} ---> ID: {config_id}")

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)

    # create a config by passing the campaign, so id and name will be identical
    cfg = analysis.add_config(campaign, name=config_name, id=config_id)
    #cfg = analysis.add_config(campaign)

    # gather campaign data
    year = campaign.x.year


    # add processes we are interested in
    process_names = [
        #"data",
        #"tt",
        #"dy",
        "h_ggf_tautau",
        #"st",
        #"ewk",
        #"vv",
    ]
    for process_name in process_names:
        # add the process
        proc = cfg.add_process(procs.get(process_name))
        # configuration of colors, labels, etc. can happen here
        if proc.is_mc:
            proc.color1 = (244, 182, 66) if proc.name == "tt" else (244, 93, 66)

    # add datasets we need to study
    dataset_names = [
        ### data ###
        #"data_mu_b",
        ### backgrounds ###
        #"dy_lep_m50_madgraph",
        "h_ggf_tautau_powheg",
        #"tt_sl_powheg",
        #"tt_dl_powheg",
        ### signals ###
    ]

    for dataset_name in dataset_names:
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))
        
        # for testing purposes, limit the number of files to 2
        
        for info in dataset.info.values():
            if limit_dataset_files: 
                info.n_files = min(info.n_files, limit_dataset_files)
            else: 
                info.n_files = min(info.n_files, 1)

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

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
    cfg.x.process_groups = {}
    
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

    # custom method and sandbox for determining dataset lfns
    cfg.x.get_dataset_lfns = None
    cfg.x.get_dataset_lfns_sandbox = None
    
    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False
    
    # lumi values in inverse pb
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2?rev=2#Combination_and_correlations
    if year == 2016:
        cfg.x.luminosity = Number(36310, {
            "lumi_13TeV_2016": 0.01j,
            "lumi_13TeV_correlated": 0.006j,
        })
    elif year == 2017:
        cfg.x.luminosity = Number(41480, {
            "lumi_13TeV_2017": 0.02j,
            "lumi_13TeV_1718": 0.006j,
            "lumi_13TeV_correlated": 0.009j,
        })
    else:  # 2018
        cfg.x.luminosity = Number(59830, {
            "lumi_13TeV_2017": 0.015j,
            "lumi_13TeV_1718": 0.002j,
            "lumi_13TeV_correlated": 0.02j,
        })
        

    # b-tag working points
    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17?rev=15
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
        },
    )
    # names of electron correction sets and working points
    # (used in the muon producer)
    cfg.x.electron_sf_names = ("UL-Electron-ID-SF", f"{year}", "wp80iso")

    # names of muon correction sets and working points
    # (used in the muon producer)
    cfg.x.muon_sf_names = ("NUM_TightRelIso_DEN_TightIDandIPCut", f"{year}_UL")
    
    # register shifts
    cfg.add_shift(name="nominal", id=0)

    cfg.x.deep_tau = DotDict.wrap({
        "tagger": "DeepTau2018v2p5",
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


    
    # tune shifts are covered by dedicated, varied datasets, so tag the shift as "disjoint_from_nominal"
    # (this is currently used to decide whether ML evaluations are done on the full shifted dataset)
    cfg.add_shift(name="tune_up", id=1, type="shape", tags={"disjoint_from_nominal"})
    cfg.add_shift(name="tune_down", id=2, type="shape", tags={"disjoint_from_nominal"})
    
    # fake jet energy correction shift, with aliases flaged as "selection_dependent", i.e. the aliases
    # affect columns that might change the output of the event selection
    cfg.add_shift(name="jec_up", id=20, type="shape")
    cfg.add_shift(name="jec_down", id=21, type="shape")
    add_shift_aliases(
        cfg,
        "jec",
        {
            "Jet.pt": "Jet.pt_{name}",
            "Jet.mass": "Jet.mass_{name}",
            "MET.pt": "MET.pt_{name}",
            "MET.phi": "MET.phi_{name}",
        },
    )
    
    # event weights due to muon scale factors
    cfg.add_shift(name="mu_up", id=10, type="shape")
    cfg.add_shift(name="mu_down", id=11, type="shape")
    add_shift_aliases(cfg, "mu", {"muon_weight": "muon_weight_{direction}"})
    
    # external files
    import os
    external_path = os.path.join(os.environ.get('HTTCP_BASE'), "httcp/data/corrections")
    cfg.x.external_files = DotDict.wrap({
        # lumi files
        "lumi": {
            "golden": (f"{external_path}/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt", "v1"),  # noqa
            "normtag": (f"{external_path}/normtag_PHYSICS.json", "v1"),
        },
        "pileup":{
            "data" : f"{external_path}/Data_PileUp_2022_preEE.root", #TODO: make a link to the common correction repo
            "mc"   : f"{external_path}/MC_PileUp_2022.root" #TODO: make a link to the common correction repo
        },
        "muon_correction" : f"{external_path}/muon_SFs_2022_preEE.root", #TODO: make a link to the common correction repo
        "tau_correction"  : f"{external_path}/tau_DeepTau2018v2p5_2022_preEE.json.gz", #TODO: make a link to the common correction repo
    })

    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    
    # columns to keep after certain steps
    from httcp.config.variables import keep_columns
    keep_columns(cfg)
    
    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight"  : [],
        "pu_weight"             : [],
        "muon_weight"           : [],
        "tau_id_sf"             : [],
    })
    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    cfg.x.versions = {
        # "cf.CalibrateEvents": "prod1",
        # "cf.SelectEvents": (lambda cls, inst, params: "prod1" if params.get("selector") == "default" else "dev1"),
        # ...
    }
    
    # channels
    cfg.add_channel(name="etau",   id=1)
    cfg.add_channel(name="mutau",  id=2)
    cfg.add_channel(name="tautau", id=4)


    # add categories using the "add_category" tool which adds auto-generated ids
    # the "selection" entries refer to names of selectors, e.g. in selection/example.py
    from httcp.config.categories import add_categories
    add_categories(cfg)
    
    # add variables
    # (the "event", "run" and "lumi" variables are required for some cutflow plotting task,
    # and also correspond to the minimal set of columns that coffea's nano scheme requires)
    # add variables
    from httcp.config.variables import add_variables
    add_variables(cfg)
    
    # add triggers
    if year == 2017:
        from httcp.config.triggers import add_triggers_2017
        add_triggers_2017(cfg)
    else:
        raise NotImplementedError(f"triggers not implemented for {year}")
    
    # add met filters
    from httcp.config.met_filters import add_met_filters

    add_met_filters(cfg)
    
    cfg.x.get_dataset_lfns = get_dataset_lfns
    
    # define a custom sandbox
    #cfg.x.get_dataset_lfns_sandbox = dev_sandbox("bash::$CF_BASE/sandboxes/cf.sh")
    
    # define custom remote fs's to look at
    #cfg.x.get_dataset_lfns_remote_fs = lambda dataset_inst: f"wlcg_fs_{cfg.campaign.x.custom['name']}"
    #cfg.x.get_dataset_lfns_remote_fs = lambda dataset_inst: f"{local_fs}"
