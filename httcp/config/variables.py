### This config is used for listing the variables used in the analysis ###

from columnflow.config_util import add_category

import order as od

from columnflow.columnar_util import EMPTY_FLOAT
from columnflow.util import DotDict
from columnflow.columnar_util import ColumnCollection

def keep_columns(cfg: od.Config) -> None:
    # columns to keep after certain steps
    cfg.x.keep_columns = DotDict.wrap({
        "cf.ReduceEvents": {
            # TauProds
            "TauProd.*",
            # general event info
            "run", "luminosityBlock", "event",
            "PV.npvs","Pileup.nTrueInt","Pileup.nPU","genWeight", "LHEWeight.originalXWGTUP",
        } | {
            f"PuppiMET.{var}" for var in [
                "pt", "phi", "significance",
                "covXX", "covXY", "covYY",
            ]
        } | {
            f"MET.{var}" for var in [
                "pt", "phi", "significance",
                "covXX", "covXY", "covYY",
            ]     
        } | {
            f"Jet.{var}" for var in [
                "pt", "eta", "phi", "mass", 
                "btagDeepFlavB", "hadronFlavour"
            ] 
        } | {
            f"Tau.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", 
                "rawDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSjet", "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu", 
                "decayMode", "decayModePNet", "genPartFlav", "rawIdx",
                "pt_no_tes", "mass_no_tes"
            ] 
        } | {
            f"Muon.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge",
		"decayMode", "pfRelIso04_all","mT", "rawIdx"
            ] 
        } | {
            f"Electron.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", 
                "decayMode", "pfRelIso03_all", "mT", "rawIdx"
            ] 
        } | {
            f"{var}_triggerd" for var in [ #Trigger variables to have a track of a particular trigger fired
                "single_electron", "cross_electron",
                "single_muon", "cross_muon",
                "cross_tau",
            ]
        } | {
            f"matched_triggerID_{var}" for var in [
                "e", "mu", "tau",
            ]
        } | {
            f"TrigObj.{var}" for var in [
                "id", "pt", "eta", "phi", "filterBits",
            ]
        } | {
            f"hcand.{var}" for var in [
                "pt","eta","phi","mass", "charge", 
                "decayMode", "rawIdx"
            ]
        } | {
            "GenTau.*", "GenTauProd.*",
        } | {
            f"hcandprod.{var}" for var in [
                "pt", "eta", "phi", "mass", "charge",
                "pdgId", "tauIdx",
            ]
        } | {ColumnCollection.ALL_FROM_SELECTOR},
        "cf.MergeSelectionMasks": {
            "normalization_weight", 
            "cutflow.*", "process_id", "category_ids",
        },
        "cf.UniteColumns": {
            "*",
        },
    })



def add_common_features(cfg: od.config) -> None:
    """
    Adds common features
    """
    cfg.add_variable(
        name="event",
        expression="event",
        binning=(1, 0.0, 1.0e9),
        x_title="Event number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="run",
        expression="run",
        binning=(1, 100000.0, 500000.0),
        x_title="Run number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="lumi",
        expression="luminosityBlock",
        binning=(1, 0.0, 5000.0),
        x_title="Luminosity block",
        discrete_x=True,
    )


def add_lepton_features(cfg: od.Config) -> None:
    """
    Adds lepton features only , ex electron_1_pt
    """
    for obj in ["Electron", "Muon", "Tau"]:
        for i in range(2):
            cfg.add_variable(
                name=f"{obj.lower()}_{i+1}_pt",
                expression=f"{obj}.pt[:,{i}]",
                null_value=EMPTY_FLOAT,
                binning=(40, 0., 200.),
                unit="GeV",
                x_title=obj + r" $p_{T}$",
            )
            cfg.add_variable(
                name=f"{obj.lower()}_{i+1}_phi",
                expression=f"{obj}.phi[:,{i}]",
                null_value=EMPTY_FLOAT,
                binning=(32, -3.2, 3.2),
                x_title=obj + r" $\phi$",
            )
            cfg.add_variable(
                name=f"{obj.lower()}_{i+1}_eta",
                expression=f"{obj}.eta[:,{i}]",
                null_value=EMPTY_FLOAT,
                binning=(25, -2.5, 2.5),
                x_title=obj + r" $\eta$",
            )
        cfg.add_variable(
            name=f"{obj.lower()}_mT",
            expression=f"{obj}.mT",
            null_value=EMPTY_FLOAT,
            binning=(40, 0.0, 200.0),
            unit="GeV",
            x_title=obj + r"$m_{T}$",
    )


def add_jet_features(cfg: od.Config) -> None:
    """
    Adds jet features only
    """
    cfg.add_variable(
        name="n_jet",
        expression="n_jet",
        binning=(11, -0.5, 10.5),
        x_title="Number of jets",
        discrete_x=True,
    )
    cfg.add_variable(
        name="jets_pt",
        expression="Jet.pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$p_{T} of all jets$",
    )
    for i in range(2):
        cfg.add_variable(
            name=f"jet_{i+1}_pt",
            expression=f"Jet.pt[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, 0.0, 400.0),
            unit="GeV",
            x_title=r"Jet $p_{T}$",
        )
        cfg.add_variable(
            name=f"jet_{i+1}_eta",
            expression=f"Jet.eta[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(30, -3.0, 3.0),
            x_title=r"Jet $\eta$",
        )
    cfg.add_variable(
        name="ht",
        # expression=lambda events: ak.sum(events.Jet.pt, axis=1),
        expression="ht",
        binning=(40, 0.0, 800.0),
        unit="GeV",
        x_title="HT",
    )
    cfg.add_variable(
        name="jet_raw_DeepJetFlavB",
        expression="Jet.btagDeepFlavB",
        null_value=EMPTY_FLOAT,
        binning=(30, 0,1),
        x_title=r"raw DeepJetFlawB",
    )


def add_highlevel_features(cfg: od.Config) -> None:    
    """
    Adds MET and other high-level features
    """
    cfg.add_variable(
        name="met",
        expression="MET.pt",
        null_value=EMPTY_FLOAT,
        binning=(20, 0.0, 100.0),
        x_title=r"MET",
    )
    cfg.add_variable(
        name="puppi_met_pt",
        expression="PuppiMET.pt",
        null_value=EMPTY_FLOAT,
        binning=(50, 0,100),
        unit="GeV",
        x_title=r"MET $p_T$",
    )
    cfg.add_variable(
        name="puppi_met_phi",
        expression="PuppiMET.phi",
        null_value=EMPTY_FLOAT,
        binning=(30, -3,3),
        x_title=r"MET $\phi$",
    )
    cfg.add_variable(
        name="hcand_mass",
        expression="hcand_obj.mass",
        null_value=EMPTY_FLOAT,
        binning=(40, 0.0, 200.0),
        unit="GeV",
        x_title=r"$m_{vis}$",
    )
    
    
    
    

def add_weight_features(cfg: od.Config) -> None:
    """
    Adds weights
    """
    cfg.add_variable(
        name="mc_weight",
        expression="mc_weight",
        binning=(200, -10, 10),
        x_title="MC weight",
    )
    cfg.add_variable(
        name="pu_weight",
        expression="pu_weight",
        null_value=EMPTY_FLOAT,
        binning=(30, 0,3),
        unit="GeV",
        x_title=r"Pileup weight",
    )

    

def add_cutflow_features(cfg: od.Config) -> None:
    """
    Adds cf features
    """
    cfg.add_variable(
        name="cf_jet1_pt",
        expression="cutflow.jet1_pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"Jet 1 $p_{T}$",
    )


def add_hcand_features(cfg: od.Config) -> None:
    """
    Adds h lepton features only
    """
    """
    for i in range(2):
        cfg.add_variable(
            name=f"hlepton_{i+1}_pt",
            expression=f"hcand_pt[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, 0., 200.),
            unit="GeV",
            x_title=f"lepton_{i+1}" + r" $p_{T}$",
        )
        cfg.add_variable(
            name=f"hlepton_{i+1}_eta",
            expression=f"hcand_eta[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(25, -2.5, 2.5),
            unit="GeV",
            x_title=f"lepton_{i+1}" + r" $\eta$",
    )
    """
    cfg.add_variable(
        name="hcand_invm",
        expression="hcand_invm",
        null_value=EMPTY_FLOAT,
        binning=(50, 0, 400),
        unit="GeV",
        x_title=r"$m_{ll}$",
    )
    cfg.add_variable(
        name="hcand_dr",
        expression="hcand_dr",
        null_value=EMPTY_FLOAT,
        binning=(40, 0, 5),
        x_title=r"$\Delta R(l,l)$",
    )
    cfg.add_variable(
        name="PhiCP_PVPV",
        expression="PhiCP_PVPV",
        null_value=EMPTY_FLOAT,
        binning=(16, 0, 6.4),
        x_title=r"$\Phi_{CP}^{PV-PV}$ (rad)",
    )
    cfg.add_variable(
        name="PhiCPGen_PVPV",
        expression="PhiCPGen_PVPV",
        null_value=EMPTY_FLOAT,
        binning=(16, 0, 6.4),
        x_title=r"$\Phi_{CP}^{PV-PV}$ (rad)",
    )
    cfg.add_variable(
        name="PhiCP_DPDP",
        expression="PhiCP_DPDP",
        null_value=EMPTY_FLOAT,
        binning=(16, 0, 6.4),
        x_title=r"$\Phi_{CP}^{DP-DP}$ (rad)",
    )
    cfg.add_variable(
        name="PhiCPGen_DPDP",
        expression="PhiCPGen_DPDP",
        null_value=EMPTY_FLOAT,
        binning=(16, 0, 6.4),
        x_title=r"$\Phi_{CP}^{DP-DP}$ (rad)",
    )

def add_test_variables(cfg: od.Config) -> None:
        cfg.add_variable(
            name="tau_pt_no_tes",
            expression="Tau.pt_no_tes",
            null_value=EMPTY_FLOAT,
            binning=(30, 25, 85),
            unit="GeV",
            x_title=r"tau $p_{T}$ (no TES)",
        )
        cfg.add_variable(
            name="mutau_mass_no_tes",
            expression="mutau_mass_no_tes",
            null_value=EMPTY_FLOAT,
            binning=(40, 0.0, 200.0),
            unit="GeV",
            x_title=r"$m_{vis}$(no TES)",
        )
    
         #single bin variables for transfer factor calculation
        cfg.add_variable(
            name="muon_eta_1bin",
            expression="Muon.eta",
            null_value=EMPTY_FLOAT,
            binning=(1, -3.0, 3.0),
            x_title=r"muon $\eta$",
        )
        cfg.add_variable(
            name="muon_pt_1bin",
            expression="Muon.pt",
            null_value=EMPTY_FLOAT,
            binning=(1, 20.0, 80.0),
            unit="GeV",
            x_title=r"muon $p_{T}$",
        )
        cfg.add_variable(
            name="muon_phi_1bin",
            expression="Muon.phi",
            null_value=EMPTY_FLOAT,
            binning=(1, -3.14159, 3.14159),
            x_title=r"muon $\varphi$",
        )
        cfg.add_variable(
            name="mutau_mass_1bin",
            expression="mutau_mass",
            null_value=EMPTY_FLOAT,
            binning=(1, 0.0, 200.0),
            unit="GeV",
            x_title=r"$m_{vis}$",
        )
        

def add_variables(cfg: od.Config) -> None:
    """
    Adds all variables to a *config*.
    """
    add_common_features(cfg)
    add_lepton_features(cfg)
    add_jet_features(cfg)
    add_highlevel_features(cfg)
    add_hcand_features(cfg)
    add_weight_features(cfg)
    add_cutflow_features(cfg)
    add_test_variables(cfg)
