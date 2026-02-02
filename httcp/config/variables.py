### This config is used for listing the variables used in the analysis ###

from columnflow.config_util import add_category

import order as od

from columnflow.columnar_util import EMPTY_FLOAT
from law.util import DotDict
from columnflow.columnar_util import ColumnCollection

from columnflow.util import maybe_import
np = maybe_import("numpy")

def keep_columns(cfg: od.Config) -> None:
    # columns to keep after certain steps
    cfg.x.keep_columns = DotDict.wrap({
        "cf.ReduceEvents": {
            # TauProds
            "TauProd.*",
            "GenPart.*",
            "GenVisTau.*",
            "GenZ.*",
            "GenVtx.*",
            # general event info
            "run", "luminosityBlock", "event",
            "PV.npvs","PV.x","PV.y","PV.z","PVBS.x","PVBS.y","PVBS.z",
            "Pileup.nTrueInt","Pileup.nPU","genWeight", "LHEWeight.originalXWGTUP", "HTXS_njets*","weight","zpt_weight","muon_weight_nom","mc_weight","tau_weight_nom",
            "LHE.Njets", "LHE.NpNLO",
        } | {
            f"PuppiMET.{var}" for var in [
                "pt", "phi", "significance",
                "covXX", "covXY", "covYY", "pt_no_jec", "eta_no_jec"
            ]
        } | {
            f"MET.{var}" for var in [
                "pt", "phi", "significance",
                "covXX", "covXY", "covYY",
            ]     
        } | {
            f"Jet.{var}" for var in [
                "pt", "eta", "phi", "mass", "jetId", 
                "btagDeepFlavB", "hadronFlavour", 
                "pt_no_jec", "phi_no_jec","eta_no_jec", 
                "mass_no_jec", "jec_no_jec_diff",
                "neEmEF","chHEF","neHEF",
                "chHEF","muEF","chEmEF",
                "neMultiplicity","chMultiplicity",
            ] 
        } | {
            f"Tau.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", 
                "rawDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSjet", "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu", 
                "decayMode", "decayModePNet", "genPartFlav", "rawIdx",
                "pt_no_tes", "mass_no_tes", "IPx", "IPy", "IPz", "ip_sig", "jetIdx", 
                "hasSV", "hasRefitSV", "IP_cov*", "refitSV*",
            ] 
        } | {
            f"Muon.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge",
                "decayMode", "pfRelIso04_all","mT", "rawIdx","IPx", "IPy", "IPz","ip_sig", "isPFcand",
                "genPartFlav","IP_cov*",
            ] 
        } | {
            f"Electron.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", 
                "decayMode", "pfRelIso03_all", "mT", "rawIdx", "IPx", "IPy", "IPz","ip_sig", "pt_no_scaling_smearing",
                "genPartFlav","IP_cov*",
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
            f"TauSpinner.weight_cp_{var}" for var in [
                "0", "0_alt", "0p25", "0p25_alt", "0p375",
                "0p375_alt", "0p5", "0p5_alt", "minus0p25", "minus0p25_alt"
            ]
        } | {
            "GenTau.*", "GenTauProd.*", "nJet", "N_b_jets", "n_jets", 
            "all_triggers_id", "triggerID_e", "triggerID_mu", "triggerID_mutau","triggerID_tau",'lead_jet.*',
            'sublead_jet.*','dijet.*'
        } | {
            f"hcandprod.{var}" for var in [
                "pt", "eta", "phi", "mass", "charge",
                "pdgId", "tauIdx"
            ]
        } | {
		"hcand_*","tau_decay_prods*", "OC_lepton_veto"
	} | {"is_b_vetoed","channel_id"} | {ColumnCollection.ALL_FROM_SELECTOR},
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
        binning=(1, 0, 2),
        x_title="Event number",
    )
    cfg.add_variable(
        name="N_events",
        expression="event",
        binning=(1, 0, 2),
        x_title="Number of events",
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
    cfg.add_variable( #PV.npvs = total number of reconstructed primary vertices https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_EGamma1_Run2024D-PromptReco-v1.html#SV
        name="npvs",
        expression="PV.npvs",
        binning=(20, 0, 100.0),
        x_title="Reconstructed primary vertices",
        null_value=EMPTY_FLOAT,
    )
    # cfg.add_variable(
    #     name="nPVBS",
    #     expression="nPVBS",
    #     binning=(20, 0, 100.0),
    #     x_title="Reconstructed primary vertices (beam-spot)",
    #     null_value=EMPTY_FLOAT,
    # )
    cfg.add_variable( #Pileup.nTrueInt
        name="pu_nTrue_Int",
        expression="Pileup.nTrueInt",
        binning=(20, 0.0, 100.0),
        x_title="True number of Pile Up per Event",
        null_value=EMPTY_FLOAT,
    )
    cfg.add_variable( #Pileup.nPU
        name="nPU",
        expression="Pileup.nPU",
        binning=(20, 0.0, 100.0),
        x_title="Number of Pile Up per Event",
        null_value=EMPTY_FLOAT,
    )

def add_lepton_features(cfg: od.Config) -> None:
    """
    Adds lepton features only , ex electron_1_pt
    """
    cfg.add_variable(
        name=f"electron_1_pt_no_scaling_smearing",
        expression="Electron.pt_no_scaling_smearing[:,0]",
        null_value=EMPTY_FLOAT,
        binning=(40, 0., 200.),
        unit="GeV",
        x_title= r" Electron $p_{T}$ no scaling or smearing",
    )
    
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
            name=f"{obj.lower()}_ip_sig",
            expression=f"{obj}.ip_sig",
            null_value=EMPTY_FLOAT,
            binning=(40, 0.0, 10),
            unit="",
            x_title=obj + r"$\frac{|IP|}{\sigma(IP)}$",
        )


def add_jet_features(cfg: od.Config) -> None:
    """
    Adds jet features only
    """
    cfg.add_variable(
        name="n_jet",
        expression="nJet",
        binning=(11, -0.5, 10.5),
        x_title="Number of jets",
        discrete_x=True,
    )
    cfg.add_variable(
        name="jets_pt_no_jec",
        expression="Jet.pt_no_jec",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$Uncorrected p_{T} of all jets$",
    )      
    cfg.add_variable(
        name="jets_pt",
        expression="Jet.pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"$p_{T} of all jets$",
    )
    cfg.add_variable(
        name="jets_eta",
        expression="Jet.eta",
        binning=(82, -5.191, 5.191),
        unit="",
        x_title="$\\eta$ of all jets",
    )
    cfg.add_variable(
        name="jets_phi",
        expression="Jet.phi",
        null_value=EMPTY_FLOAT,
        binning=(72, -np.pi, np.pi),
        x_title="$\\phi$ of all jets",
    )  
        
    cfg.add_variable(
        name="n_j",
        expression="n_jets",
        binning=(4, -0.5, 3.5),
        discrete_x=True,
        x_title="N_jets_pT_20_eta_4_7_Tight",
    )
    
    cfg.add_variable(
        name="N_jets_pT_20_eta_2_5_Tight",
        expression="n_jets_tag",
        binning=(11, -0.5, 10.5),
        discrete_x=True,
        x_title="N_jets_pT_20_eta_2_5_Tight",
    )
    cfg.add_variable(
        name="N_b_jets",
        expression="N_b_jets",
        binning=(11, -0.5, 10.5),
        discrete_x=True,
        x_title="N_b_jets",
    )
    cfg.add_variable(
        name="leading_jet_pt",
        expression="lead_jet.pt",
        null_value=EMPTY_FLOAT,
        binning=(30, 30.0, 330.0),
        unit="GeV",
        x_title=r"Leading jet $p_{T}$",
    )        
    cfg.add_variable(
        name="subleading_jet_pt",
        expression="sublead_jet.pt",
        null_value=EMPTY_FLOAT,
        binning=(25, 30.0, 280.0),
        unit="GeV",
        x_title=r"Subleading jet $p_{T}$",
    )
    cfg.add_variable(
        name="leading_jet_eta",
        expression="lead_jet.eta",
        null_value=EMPTY_FLOAT,
        binning=(47, -4.7, 4.7),
        x_title="Leading Jet $\\eta$",
    ) 
    cfg.add_variable(
        name="subleading_jet_eta",
        expression="sublead_jet.eta",
        null_value=EMPTY_FLOAT,
        binning=(47, -4.7, 4.7),
        x_title="Subleading Jet $\\eta$",
    ) 
    cfg.add_variable(
        name="leading_jet_phi",
        expression="lead_jet.phi",
        null_value=EMPTY_FLOAT,
        binning=(32, -3.2, 3.2),
        x_title="Leading Jet $\\phi$",
    )  
    cfg.add_variable(
        name="subleading_jet_phi",
        expression="sublead_jet.phi",
        null_value=EMPTY_FLOAT,
        binning=(32, -3.2, 3.2),
        x_title="Subleading Jet $\\phi$",
    ) 
    cfg.add_variable(
        name="dijet_delta_eta",
        expression="dijet.deltaeta",
        null_value=EMPTY_FLOAT,
        binning=(20,-6,6),
        x_title="$\\Delta \\eta_{jj}$",
    ) 
    cfg.add_variable(
        name="dijet_delta_phi",
        expression="dijet.deltaphi",
        null_value=EMPTY_FLOAT,
        binning=(20,-6,6),
        x_title="$\\Delta \\phi_{jj}$",
    ) 
    cfg.add_variable(
        name="dijet_pt",
        expression="dijet.pt",
        null_value=EMPTY_FLOAT,
        binning=(40, 0.0, 400.0),
        x_title="$pT_{jj}$",
    ) 
    cfg.add_variable(
        name="dijet_delta_r",
        expression="dijet.delta_r",
        null_value=EMPTY_FLOAT,
        binning=(20,0,5),
        x_title="$\\Delta R_{jj}$",
    ) 
    cfg.add_variable(
        name="mjj",
        expression="dijet.mass",
        null_value=EMPTY_FLOAT,
        binning=(40, 10.0, 410.0),
        unit="GeV",
        x_title=r"$m_{jj}$",
    )      
    cfg.add_variable(
        name="jet_jec_no_jec_diff",
        expression="Jet.jec_no_jec_diff",
        null_value=EMPTY_FLOAT,
        binning=(20,-10,10),
        x_title=r"$Jet_{jec} $p_{T} - $Jet_{no jec} $p_{T}",
    )
    cfg.add_variable(
        name="ht",
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

def add_hcp_bdt_output(cfg: od.Config) -> None:    
    for the_var in ['gtau','higgs','fake']:
        cfg.add_variable(
                name=f"bdt_raw_score_{the_var}",
                expression=f"bdt_raw_score_{the_var}",
                binning=(30, 0., 1.),
                x_title=f"raw BDT score for {the_var}",
            )
    cfg.add_variable(
                name=f"bdt_cat",
                expression=f"bdt_cat",
                binning=(4,-0.5,3.5),
                discrete_x=True,
                x_title=f"BDT class",
            )
    cfg.add_variable(
            name=f"bdt_raw_score_4bins_gtau",
            expression=f"bdt_raw_score_gtau",
            binning=[0,0.5,0.6,0.7,1],
            x_title=f"raw BDT score for gtau",
        )
    cfg.add_variable(
            name=f"bdt_raw_score_4bins_fake",
            expression=f"bdt_raw_score_fake",
            binning=[0,0.5,0.6,0.7,1],
            x_title=f"raw BDT score for fake",
        )
    


def add_highlevel_features(cfg: od.Config) -> None:    
    """
    Adds MET and other high-level features
    """
    bin_factor = 1 # set to 1/2 for broader bins
    cfg.add_variable(
        name="met",
        expression="MET.pt",
        null_value=EMPTY_FLOAT,
        binning=(20, 0.0, 100.0),
        x_title=r"MET",
    )
    cfg.add_variable(
        name="puppi_met_pt_no_jec",
        expression="PuppiMET.pt_no_jec",
        null_value=EMPTY_FLOAT,
        binning=(50, 0,100),
        unit="GeV",
        x_title=r"Uncorrected PuppiMET $p_T$",
    )
    cfg.add_variable(
        name="puppi_met_phi_no_jec",
        expression="PuppiMET.phi_no_jec",
        null_value=EMPTY_FLOAT,
        binning=(32, -3.2,3.2),
        x_title=r"Uncorrected PuppiMET $\phi$",
    )
    cfg.add_variable(
        name="puppi_met_pt",
        expression="PuppiMET.pt",
        null_value=EMPTY_FLOAT,
        binning=(bin_factor*50, 0,100),
        unit="GeV",
        x_title=r"PUPPI MET $p_T$",
    )
    
    cfg.add_variable(
        name="puppi_met_pt_coarse",
        expression="PuppiMET.pt",
        null_value=EMPTY_FLOAT,
        binning=(40, 0,200),
        unit="GeV",
        x_title=r"PUPPI MET $p_T$",
    )
    
    cfg.add_variable(
        name="puppi_met_phi",
        expression="PuppiMET.phi",
        null_value=EMPTY_FLOAT,
        binning=(bin_factor*32, -3.2,3.2),
        x_title=r"PUPPI MET $\phi$",
    ) 
    cfg.add_variable(
        name="puppi_met_pt_no_recoil",
        expression="PuppiMET.pt_no_recoil",
        null_value=EMPTY_FLOAT,
        binning=(50, 0,100),
        unit="GeV",
        x_title=r"PUPPI MET $p_T$ (no recoil corr)",
    )
    cfg.add_variable(
        name="puppi_met_phi_no_recoil",
        expression="PuppiMET.phi_no_recoil",
        null_value=EMPTY_FLOAT,
        binning=(32, -3.2,3.2),
        x_title=r"PUPPI MET $\phi$ (no recoil corr)",
    ) 
    cfg.add_variable(
        name="puppi_met_recoil_effect",
        expression="PuppiMET.pt_recoil_effect",
        null_value=EMPTY_FLOAT,
        binning=(50, -10,10),
        unit="GeV",
        x_title=r"PUPPI MET $p_T^{recoil corr}$ - $p_T^{no corr}$ ",
    )
    cfg.add_variable(
        name="pion_E_split",
        expression="pion_E_split",
        null_value=EMPTY_FLOAT,
        binning=(50, 0,1.2),
        unit="",
        x_title=r"$\frac{|E(\pi^{\pm}) - E(\pi^{0})|}{E(\pi^{\pm}) + E(\pi^{0})}$",
    ) 
    
    
    
    

def add_weight_features(cfg: od.Config) -> None:
    """
    Adds weights
    """
    cfg.add_variable(
        name="mc_weight",
        expression="mc_weight",
        binning=(20, -2, 2),
        x_title="MC weight",
    )
    cfg.add_variable(
        name="pu_weight",
        expression="pu_weight",
        null_value=EMPTY_FLOAT,
        binning=(30, 0,3),
        unit="",
        x_title=r"Pileup weight",
    )
    
    cfg.add_variable(
        name="muon_weight",
        expression="muon_weight_nom",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.5,1.5),
        unit="",
        x_title=r"muon weight",
    )
    
    cfg.add_variable(
        name="tau_weight",
        expression="tau_weight_nom",
        null_value=EMPTY_FLOAT,
        binning=(50, 0.5,1.5),
        unit="",
        x_title=r"tau weight",
    )
    
    for var in ["0", "0p25", "0p375", "0p5", "minus0p25"]:
        
        angle = float(var.replace("minus","-").replace("p", "."))*180
        cfg.add_variable(
        name=f"TauSpinner_weight_cp_{var}",
        expression=f"TauSpinner.weight_cp_{var}",
        null_value=EMPTY_FLOAT,
        binning=(60, -3,3),
        unit="",
        x_title=fr"Tau spinner weight $\Delta \phi$=${angle}^{{\circ}}$",
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


def phi_cp_variables(cfg: od.Config) -> None:
    n_bins_phi_cp = 8
    cfg.add_variable(
        name="phi_cp_incl",
        expression="phi_cp_incl",
        null_value=EMPTY_FLOAT,
        binning=(n_bins_phi_cp, 0, 2*np.pi),
        x_title=rf"$\varphi_{{CP}} (rad)",
    )
    for the_ch in ['mu_pi','mu_pi_gen','mu_rho','mu_rho_gen',
                   'mu_a1_1pr','mu_a1_1pr_gen','mu_a1_3pr_dp','mu_a1_3pr_dp_gen',
                   'mu_a1_3pr_pv','mu_a1_3pr_pv_mtt','mu_a1_3pr_pv_gef','mu_a1_3pr_pv_gen']:

        if the_ch.startswith("mu_a1_3pr"):
            suffix = the_ch.split("_")
            if len(suffix)==4: suffix.append('reco')
            title_str = fr"\mu a_1, 3pr ({suffix[-2]} {suffix[-1]})"

        elif the_ch.startswith("mu_a1_1pr"):
            title_str = r"\mu a_1, 1pr" + (" (gen)" if the_ch.endswith("_gen") else "")

        elif the_ch.startswith(("mu_rho", "mu_pi")):
            had = the_ch.split("_")[1]
            title_str = fr"\mu \{had}" + (" (gen)" if the_ch.endswith("_gen") else "")


        if the_ch.startswith('mu_rho') or the_ch.startswith('mu_a1_3pr'):
            n_bins_phi_cp = 10
        else: 
            n_bins_phi_cp = 8
        cfg.add_variable(
            name=f"phi_cp_{the_ch}",
            expression=f"phi_cp_{the_ch}",
            null_value=EMPTY_FLOAT,
            binning=(n_bins_phi_cp, 0, 2*np.pi),
            x_title=rf"$\varphi_{{CP}} [{title_str}]$ (rad)",
        )
        # cfg.add_variable(
        #     name=f"phi_cp_{the_ch}_reg1",
        #     expression=f"phi_cp_{the_ch}_reg1",
        #     null_value=EMPTY_FLOAT,
        #     binning=(20, 0, 2*np.pi),
        #     x_title=rf"$\varphi_{{CP}} [{title_str}], \alpha < \pi/4$ (rad)",
        # )
        # cfg.add_variable(
        #     name=f"phi_cp_{the_ch}_reg2",
        #     expression=f"phi_cp_{the_ch}_reg2",
        #     null_value=EMPTY_FLOAT,
        #     binning=(20, 0, 2*np.pi),
        #     x_title=rf"$\varphi_{{CP}} [{title_str}], \alpha \geq \pi/4$ (rad)",
        # )
        # # 2-bin histograms
        # cfg.add_variable(
        #     name=f"phi_cp_{the_ch}_2bin",
        #     expression=f"phi_cp_{the_ch}_2bin",
        #     null_value=EMPTY_FLOAT,
        #     binning=(2, 0, 2*np.pi), 
        #     x_title=rf"$\varphi_{{CP}} [{title_str}]$ (rad)",
        # )
        # cfg.add_variable(
        #     name=f"phi_cp_{the_ch}_reg1_2bin",
        #     expression=f"phi_cp_{the_ch}_reg1_2bin",
        #     null_value=EMPTY_FLOAT,
        #     binning=(2, 0, 2*np.pi),
        #     x_title=rf"$\varphi_{{CP}} [{title_str}], \alpha < \pi/4$ (rad)",
        # )
        # cfg.add_variable(
        #     name=f"phi_cp_{the_ch}_reg2_2bin",
        #     expression=f"phi_cp_{the_ch}_reg2_2bin",
        #     null_value=EMPTY_FLOAT,
        #     binning=(2, 0, 2*np.pi),
        #     x_title=rf"$\varphi_{{CP}} [{title_str}], \alpha \geq \pi/4$ (rad)",
        # )
        # cfg.add_variable(
        #     name=f"alpha_{the_ch}",
        #     expression=f"alpha_{the_ch}",
        #     null_value=EMPTY_FLOAT,
        #     binning=(6, 0, np.pi/2),
        #     x_title=rf"$ \alpha [{title_str}] $(rad)",
        # )

def add_dilepton_features(cfg: od.Config) -> None:
    bin_factor = 1/2
    channels = cfg.channels.names()
    ch_objects = DotDict.wrap({
        'etau'  : {'lep0':'Electron',
                   'lep1':'Tau'    },
        'mutau' : {'lep0':'Muon'    ,
                   'lep1':'Tau'    },
        'emu'   : {'lep0':'Electron',
                   'lep1':'Muon'   },
        'tautau': {'lep0':'Tau'     ,
                   'lep1':'Tau'    },
    })
    bin_split_factor = 4 #Used to define histograms for kinematic variables with finer binning
    for ch_str in channels:
        lep0_str = ch_str.replace('tau','')
        
        cfg.add_variable(
                name=f"{ch_str}_mvis",
                expression=f"hcand_{ch_str}.mass",
                null_value=EMPTY_FLOAT,
                binning=(bin_factor*40, 0.0, 200.0),
                unit="GeV",
                x_title=r"$m_{vis}$",
            )
        cfg.add_variable(
                name=f"{ch_str}_mvis_coarse",
                expression=f"hcand_{ch_str}.mass",
                null_value=EMPTY_FLOAT,
                binning=(50, 0.0, 200.0),
                unit="GeV",
                x_title=r"$m_{vis}$",
            )
        cfg.add_variable(
                name=f"{ch_str}_fastMTT_mass",
                expression=f"hcand_{ch_str}.fastMTT.mass",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=r"$m_{FastMTT}$",
            )
        cfg.add_variable(
                name=f"{ch_str}_mvis_fine",
                expression=f"hcand_{ch_str}.mass",
                null_value=EMPTY_FLOAT,
                binning=(115, 15, 130),
                unit="GeV",
                x_title=r"$m_{vis}$",
            )
        if ch_str != 'tautau':
            cfg.add_variable(
                name=f"{ch_str}_mt0",
                expression=f"hcand_{ch_str}.mt_0",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=rf"$m_{{T}}(\{lep0_str},\mathrm{{MET}})$",
            )
            cfg.add_variable(
                name=f"{ch_str}_mt1",
                expression=f"hcand_{ch_str}.mt_1",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=r"$m_{T}(\tau,\mathrm{MET})$",
            )
            cfg.add_variable(
                name=f"{ch_str}_mt_ll",
                expression=f"hcand_{ch_str}.mt_ll",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=r"$m_{T}(\ell\ell,\mathrm{MET})$",
            )
            cfg.add_variable(
                name=f"{ch_str}_mt_vis",
                expression=f"hcand_{ch_str}.mt_vis",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=rf"$m_{{T}}(\{lep0_str},\tau)$",
            )
            cfg.add_variable(
                name=f"{ch_str}_mt_tot",
                expression=f"hcand_{ch_str}.mt_tot",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=r"$m_{T}(\mathrm{tot})$",
            )
            cfg.add_variable(
                name=f"{ch_str}_pt_vis",
                expression=f"hcand_{ch_str}.pt_vis",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=r"$p_{T\mathrm{vis}}$",
            )
            
            cfg.add_variable(
                name=f"{ch_str}_pt_vis_coarse",
                expression=f"hcand_{ch_str}.pt_vis",
                null_value=EMPTY_FLOAT,
                binning=(50, 0.0, 400.0),
                unit="GeV",
                x_title=r"$p_{T\mathrm{vis}}$",
            )
             
            cfg.add_variable(
                name=f"{ch_str}_pt_ll",
                expression=f"hcand_{ch_str}.pt_ll",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=r"$p_{T}(\ell\ell)}$",
            )
            
        cfg.add_variable(
            name=f"{ch_str}_delta_r",
            expression=f"hcand_{ch_str}.delta_r",
            null_value=EMPTY_FLOAT,
            binning=(bin_factor*40, 0, 4),
            x_title=r"$\Delta R(\ell,\ell)$",
        )
        cfg.add_variable(
            name=f"{ch_str}_delta_eta",
            expression=f"hcand_{ch_str}.delta_eta",
            null_value=EMPTY_FLOAT,
            binning=(40, -5, 5),
            x_title=r"$\Delta \eta(\ell,\ell)$",
        )
        cfg.add_variable(
            name=f"{ch_str}_pt",
            expression=f"hcand_{ch_str}.pt",
            null_value=EMPTY_FLOAT,
            binning=(bin_factor*40, 0.0, 200.0),
            unit="GeV/c",
            x_title=r"$p_{T}(\ell\ell)$",
        )
        cfg.add_variable(
            name=f"{ch_str}_delta_phi_0_met",
            expression=f"hcand_{ch_str}.delta_phi_0_met",
            null_value=EMPTY_FLOAT,
            binning=(32, -3.3, 3.3),
            x_title=rf"$\Delta \phi(\{lep0_str},\mathrm{{MET}})$",
        ) 
        cfg.add_variable(
            name=f"{ch_str}_delta_phi_1_met",
            expression=f"hcand_{ch_str}.delta_phi_1_met",
            null_value=EMPTY_FLOAT,
            binning=(32, -3.3, 3.3),
            x_title=rf"$\Delta \phi(\tau,\mathrm{{MET}})$",
        ) 
        cfg.add_variable(
            name=f"{ch_str}_delta_phi",
            expression=f"hcand_{ch_str}.delta_phi",
            null_value=EMPTY_FLOAT,
            binning=(32, -3.3, 3.3),
            x_title=rf"$\Delta \phi(\{lep0_str},\tau)$",
        ) 
        
        for lep in ['lep0','lep1']:
            if ch_str != 'tautau': lep_str = ch_objects[ch_str][lep].lower()
            else: lep_str = f'tau {lep[3:]}'
            cfg.add_variable(
                name=f"{ch_str}_{lep}_pt",
                expression=f"hcand_{ch_str}.{lep}.pt",
                null_value=EMPTY_FLOAT,
                binning=(bin_factor*40, 20, 100.),
                unit="GeV",
                x_title= rf"{lep_str} $p_{{T}}$",
            )
            
            cfg.add_variable(
                name=f"{ch_str}_{lep}_pt_coarse",
                expression=f"hcand_{ch_str}.{lep}.pt",
                null_value=EMPTY_FLOAT,
                binning=(40, 0, 200.),
                unit="GeV",
                x_title= rf"{lep_str} $p_{{T}}$",
            )
            
            cfg.add_variable(
                name=f"{ch_str}_{lep}_eta",
                expression=f"hcand_{ch_str}.{lep}.eta",
                null_value=EMPTY_FLOAT,
                binning=(30, -3, 3),
                x_title=rf"{lep_str} $\eta$",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_phi",
                expression=f"hcand_{ch_str}.{lep}.phi",
                null_value=EMPTY_FLOAT,
                binning=(32, -3.3, 3.3),
                x_title=rf"{lep_str} $\phi$",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_mass",
                expression=f"hcand_{ch_str}.{lep}.mass",
                null_value=EMPTY_FLOAT,
                binning=(bin_factor*30, 0, 3),
                unit="GeV",
                x_title=f"{lep_str} mass",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_decayModePNet",
                expression=f"hcand_{ch_str}.{lep}.decayModePNet",
                null_value=EMPTY_FLOAT,
                binning=(12,0,12),
                unit="",
                x_title=rf"{lep_str} PNet decay mode",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_decayMode",
                expression=f"hcand_{ch_str}.{lep}.decayMode",
                null_value=EMPTY_FLOAT,
                binning=(12,0,12),
                unit="",
                x_title=rf"{lep_str} HPS decay mode",
            )
            cfg.add_variable(
                name='_'.join((ch_str,lep,'ip_sig', 'nom')),
                expression=f"hcand_{ch_str}.{lep}.ip_sig_nom",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 10),
                unit="",
                x_title= rf"{lep_str} $\frac{{|IP|}}{{\sigma(IP)}}$ nom",
            )
            cfg.add_variable(
                name='_'.join((ch_str,lep,'ip_sig')),
                expression=f"hcand_{ch_str}.{lep}.ip_sig",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 10),
                unit="",
                x_title= rf"{lep_str} $\frac{{|IP|}}{{\sigma(IP)}}$",
            )
            for coord in ['x','y','z']:
                cfg.add_variable(
                    name=f"{lep}_IP{coord}",
                    expression=f"hcand_{ch_str}.{lep}.IP{coord}",
                    null_value=EMPTY_FLOAT,
                    binning=(30, -0.015, 0.015),
                    unit="",
                    x_title= rf"{lep_str} $IP_{coord}$",
                )
                if lep == 'lep1':
                    if coord == 'z':
                        cfg.add_variable(
                            name=f"PV_{coord}",
                            expression=f"PV.{coord}",
                            null_value=EMPTY_FLOAT,
                            binning=(40, -20.0, 20.0),
                            unit="cm",
                            x_title= rf"$PV_{coord}$",
                        )
                        cfg.add_variable(
                            name=f"PVBS_{coord}",
                            expression=f"PVBS.{coord}",
                            null_value=EMPTY_FLOAT,
                            binning=(40, -20.0, 20.0),
                            unit="cm",
                            x_title= rf"$PV_{coord}$ (beam-spot)",
                        )
                        cfg.add_variable(
                            name=f"{ch_str}_{lep}_refitSV{coord}",
                            expression=f"hcand_{ch_str}.{lep}.refitSV{coord}",
                            null_value=EMPTY_FLOAT,
                            binning=(40, -20.0, 20.0),
                            unit="cm",
                            x_title= rf"{lep_str} $refited SV_{coord}$",
                        )
                    else:
                        cfg.add_variable(
                            name=f"PV_{coord}",
                            expression=f"PV.{coord}",
                            null_value=EMPTY_FLOAT,
                            binning=(40, -0.25, 0.25),
                            unit="cm",
                            x_title= rf"$PV_{coord}$",
                        )
                        cfg.add_variable(
                            name=f"PVBS_{coord}",
                            expression=f"PVBS.{coord}",
                            null_value=EMPTY_FLOAT,
                            binning=(40, -0.25, 0.25),
                            unit="cm",
                            x_title= rf"$PV_{coord}$ (beam-spot)",
                        )
                        cfg.add_variable(
                            name=f"{ch_str}_{lep}_refitSV{coord}",
                            expression=f"hcand_{ch_str}.{lep}.refitSV{coord}",
                            null_value=EMPTY_FLOAT,
                            binning=(40, -1.0, 1.0),
                            unit="cm",
                            x_title= rf"{lep_str} $refited SV_{coord}$",
                        )
            if lep == 'lep0':
                if ch_str == 'mutau':
                    cfg.add_variable(
                        name=f"{ch_str}_{lep}_iso",
                        expression=f"hcand_{ch_str}.{lep}.pfRelIso04_all",
                        null_value=EMPTY_FLOAT,
                        binning=(40,0, 0.3),
                        x_title=rf"{lep_str} rel. iso",
                    )
                    cfg.add_variable(
                        name=f"{ch_str}_{lep}_iso_full_range",
                        expression=f"hcand_{ch_str}.{lep}.pfRelIso04_all",
                        null_value=EMPTY_FLOAT,
                        binning=(30,0, 0.6),
                        x_title=rf"{lep_str} rel. iso",
                    )
                elif ch_str == 'etau':
                    cfg.add_variable(
                        name=f"{ch_str}_{lep}_iso",
                        expression=f"hcand_{ch_str}.{lep}.pfRelIso03_all",
                        null_value=EMPTY_FLOAT,
                        binning=(40,0,0.6),
                        x_title=rf"{lep_str} rel. iso",
                    )
                    cfg.add_variable(
                        name=f"{ch_str}_{lep}_iso_full_range",
                        expression=f"hcand_{ch_str}.{lep}.pfRelIso03_all",
                        null_value=EMPTY_FLOAT,
                        binning=(40,0,2.4),
                        x_title=rf"{lep_str} rel. iso",
                    )
            
            for proj in ['x','y','z']:
                cfg.add_variable(
                    name=f"{lep}_IP{proj}_qm",
                    expression=f"hcand_{ch_str}.{lep}.IP{proj}_qm",
                    null_value=EMPTY_FLOAT,
                    binning=(30, -0.015, 0.015),
                    unit="",
                    x_title= rf"{lep_str} $IP_{proj} corrected$",
                )

            #Variables with finer binning
            cfg.add_variable(
                name=f"{ch_str}_{lep}_pt_fine_binning",
                expression=f"hcand_{ch_str}.{lep}.pt",
                null_value=EMPTY_FLOAT,
                binning=(30*bin_split_factor, 20, 80.),
                unit="GeV",
                x_title= rf"{lep_str} $p_{{T}}$",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_eta_fine_binning",
                expression=f"hcand_{ch_str}.{lep}.eta",
                null_value=EMPTY_FLOAT,
                binning=(32*bin_split_factor, -3.2, 3.2),
                x_title=rf"{lep_str} $\eta$",
            )
            cfg.add_variable(
                name=f"{ch_str}_{lep}_phi_fine_binning",
                expression=f"hcand_{ch_str}.{lep}.phi",
                null_value=EMPTY_FLOAT,
                binning=(32*bin_split_factor, -3.2, 3.2),
                x_title=rf"{lep_str} $\phi$",
            )

            ## FastMTT variables

            cfg.add_variable(
                name=f"hcand_{ch_str}_fastMTT_{lep}_px",
                expression=f"hcand_{ch_str}.fastMTT.{lep}.px",
                null_value=EMPTY_FLOAT,
                binning=(42, -10., 200.),
                unit="GeV",
                x_title=f"{lep} " + r"$p_{x}^{fastMTT}$",
            )
            cfg.add_variable(
                name=f"hcand_{ch_str}_fastMTT_{lep}_py",
                expression=f"hcand_{ch_str}.fastMTT.{lep}.py",
                null_value=EMPTY_FLOAT,
                binning=(42, -10., 200.),
                unit="GeV",
                x_title=f"{lep} " + r"$p_{y}^{fastMTT}$",
            )
            cfg.add_variable(
                name=f"hcand_{ch_str}_fastMTT_{lep}_pz",
                expression=f"hcand_{ch_str}.fastMTT.{lep}.pz",
                null_value=EMPTY_FLOAT,
                binning=(42, -10., 200.),
                unit="GeV",
                x_title=f"{lep} " + r"$p_{z}^{fastMTT}$",
            )
            cfg.add_variable(
                name=f"hcand_{ch_str}_fastMTT_{lep}_pt",
                expression=f"hcand_{ch_str}.fastMTT.{lep}.pt",
                null_value=EMPTY_FLOAT,
                binning=(40, 0., 200.),
                unit="GeV",
                x_title=f"{lep} " + r"$p_{T}^{fastMTT}$",
            )
            cfg.add_variable(
                name=f"hcand_{ch_str}_fastMTT_{lep}_eta",
                expression=f"hcand_{ch_str}.fastMTT.{lep}.eta",
                null_value=EMPTY_FLOAT,
                binning=(25, -3.0, 3.0),
                unit="GeV",
                x_title=f"{lep} " + r"$\eta^{fastMTT}$",
            )
            cfg.add_variable(
                name=f"hcand_{ch_str}_fastMTT_{lep}_phi",
                expression=f"hcand_{ch_str}.fastMTT.{lep}.phi",
                null_value=EMPTY_FLOAT,
                binning=(32, -3.2, 3.2),
                unit="GeV",
                x_title=f"{lep} " + r"$\phi^{fastMTT}$",
            )
            cfg.add_variable(
                name=f"hcand_{ch_str}_fastMTT_{lep}_mass",
                expression=f"hcand_{ch_str}.fastMTT.{lep}.mass",
                null_value=EMPTY_FLOAT,
                binning=(50, 0.01, 3.0),
                unit="GeV",
                x_title=f"{lep} " + r"$m^{fastMTT}$",
            )
        cfg.add_variable(
                name=f"hcand_{ch_str}_fastMTT_mass",
                expression=f"hcand_{ch_str}.fastMTT.mass",
                null_value=EMPTY_FLOAT,
                binning=(40, 0.0, 200.0),
                unit="GeV",
                x_title=r"$mass^{fastMTT}$",
        )
            
            ## FastMTT variables
            pretty_label = {
                "": "unconstrained",
                "BW": "mass-window",
                "cons": "Breit Wigner",
            }
            for prefix in ['', '_BW', '_cons']:
                key = prefix.strip("_")
                label = pretty_label[key]

                cfg.add_variable(
                    name=f"hcand_{ch_str}_fastMTT{prefix}_{lep}_px",
                    expression=f"hcand_{ch_str}.fastMTT{prefix}.{lep}.px",
                    null_value=EMPTY_FLOAT,
                    binning=(42, -10., 200.),
                    unit="GeV",
                    x_title=f"{lep} " + rf"$p_{{x}}$ (fastMTT: {label})"
                )
                cfg.add_variable(
                    name=f"hcand_{ch_str}_fastMTT{prefix}_{lep}_py",
                    expression=f"hcand_{ch_str}.fastMTT{prefix}.{lep}.py",
                    null_value=EMPTY_FLOAT,
                    binning=(42, -10., 200.),
                    unit="GeV",
                    x_title=f"{lep} " + rf"$p_{{y}}$ (fastMTT: {label})"
                )
                cfg.add_variable(
                    name=f"hcand_{ch_str}_fastMTT{prefix}_{lep}_pz",
                    expression=f"hcand_{ch_str}.fastMTT{prefix}.{lep}.pz",
                    null_value=EMPTY_FLOAT,
                    binning=(42, -10., 200.),
                    unit="GeV",
                    x_title=f"{lep} " + rf"$p_{{z}}$ (fastMTT: {label})"
                )
                cfg.add_variable(
                    name=f"hcand_{ch_str}_fastMTT{prefix}_{lep}_pt",
                    expression=f"hcand_{ch_str}.fastMTT{prefix}.{lep}.pt",
                    null_value=EMPTY_FLOAT,
                    binning=(bin_factor*40, 0., 200.),
                    unit="GeV",
                    x_title=f"{lep} " + rf"$p_{{T}}$ (fastMTT: {label})"
                )
                cfg.add_variable(
                    name=f"hcand_{ch_str}_fastMTT{prefix}_{lep}_eta",
                    expression=f"hcand_{ch_str}.fastMTT{prefix}.{lep}.eta",
                    null_value=EMPTY_FLOAT,
                    binning=(25, -3.0, 3.0),
                    unit="GeV",
                    x_title=f"{lep} " + rf"$\eta$ (fastMTT: {label})"
                )
                cfg.add_variable(
                    name=f"hcand_{ch_str}_fastMTT{prefix}_{lep}_phi",
                    expression=f"hcand_{ch_str}.fastMTT{prefix}.{lep}.phi",
                    null_value=EMPTY_FLOAT,
                    binning=(32, -3.2, 3.2),
                    unit="GeV",
                    x_title=f"{lep} " + rf"$\phi$ (fastMTT: {label})"
                )
                cfg.add_variable(
                    name=f"hcand_{ch_str}_fastMTT{prefix}_{lep}_mass",
                    expression=f"hcand_{ch_str}.fastMTT{prefix}.{lep}.mass",
                    null_value=EMPTY_FLOAT,
                    binning=(50, 0.01, 3.0),
                    unit="GeV",
                    x_title=f"{lep} " + rf"$m$ (fastMTT: {label})"
               )
        pretty_label = {
                "": "unconstrained",
                "BW": "mass-window",
                "cons": "Breit Wigner",
            }
        for prefix in ['', '_BW', '_cons']:
            key = prefix.strip("_")
            label = pretty_label[key]

            cfg.add_variable(
               name=f"hcand_{ch_str}_fastMTT{prefix}_mass",
               expression=f"hcand_{ch_str}.fastMTT{prefix}.mass",
               null_value=EMPTY_FLOAT,
               binning=(bin_factor*40, 0.0, 200.0),
               unit="GeV",
               x_title=rf"$m_{{\text{{fastMTT}}}}$ ({label})",
            )
        cfg.add_variable(
            name=f"theta_gj_mu_a1_3pr_pv_gef",
            expression=f"theta_gj_mu_a1_3pr_pv_gef",
            null_value=EMPTY_FLOAT,
            binning=(40, 0, 0.1),
            x_title=r"$\theta_{GJ}$",
        )
        cfg.add_variable(
            name=f"theta_max_mu_a1_3pr_pv_gef",
            expression=f"theta_max_mu_a1_3pr_pv_gef",
            null_value=EMPTY_FLOAT,
            binning=(40, 0, 0.1),
            x_title=r"$\theta_{GJ}^{max}$",
        )
        cfg.add_variable(
            name=f"theta_rot_mu_a1_3pr_pv_gef",
            expression=f"theta_rot_mu_a1_3pr_pv_gef",
            null_value=EMPTY_FLOAT,
            binning=(40, 0, 0.1),
            x_title=r"$\theta_{GJ}^{rotated}$",
        )
        for coord in ['x','y','z']:
            cfg.add_variable(
                name=f"dsvpv_{coord}",
                expression=f"dsvpv_{coord}",
                null_value=EMPTY_FLOAT,
                binning=(40, -1.2, 1.2),
                unit="cm",
                x_title=rf"$d(SV^{{\text{{refit}}}},PV)_{coord}$",
            )
        cfg.add_variable(
            name=f"dsvpv_mag",
            expression=f"dsvpv_mag",
            null_value=EMPTY_FLOAT,
            binning=(40, 0, 1.5),
            unit="cm",
            x_title=r"$|d(SV^{refit},PV)|$",
        )
        
        
        



def add_genlep_features(cfg: od.Config) -> None:
    bin_factor = 1/2
    channels = cfg.channels.names()
    ch_objects = DotDict.wrap({
        'etau'  : {'lep0':'Electron',
                   'lep1':'Tau'    },
        'mutau' : {'lep0':'Muon'    ,
                   'lep1':'Tau'    },
        'emu'   : {'lep0':'Electron',
                   'lep1':'Muon'   },
        'tautau': {'lep0':'Tau'     ,
                   'lep1':'Tau'    },
    })
    for ch_str in channels:
        for lep in ['lep0','lep1']:
            if ch_str != 'tautau': lep_str = ch_objects[ch_str][lep].lower()
            cfg.add_variable(
                name=f"gen_{lep}_pt",
                expression=f"gen_lep.{lep}.pt",
                null_value=EMPTY_FLOAT,
                binning=(bin_factor*40, 20, 100.),
                unit="GeV",
                x_title= rf"Gen {lep_str} $p_{{T}}$",
            )
            cfg.add_variable(
                name=f"gen_{lep}_pt_coarse",
                expression=f"gen_lep.{lep}.pt",
                null_value=EMPTY_FLOAT,
                binning=(40, 0, 200.),
                unit="GeV",
                x_title= rf"Gen {lep_str} $p_{{T}}$",
            )
            cfg.add_variable(
                name=f"gen_{lep}_eta",
                expression=f"gen_lep.{lep}.eta",
                null_value=EMPTY_FLOAT,
                binning=(30, -3, 3),
                x_title=rf"Gen {lep_str} $\eta$",
            )
            cfg.add_variable(
                name=f"gen_{lep}_phi",
                expression=f"gen_lep.{lep}.phi",
                null_value=EMPTY_FLOAT,
                binning=(32, -3.3, 3.3),
                x_title=rf"Gen {lep_str} $\phi$",
            )
            for coord in ['x','y','z']:
                cfg.add_variable(
                    name=f"gen_{lep}_IP{coord}",
                    expression=f"gen_lep.{lep}.IP{coord}",
                    null_value=EMPTY_FLOAT,
                    binning=(40, -0.004, 0.004),
                    unit="cm",
                    x_title= rf"Gen {lep_str} $IP_{coord}$",
                )

def add_ff_variables(cfg: od.Config) -> None:
    for the_name, ax in cfg.x.fake_factor_method.axes.items():
        cfg.add_variable(name = the_name, expression = '.'.join(ax.var_route))

def add_variables(cfg: od.Config) -> None:
    """
    Adds all variables to a *config*.
    """
    add_common_features(cfg)
    add_lepton_features(cfg)
    add_jet_features(cfg)
    add_highlevel_features(cfg)
    phi_cp_variables(cfg)
    add_weight_features(cfg)
    add_cutflow_features(cfg)
    add_dilepton_features(cfg)
    add_genlep_features(cfg)
    add_hcp_bdt_output(cfg)
    add_ff_variables(cfg)
   