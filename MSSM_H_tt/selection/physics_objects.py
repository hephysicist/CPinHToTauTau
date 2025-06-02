# coding: utf-8

"""
Exemplary selection methods.
"""

from collections import defaultdict

from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.util import sorted_indices_from_mask
from columnflow.util import maybe_import, DotDict
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional

from MSSM_H_tt.util import IF_NANO_V9, IF_NANO_V11
from MSSM_H_tt.util import getGenTauDecayMode

np = maybe_import("numpy")
ak = maybe_import("awkward")


# ------------------------------------------------------------------------------------------------------- #
# Muon Selection
# Reference:
#   https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
#   http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
# ------------------------------------------------------------------------------------------------------- #
@selector(
    uses={
        # Muon nano columns
        f"Muon.{var}" for var in [
            "pt", "eta", "phi", "dxy", "dz", "mediumId", "looseId", 
            "pfRelIso04_all", "isGlobal", "isPFcand", "pdgId",'sip3d',
            #"isTracker",
        ]
    },
    produces={
        f"Muon.{var}" for var in [
            "rawIdx","ip_sig"
        ]
    },
    exposed=False,
)
def muon_selection(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult, ak.Array, ak.Array, ak.Array]:
    """
    Muon selection returning two sets of indidces for default and veto muons.
    
    References:
      - Isolation working point: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=59
      - ID und ISO : https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017?rev=15
    """
    # setting two new columns for the muons
    events = set_ak_column(events, "Muon.rawIdx",    ak.local_index(events.Muon))
    events = set_ak_column(events, "Muon.ip_sig", events.Muon.sip3d)

    good_selections = {
        "muon_pt_26"          : events.Muon.pt > 26,
        "muon_eta_2p4"        : abs(events.Muon.eta) < 2.4,
        "mediumID"            : events.Muon.mediumId == 1,
        "muon_dxy_0p045"      : abs(events.Muon.dxy) < 0.045,
        "muon_dz_0p2"         : abs(events.Muon.dz) < 0.2,
        #"muon_iso_0p15"       : events.Muon.pfRelIso04_all < 0.15
    }
    mu_emu_selections = {
        "muon_pt_15"          : events.Muon.pt > 15,
        "muon_eta_2p4"        : abs(events.Muon.eta) < 2.4,
        "mediumID"            : events.Muon.mediumId == 1,
        "muon_dxy_0p045"      : abs(events.Muon.dxy) < 0.045,
        "muon_dz_0p2"         : abs(events.Muon.dz) < 0.2,
#        "muon_iso_0p2"        : events.Muon.pfRelIso04_all < 0.2, remove the cut, need it to define ABCD regions
    }
    single_veto_selections = {
        "muon_pt_10"          : events.Muon.pt > 10,
        "muon_eta_2p4"        : abs(events.Muon.eta) < 2.4,
        "mediumID"            : events.Muon.mediumId == 1,
        "muon_iso_0p3"        : events.Muon.pfRelIso04_all < 0.3
    }
    second_mu_veto_emu_selections = {
        "muon_pt_10"          : events.Muon.pt > 10,
        "muon_eta_2p4"        : abs(events.Muon.eta) < 2.4,
        "mediumID"            : events.Muon.mediumId == 1,
        "muon_dxy_0p045"      : abs(events.Muon.dxy) < 0.045,
        "muon_dz_0p2"         : abs(events.Muon.dz) < 0.2,
        "muon_iso_0p3"        : events.Muon.pfRelIso04_all < 0.3
    }
    OC_veto_selections = {
        "muon_pt_15"          : events.Muon.pt > 15,
        "muon_eta_2p4"        : abs(events.Muon.eta) < 2.4,
        "muon_dxy_0p045"      : abs(events.Muon.dxy) < 0.045,
        "muon_dz_0p2"         : abs(events.Muon.dz) < 0.2,
        "looseID"             : events.Muon.looseId == 1,
        "muon_iso_0p3"        : events.Muon.pfRelIso04_all < 0.3
    }

    # pt sorted indices for converting masks to indices
    sorted_indices          = ak.argsort(events.Muon.pt, axis=-1, ascending=False)    
    muon_mask               = ak.local_index(events.Muon.pt) >= 0
    good_muon_mask          = muon_mask
    mu_emu_mask             = muon_mask
    single_veto_muon_mask   = muon_mask
    second_mu_veto_emu_mask = muon_mask
    OC_veto_muon_mask       = muon_mask
    selection_steps = {}

    selection_steps = {"Starts with": good_muon_mask}
    for cut in good_selections.keys():
        good_muon_mask = good_muon_mask & ak.fill_none(good_selections[cut], False)
    
    for cut in mu_emu_selections.keys():
        mu_emu_mask = mu_emu_mask & ak.fill_none(mu_emu_selections[cut], False)
        selection_steps[cut] = mu_emu_mask
        
    for cut in single_veto_selections.keys():
        single_veto_muon_mask = single_veto_muon_mask & ak.fill_none(single_veto_selections[cut], False)
        
    for cut in second_mu_veto_emu_selections.keys():
        second_mu_veto_emu_mask = second_mu_veto_emu_mask & ak.fill_none(second_mu_veto_emu_selections[cut], False) 
        
    for cut in OC_veto_selections.keys():
        OC_veto_muon_mask = OC_veto_muon_mask & ak.fill_none(OC_veto_selections[cut], False)
        
    # convert to sorted indices
    good_muon_indices          = sorted_indices[good_muon_mask[sorted_indices]]
    good_muon_indices          = ak.values_astype(good_muon_indices, np.int32)
    
    mu_emu_indices             = sorted_indices[mu_emu_mask[sorted_indices]]
    mu_emu_indices             = ak.values_astype(mu_emu_indices, np.int32)

    single_veto_muon_indices   = sorted_indices[single_veto_muon_mask[sorted_indices]]
    single_veto_muon_indices   = ak.values_astype(single_veto_muon_indices, np.int32)
    
    second_mu_veto_emu_indices = sorted_indices[second_mu_veto_emu_mask[sorted_indices]]
    second_mu_veto_emu_indices = ak.values_astype(second_mu_veto_emu_indices, np.int32)
    
    OC_veto_muon_indices       = sorted_indices[OC_veto_muon_mask[sorted_indices]]
    OC_veto_muon_indices       = ak.values_astype(OC_veto_muon_indices, np.int32)

    return events, SelectionResult(
        objects={
            "Muon": {
                "Muon"             : good_muon_indices,
                "Ele_emu"          : mu_emu_indices,
                "SingleVetoMuon"   : single_veto_muon_indices,
                "SecondMuVeto_emu" : second_mu_veto_emu_indices,
                "OCVetoMuon"       : OC_veto_muon_indices,
            },
        },
        aux=selection_steps,
    ), good_muon_indices, mu_emu_indices, single_veto_muon_indices, second_mu_veto_emu_indices, OC_veto_muon_indices


# ------------------------------------------------------------------------------------------------------- #
# Electron Selection
# Reference:
#   https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
#   http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
# ------------------------------------------------------------------------------------------------------- #
@selector(
    uses={
        # Electron nano columns
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass", "Electron.dxy", "Electron.dz", "Electron.pdgId",
        "Electron.pfRelIso03_all", "Electron.convVeto", #"lostHits",
        "Electron.mvaIso_WP80", "Electron.mvaIso_WP90", "Electron.mvaNoIso_WP90",
        "Electron.cutBased",'Electron.sip3d','Electron.lostHits','Electron.convVeto',
    },
    produces={
        f"Electron.{var}" for var in [
            "rawIdx", "ip_sig"
        ]
    },
    exposed=False,
)
def electron_selection(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult, ak.Array, ak.Array, ak.Array]:
    """
    Electron selection returning two sets of indidces for default and veto muons.
    
    References:
      - https://twiki.cern.ch/twiki/bin/view/CMS/EgammaNanoAOD?rev=4
    """
    # adding two new fields
    events = set_ak_column(events, "Electron.rawIdx",    ak.local_index(events.Electron))
    events = set_ak_column(events, "Electron.ip_sig", events.Electron.sip3d)

    # >= nano v10
    mva_iso_wp80 = events.Electron.mvaIso_WP80
    mva_iso_wp90 = events.Electron.mvaIso_WP90
    mva_noniso_wp90 = events.Electron.mvaNoIso_WP90

    number_of_missing_hits = events.Electron.lostHits
    
    good_selections = {
        "electron_pt_33"          : events.Electron.pt > 33, #We are using HLT_Ele30 so this cut is needed
        "electron_eta_2p1"        : abs(events.Electron.eta) < 2.1,
        "electron_dxy_0p045"      : abs(events.Electron.dxy) < 0.045,
        "electron_dz_0p2"         : abs(events.Electron.dz) < 0.2,
        #"electron_pfRelIso03_all" : events.Electron.pfRelIso03_all < 0.15,
        "electron_mva_iso_wp90"   : mva_iso_wp90 == 1,
    }
    ele_emu_selections = {
        "electron_pt_15"          : events.Electron.pt > 15,
        "electron_eta_2p5"        : abs(events.Electron.eta) < 2.5,
        "electron_dxy_0p045"      : abs(events.Electron.dxy) < 0.045,
        "electron_dz_0p2"         : abs(events.Electron.dz) < 0.2,
        #"electron_pfRelIso03_all" : events.Electron.pfRelIso03_all < 0.15,
        "electron_mva_iso_wp90"   : mva_iso_wp90 == 1,
        "electron_missing_hits"   : number_of_missing_hits < 2,
        "electron_pass_conVeto"   : events.Electron.convVeto == 1,
    }
    single_veto_selections = {
        "electron_pt_10"          : events.Electron.pt > 10,
        "electron_eta_2p5"        : abs(events.Electron.eta) < 2.5,
        "electron_mva_iso_wp90"   : mva_iso_wp90 == 1,
        "electron_pfRelIso03_all" : events.Electron.pfRelIso03_all < 0.3,
    }
    second_ele_veto_emu_selections = {
        "electron_pt_10"          : events.Electron.pt > 10,
        "electron_eta_2p5"        : abs(events.Electron.eta) < 2.5,
        "electron_dxy_0p045"      : abs(events.Electron.dxy) < 0.045,
        "electron_dz_0p2"         : abs(events.Electron.dz) < 0.2,
        "electron_pfRelIso03_all" : events.Electron.pfRelIso03_all < 0.3,
        "electron_mva_noniso_wp90": mva_noniso_wp90 == 1,
        "electron_missing_hits"   : number_of_missing_hits < 2,
        "electron_pass_conVeto"   : events.Electron.convVeto == 1,
    }
    OC_veto_selections = {
        "electron_pt_15"          : events.Electron.pt > 15,
        "electron_eta_2p5"        : abs(events.Electron.eta) < 2.5,
        "electron_dxy_0p045"      : abs(events.Electron.dxy) < 0.045,
        "electron_dz_0p2"         : abs(events.Electron.dz) < 0.2,
        "electron_cutBased"       : events.Electron.cutBased == 1,
        "electron_pfRelIso03_all" : events.Electron.pfRelIso03_all < 0.3,
    }

    # pt sorted indices for converting masks to indices
    sorted_indices            = ak.argsort(events.Electron.pt, axis=-1, ascending=False)    
    electron_mask             = ak.local_index(events.Electron.pt) >= 0
    good_electron_mask        = electron_mask
    ele_emu_mask              = electron_mask
    single_veto_electron_mask = electron_mask
    second_ele_veto_emu_mask  = electron_mask
    OC_veto_electron_mask     = electron_mask
    selection_steps = {}

    selection_steps = {"Starts with": good_electron_mask}
    for cut in good_selections.keys():
        good_electron_mask = good_electron_mask & ak.fill_none(good_selections[cut], False)
        
        
    for cut in ele_emu_selections.keys():
        ele_emu_mask = ele_emu_mask & ak.fill_none(ele_emu_selections[cut], False)
        selection_steps[cut] = ele_emu_mask
        
    for cut in single_veto_selections.keys():
        single_veto_electron_mask = single_veto_electron_mask & ak.fill_none(single_veto_selections[cut], False)
    
    for cut in second_ele_veto_emu_selections.keys():
        second_ele_veto_emu_mask = second_ele_veto_emu_mask & ak.fill_none(second_ele_veto_emu_selections[cut], False)    
    
    for cut in OC_veto_selections.keys():
        OC_veto_electron_mask = OC_veto_electron_mask & ak.fill_none(OC_veto_selections[cut], False)
        
    # convert to sorted indices
    good_electron_indices        = sorted_indices[good_electron_mask[sorted_indices]]
    good_electron_indices        = ak.values_astype(good_electron_indices, np.int32)

    ele_emu_indices              = sorted_indices[ele_emu_mask[sorted_indices]]
    ele_emu_indices              = ak.values_astype(ele_emu_indices, np.int32)

    single_veto_electron_indices = sorted_indices[single_veto_electron_mask[sorted_indices]]
    single_veto_electron_indices = ak.values_astype(single_veto_electron_indices, np.int32)

    second_ele_veto_emu_indices  = sorted_indices[second_ele_veto_emu_mask[sorted_indices]]
    second_ele_veto_emu_indices  = ak.values_astype(second_ele_veto_emu_indices, np.int32)
    
    OC_veto_electron_indices     = sorted_indices[OC_veto_electron_mask[sorted_indices]]
    OC_veto_electron_indices     = ak.values_astype(OC_veto_electron_indices, np.int32)


    return events, SelectionResult(
        objects={
            "Electron": {
                "Electron"          : good_electron_indices,
                "Ele_emu"           : ele_emu_indices,
                "SingleVetoElectron": single_veto_electron_indices,
                "SecondEleVeto_emu" : second_ele_veto_emu_indices,
                "OCVetoElectron"    : OC_veto_electron_indices,
            },
        },
        aux=selection_steps,
    ), good_electron_indices, ele_emu_indices, single_veto_electron_indices, second_ele_veto_emu_indices, OC_veto_electron_indices



# ------------------------------------------------------------------------------------------------------- #
# Tau Selection
# Reference:
#   https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
#   http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
# ------------------------------------------------------------------------------------------------------- #
@selector(
    uses={
        # Tau nano columns
        f"Tau.{var}" for var in [
            "pt", "eta", "phi", "dz",
            "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu", "idDeepTau2018v2p5VSjet",
            "decayMode","decayModePNet", 'ipLengthSig'
        ]
    },
    produces={
        f"Tau.{var}" for var in [
            "rawIdx","ip_sig"
        ]
    },
    exposed=False,
)
def tau_selection(
        self: Selector,
        events: ak.Array,
        electron_indices: Optional[ak.Array]=None,
        muon_indices    : Optional[ak.Array]=None,
        **kwargs
) -> tuple[ak.Array, SelectionResult, ak.Array]:
    """
    Tau selection returning two sets of indidces for default and veto muons.
    
    References:
      - 
    """
    tau_local_indices = ak.local_index(events.Tau)
    events = set_ak_column(events, "Tau.rawIdx", tau_local_indices)
    events = set_ak_column(events, "Tau.ip_sig", events.Tau.ipLengthSig)
   

    # https://cms-nanoaod-integration.web.cern.ch/integration/cms-swmaster/data106Xul17v2_v10_doc.html#Tau
    
    # Definition of the working points of the DeepTau 
    deep_tau_vs_e_jet_wps = self.config_inst.x.deep_tau.vs_e_jet_wps
    deep_tau_vs_mu_wps = self.config_inst.x.deep_tau.vs_mu_wps

    good_selections = {
        "tau_pt_30"     : events.Tau.pt > 30,
        "tau_eta_2p3"   : abs(events.Tau.eta) < 2.3,
        "tau_dz_0p2"    : abs(events.Tau.dz) < 0.2,
        "DeepTauVSjet"  : events.Tau.idDeepTau2018v2p5VSjet >= deep_tau_vs_e_jet_wps["VVVLoose"], 
        "DeepTauVSe"    : events.Tau.idDeepTau2018v2p5VSe   >= deep_tau_vs_e_jet_wps["VVVLoose"],
        "DeepTauVSmu"   : events.Tau.idDeepTau2018v2p5VSmu  >= deep_tau_vs_mu_wps["VLoose"],
        "DecayMode"     : ((events.Tau.decayMode == 0) 
                           | (events.Tau.decayMode == 1)
                           | (events.Tau.decayMode == 10)
                           | (events.Tau.decayMode == 11))
    }
    # pt sorted indices for converting masks to indices
    sorted_indices = ak.argsort(events.Tau.pt, axis=-1, ascending=False)
    tau_mask  = ak.ones_like(events.Tau.pt, dtype=np.bool_)
    
    for cut in good_selections.keys():
        tau_mask = tau_mask & good_selections[cut] 
        
    # convert to sorted indices
    selected_tau_idx = sorted_indices[tau_mask[sorted_indices]]
    selected_tau_idx = ak.values_astype(selected_tau_idx, np.int32)

    return events, selected_tau_idx


# ------------------------------------------------------------------------------------------------------- #
# Jet Selection
# Reference:
#   https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
#   http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
# ------------------------------------------------------------------------------------------------------- #
@selector(
    uses={ f"Jet.{var}" for var in 
        [
            "pt", "eta", "phi", "mass",
            "jetId", "btagDeepFlavB",
        ]} | {optional("Jet.puId")} ,
    exposed=False,
)
def jet_selection(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    """
    This function vetoes b-jets with sufficiently high pt and incide eta region of interest
    """
    year = self.config_inst.campaign.x.year
    is_run2 = (self.config_inst.campaign.x.year in [2016,2017,2018])
    
   
    # nominal selection
    good_selections = {
        "jet_pt_30"               : events.Jet.pt > 30.0,
        "jet_eta_2.4"             : abs(events.Jet.eta) < 2.4,
        #"jet_id"                  : events.Jet.jetId == 0b110,  # Jet ID flag: bit2 is tight, bit3 is tightLepVeto 
    }
    
    #if is_run2: 
    #    good_selections["jet_puId"] = ((events.Jet.pt >= 50.0) | (events.Jet.puId)) #For the Run2 there was
    
    # b-tagged jets, tight working point
    
    jet_mask  = ak.local_index(events.Jet.pt) >= 0 #Create a mask filled with ones
    selection_steps = {}

    selection_steps = {"Starts with": jet_mask}
    for cut in good_selections.keys():
        jet_mask = jet_mask & ak.fill_none(good_selections[cut], False)
        selection_steps[cut] = jet_mask
    
    sorted_indices = ak.argsort(events.Jet.pt, axis=-1, ascending=False)
    good_jet_indices = sorted_indices[jet_mask[sorted_indices]]
    good_jet_indices = ak.values_astype(good_jet_indices, np.int32)

    # b-tagged jets, tight working point
    btag_wp = self.config_inst.x.btag_working_points[year].deepjet.medium
    b_jet_mask = jet_mask & (events.Jet.btagDeepFlavB >= btag_wp)
    selection_steps["btag"] = ak.fill_none(b_jet_mask, False)

    # bjet veto
    bjet_veto = ak.sum(b_jet_mask, axis=1) == 0

    return events, SelectionResult(
        steps = {
            "b_veto": bjet_veto,
        }, 
        objects = {
            "Jet": {
                "Jet": good_jet_indices,
            },
        },
        aux = selection_steps,
    )
    


# ------------------------------------------------------------------------------------------------------- #
# GenTau Selection
# Reference:
#   GenPart =  ['eta', 'mass', 'phi', 'pt', 'genPartIdxMother', 'pdgId', 'status', 'statusFlags', 
#               'genPartIdxMotherG', 'distinctParentIdxG', 'childrenIdxG', 'distinctChildrenIdxG', 
#               'distinctChildrenDeepIdxG']
# ------------------------------------------------------------------------------------------------------- #
@selector(
    uses={
        "GenPart.*",
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass",
    },
    produces={
        'GenPart.rawIdx',
        'GenTau.rawIdx', 'GenTau.eta', 'GenTau.mass', 'GenTau.phi', 'GenTau.pt', 'GenTau.pdgId', 'GenTau.decayMode', 'GenTau.charge',
        'GenTauProd.rawIdx', 'GenTauProd.eta', 'GenTauProd.mass', 'GenTauProd.phi', 'GenTauProd.pt', 'GenTauProd.pdgId',
    },
    mc_only=True,
    exposed=False,
)
def gentau_selection(
        self: Selector,
        events: ak.Array,
        match: bool,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    """
    Selecting the generator level taus only, no martching here
    select the gen tau decay products as well
    
    References:
      - 
    """
    genpart_indices = ak.local_index(events.GenPart.pt)
    events = set_ak_column(events, "GenPart.rawIdx", genpart_indices)

    # masks to select gen tau+ and tau-    
    good_selections = {
        "genpart_pdgId"           : np.abs(events.GenPart.pdgId) == 15,
        "genpart_status"          : events.GenPart.status == 2,
        "genpart_status_flags"    : events.GenPart.hasFlags(["isPrompt", "isFirstCopy"]),
        "genpart_pt_10"           : events.GenPart.pt > 10.0,
        "genpart_eta_2p3"         : np.abs(events.GenPart.eta) < 2.3,
        "genpart_momid_25"        : events.GenPart[events.GenPart.distinctParent.genPartIdxMother].pdgId == 25,
        "genpart_mom_status_22"   : events.GenPart[events.GenPart.distinctParent.genPartIdxMother].status == 22,
    }
    
    gen_mask  = genpart_indices >= 0
    good_gen_mask = gen_mask

    selection_steps = {"Starts with": good_gen_mask}
    for cut in good_selections.keys():
        good_gen_mask = good_gen_mask & ak.fill_none(good_selections[cut], False)
        selection_steps[cut] = good_gen_mask

    gentau_indices = genpart_indices[good_gen_mask]
    
    gentaus = ak.with_name(events.GenPart[gentau_indices], "PtEtaPhiMLorentzVector")
    hcands  = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")

    # check the matching
    matched_gentaus = gentaus
    if match:
        matched_gentaus = hcands.nearest(gentaus, threshold=0.5)
        
    # nearest method can include None if a particle is not matched
    # so, taking care of the none values before adding it as a new column
    is_none = ak.sum(ak.is_none(matched_gentaus, axis=1), axis=1) > 0
    matched_gentaus = ak.where(is_none, gentaus[:,:0], matched_gentaus)
    has_full_match           = ~is_none
    has_two_matched_gentaus  = has_full_match & ak.fill_none(ak.num(matched_gentaus.rawIdx, axis=1) == 2, False)
    gentaus_of_opposite_sign = ak.fill_none(ak.sum(matched_gentaus.pdgId, axis=1) == 0, False)

    # Get gentau decay products 
    # hack: _apply_global_index [todo: https://github.com/columnflow/columnflow/discussions/430]
    # get decay modes for the GenTaus
    decay_gentau_indices = matched_gentaus.distinctChildrenIdxG
    decay_gentaus = events.GenPart._apply_global_index(decay_gentau_indices)
    gentaus_dm = getGenTauDecayMode(decay_gentaus)

    mask_genmatchedtaus_1 = ak.fill_none(ak.firsts((( gentaus_dm[:,:1] == -2)
                                                    |(gentaus_dm[:,:1] == -1)
                                                    |(gentaus_dm[:,:1] == 0) 
                                                    |(gentaus_dm[:,:1] == 1)
                                                    |(gentaus_dm[:,:1] == 2)
                                                    |(gentaus_dm[:,:1] == 10)
                                                    |(gentaus_dm[:,:1] == 11)), axis=1), False) # ele/had
    mask_genmatchedtaus_2 = ak.fill_none(ak.firsts((( gentaus_dm[:,1:2] == 0) 
                                                    |(gentaus_dm[:,1:2] == 1)
                                                    |(gentaus_dm[:,1:2] == 2)
                                                    |(gentaus_dm[:,1:2] == 10)
                                                    |(gentaus_dm[:,1:2] == 11)), axis=1), False) # had only

    mask_genmatchedtaus   = mask_genmatchedtaus_1 & mask_genmatchedtaus_2
    

    # check decaymodes
    # make sure that the decay mode of hcand is the same as decay mode of GenTau
    dm_match_evt_mask = ak.num(hcands.decayMode, axis=1) == 2
    if match:
        has_2         = (ak.num(hcands.decayMode, axis=1) == 2) & (ak.num(gentaus_dm, axis=1) == 2)
        _gentaus_dm   = ak.where(has_2, gentaus_dm, gentaus_dm[:,:0])
        _hcands_dm    = ak.where(has_2, hcands.decayMode, hcands.decayMode[:,:0])
        dm_match_mask = _hcands_dm == _gentaus_dm
        dm_match_evt_mask = ak.sum(dm_match_mask, axis=1) == 2
    #1/0
    

    # creating a proper array to save it as a new column
    dummy_decay_gentaus = decay_gentaus[:,:0][:,None]
    decay_1             = decay_gentaus[:,:1]
    decay_1             = ak.where(ak.num(decay_1, axis=1) > 0, decay_1, dummy_decay_gentaus)
    decay_2             = decay_gentaus[:,1:2]
    decay_2             = ak.where(ak.num(decay_2, axis=1) > 0, decay_2, dummy_decay_gentaus)
    decay_gentaus       = ak.concatenate([decay_1, decay_2], axis=1)
    
    # WARNING: Not a smart way to convert ak.Array -> List -> ak.Array
    # Must use ak.enforce_type: NOT WORKING here, but it was good for hcand selection
    events = set_ak_column(events, "GenTau",           ak.Array(ak.to_list(matched_gentaus)))
    events = set_ak_column(events, "GenTau.decayMode", gentaus_dm)
    events = set_ak_column(events, "GenTau.mass",      ak.ones_like(events.GenTau.mass) * 1.777)
    events = set_ak_column(events, "GenTau.charge",    ak.where(events.GenTau.pdgId > 0,
                                                                -1, 
                                                                ak.where(events.GenTau.pdgId < 0, 
                                                                         1, 
                                                                         0)))
    events = set_ak_column(events, "GenTauProd",       ak.Array(ak.to_list(decay_gentaus)))

    #1/0

    return events, SelectionResult(
        steps = {
            "has two matched gentaus"  : has_two_matched_gentaus,
            "gentaus of opposite sign" : gentaus_of_opposite_sign,
            "valid decay products"     : mask_genmatchedtaus,
            "gen DMs same as hcands"   : dm_match_evt_mask,
        },
        aux = selection_steps,
    )
